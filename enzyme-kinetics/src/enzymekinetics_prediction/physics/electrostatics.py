"""pH-aware PQR preparation, explicit screened Coulomb, and APBS/OpenDX."""

from pathlib import Path
import os
import shutil
import subprocess

import numpy as np
from scipy.constants import Avogadro, Boltzmann, elementary_charge, epsilon_0
from scipy.interpolate import RegularGridInterpolator

from .artifacts import artifact, digest_file, package_version
from .kinetics import positive
from enzymekinetics_prediction.paths import VENDOR_DIR


def prepare_pqr(pdb, output, ph, timeout=180):
    if not np.isfinite(ph) or not 0 <= ph <= 14:
        raise ValueError("pH must be within [0,14]")
    executable = shutil.which("pdb2pqr")
    if executable is None:
        raise RuntimeError(
            "pdb2pqr executable not installed; install the physics extra"
        )
    pdb, output = Path(pdb).resolve(), Path(output).resolve()
    if not pdb.is_file():
        raise FileNotFoundError(pdb)
    output.parent.mkdir(parents=True, exist_ok=True)
    command = [
        executable,
        "--ff=AMBER",
        "--titration-state-method=propka",
        f"--with-ph={ph}",
        "--drop-water",
        "--keep-chain",
        "--whitespace",
        str(pdb),
        str(output),
    ]
    completed = subprocess.run(
        command, capture_output=True, text=True, timeout=timeout, check=False
    )
    output.with_suffix(".stdout.log").write_text(completed.stdout)
    output.with_suffix(".stderr.log").write_text(completed.stderr)
    if completed.returncode or not output.is_file():
        raise RuntimeError(f"PDB2PQR failed: see {output.with_suffix('.stderr.log')}")
    return output


def read_pqr(path):
    positions, charges, radii = [], [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            fields = line.split()
            if len(fields) < 10:
                raise ValueError("Malformed PQR atom record")
            x, y, z, charge, radius = map(float, fields[-5:])
            positions.append([x, y, z])
            charges.append(charge)
            radii.append(radius)
    if not positions:
        raise ValueError("PQR contains no atoms")
    xyz, q, r = np.asarray(positions), np.asarray(charges), np.asarray(radii)
    if not (
        np.isfinite(xyz).all() and np.isfinite(q).all() and np.isfinite(r).all()
    ) or np.any(r < 0):
        raise ValueError("Invalid PQR coordinates, charges or radii")
    return xyz, q, r


def screened_coulomb(
    positions_A,
    charges_e,
    query_A,
    *,
    temperature_K=298.15,
    ionic_strength_M=0.15,
    dielectric=78.5,
):
    """Homogeneous Debye-Hückel approximation; not APBS and not a dielectric boundary solver."""
    positive(temperature_K, "temperature_K")
    positive(dielectric, "dielectric")
    positive(ionic_strength_M, "ionic_strength_M", allow_zero=True)
    xyz, charges, points = (
        np.asarray(positions_A, float),
        np.asarray(charges_e, float),
        np.atleast_2d(query_A).astype(float),
    )
    if (
        xyz.ndim != 2
        or xyz.shape[1] != 3
        or charges.shape != (len(xyz),)
        or points.shape[1] != 3
    ):
        raise ValueError("Inconsistent electrostatics arrays")
    if not all(np.isfinite(a).all() for a in (xyz, charges, points)):
        raise ValueError("Nonfinite electrostatics inputs")
    kappa = np.sqrt(
        2
        * Avogadro
        * ionic_strength_M
        * 1000
        * elementary_charge**2
        / (epsilon_0 * dielectric * Boltzmann * temperature_K)
    )
    potential, fields = [], []
    prefactor = elementary_charge / (4 * np.pi * epsilon_0 * dielectric)
    for point in points:
        delta = (point - xyz) * 1e-10
        distance = np.linalg.norm(delta, axis=1)
        if np.any(distance < 1e-12):
            raise ValueError(
                "Query coincides with a source charge; self fields must be excluded explicitly"
            )
        screen = np.exp(-kappa * distance)
        potential.append(float(prefactor * np.sum(charges * screen / distance)))
        fields.append(
            (
                prefactor
                * np.sum(
                    (charges * screen * (1 + kappa * distance) / distance**3)[:, None]
                    * delta,
                    axis=0,
                )
                / 1e10
            ).tolist()
        )
    return np.asarray(potential), np.asarray(fields)


def pqr_observables(
    path, query_A, *, ph, temperature_K=298.15, ionic_strength_M=0.15, dielectric=78.5
):
    xyz, charges, radii = read_pqr(path)
    potential, field = screened_coulomb(
        xyz,
        charges,
        query_A,
        temperature_K=temperature_K,
        ionic_strength_M=ionic_strength_M,
        dielectric=dielectric,
    )
    return artifact(
        "screened-coulomb",
        {"pqr": str(path), "pqr_sha256": digest_file(path)},
        {
            "pH": ph,
            "temperature_K": temperature_K,
            "ionic_strength_M": ionic_strength_M,
            "dielectric": dielectric,
            "query_A": np.atleast_2d(query_A).tolist(),
            "pdb2pqr_version": package_version("pdb2pqr"),
        },
        {
            "potential_V": potential.tolist(),
            "field_V_per_A": field.tolist(),
            "net_charge_e": float(charges.sum()),
        },
        {
            "potential_V": "Volt",
            "field_V_per_A": "Volt/Angstrom",
            "net_charge_e": "elementary_charge",
        },
        limitations=[
            "Homogeneous screened Coulomb approximation; not APBS/Poisson-Boltzmann",
            "Protein PQR does not include unknown ligand charges or a reaction transition state",
        ],
    )


def read_dx(path, temperature_K):
    positive(temperature_K, "temperature_K")
    header, marker, numerical = Path(path).read_text().partition("data follows")
    if not marker:
        raise ValueError("Incomplete OpenDX grid")
    counts, origin, deltas = None, None, []
    for line in header.splitlines():
        tokens = line.split()
        if "gridpositions" in tokens and "counts" in tokens:
            counts = tuple(map(int, tokens[-3:]))
        elif tokens[:1] == ["origin"]:
            origin = np.asarray(tokens[1:4], float)
        elif tokens[:1] == ["delta"]:
            deltas.append(list(map(float, tokens[1:4])))
    if (
        counts is None
        or origin is None
        or len(deltas) != 3
        or any(c < 2 for c in counts)
    ):
        raise ValueError("Incomplete OpenDX grid")
    expected = int(np.prod(counts))
    if expected > 100_000_000:
        raise ValueError("OpenDX grid exceeds this reader's resource limit")
    try:
        values = np.fromstring(numerical, dtype=float, sep=" ", count=expected)
    except ValueError:
        raise ValueError("Malformed OpenDX numerical data") from None
    if len(values) != expected:
        raise ValueError("Incomplete OpenDX grid")
    delta = np.asarray(deltas)
    if not np.allclose(delta, np.diag(np.diag(delta))) or np.any(np.diag(delta) <= 0):
        raise ValueError("Only axis-aligned increasing OpenDX grids are supported")
    potential_V = (
        np.asarray(values).reshape(counts)
        * Boltzmann
        * temperature_K
        / elementary_charge
    )
    if not np.isfinite(potential_V).all():
        raise ValueError("Nonfinite OpenDX grid")
    axes = tuple(origin[i] + np.arange(counts[i]) * delta[i, i] for i in range(3))
    field = [-g for g in np.gradient(potential_V, *np.diag(delta), edge_order=1)]
    return axes, potential_V, field


def sample_dx(path, query_A, temperature_K):
    axes, potential, field = read_dx(path, temperature_K)
    points = np.atleast_2d(query_A)
    phi = RegularGridInterpolator(axes, potential, bounds_error=True)(points)
    vector = np.column_stack(
        [RegularGridInterpolator(axes, a, bounds_error=True)(points) for a in field]
    )
    return phi, vector


def run_apbs(
    pqr,
    output,
    query_A,
    *,
    temperature_K=298.15,
    ionic_strength_M=0.15,
    dielectric=78.5,
    timeout=180,
    threads=1,
    grid_points=65,
    focus_length_A=None,
):
    positive(temperature_K, "temperature_K")
    positive(ionic_strength_M, "ionic_strength_M", allow_zero=True)
    positive(dielectric, "dielectric")
    if (
        not isinstance(grid_points, int)
        or not 33 <= grid_points <= 257
        or (grid_points - 1) % 32
    ):
        raise ValueError("APBS grid_points must be 32*n+1 within [33,257]")
    if not isinstance(threads, int) or threads < 1:
        raise ValueError("Positive APBS thread count required")
    executable = os.environ.get("APBS_EXECUTABLE") or shutil.which("apbs")
    if executable is None:
        local = VENDOR_DIR / "apbs" / "APBS-3.5.0.Linux" / "bin" / "apbs"
        if local.is_file() and os.access(local, os.X_OK):
            executable = str(local)
    if executable is None:
        raise RuntimeError(
            "APBS executable unavailable; explicitly choose screened-coulomb if that approximation is intended"
        )
    xyz, _, _ = read_pqr(pqr)
    query = np.atleast_2d(query_A).astype(float)
    if query.shape[1] != 3 or not np.isfinite(query).all():
        raise ValueError("Expected finite query coordinates in Angstrom")
    combined = np.vstack([xyz, query])
    center = (combined.max(axis=0) + combined.min(axis=0)) / 2
    fine = np.maximum(np.ptp(combined, axis=0), 20) + 20
    coarse = fine * 1.7
    fine_center = center.copy()
    if focus_length_A is not None:
        positive(focus_length_A, "focus_length_A")
        if focus_length_A < 10:
            raise ValueError("Focused box must be at least 10 Angstrom wide")
        fine_center = (query.max(axis=0) + query.min(axis=0)) / 2
        fine = np.full(3, focus_length_A, dtype=float)
        if np.any(np.max(np.abs(query - fine_center), axis=0) >= fine / 2 - 1):
            raise ValueError(
                "Query points must lie inside the focused box with a margin"
            )
        coarse = np.maximum(coarse, 2 * np.abs(fine_center - center) + fine + 20)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(pqr, output / "input.pqr")
    text = f"""read
 mol pqr input.pqr
end
elec
 mg-auto
 dime {grid_points} {grid_points} {grid_points}
 cglen {" ".join(map(str, coarse))}
 fglen {" ".join(map(str, fine))}
 cgcent {" ".join(map(str, center))}
 fgcent {" ".join(map(str, fine_center))}
 mol 1
 lpbe
 bcfl sdh
 pdie 2.0
 sdie {dielectric}
 srfm smol
 chgm spl2
 sdens 10.0
 srad 1.4
 swin 0.3
 temp {temperature_K}
 ion charge 1 conc {ionic_strength_M} radius 2.0
 ion charge -1 conc {ionic_strength_M} radius 2.0
 calcenergy no
 calcforce no
 write pot dx potential
end
quit
"""
    (output / "apbs.in").write_text(text)
    result = subprocess.run(
        [executable, "apbs.in"],
        cwd=output,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
        env={
            **os.environ,
            "OMP_NUM_THREADS": str(threads),
            "OPENBLAS_NUM_THREADS": str(threads),
        },
    )
    (output / "apbs.log").write_text(result.stdout + "\n" + result.stderr)
    if result.returncode:
        raise RuntimeError(f"APBS failed; see {output / 'apbs.log'}")
    dx = output / "potential.dx"
    phi, field = sample_dx(dx, query, temperature_K)
    return artifact(
        "APBS",
        {
            "pqr_sha256": digest_file(pqr),
            "dx_sha256": digest_file(dx),
            "executable_sha256": digest_file(shutil.which(executable) or executable),
        },
        {
            "temperature_K": temperature_K,
            "ionic_strength_M": ionic_strength_M,
            "grid": [grid_points] * 3,
            "threads": threads,
            "coarse_lengths_A": coarse.tolist(),
            "fine_lengths_A": fine.tolist(),
            "coarse_center_A": center.tolist(),
            "fine_center_A": fine_center.tolist(),
            "dielectric_solvent": dielectric,
            "dielectric_protein": 2.0,
            "query_A": query.tolist(),
        },
        {"potential_V": phi.tolist(), "field_V_per_A": field.tolist()},
        {"potential_V": "Volt", "field_V_per_A": "Volt/Angstrom"},
        limitations=[
            "Linear Poisson-Boltzmann with specified dielectric boundary and ion assumptions; grid convergence requires separate checks"
        ],
    )
