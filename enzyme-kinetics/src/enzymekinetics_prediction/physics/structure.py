"""P2Rank execution and normalized, provenance-preserving pocket features."""

from collections import Counter
import os
from pathlib import Path
import re
import subprocess

import numpy as np
import pandas as pd

from .artifacts import artifact, digest_file

CATALYTIC = ("HIS", "ASP", "GLU", "SER", "CYS", "TYR", "LYS", "ARG")


def protein_id(path):
    name = Path(path).name.replace("_predictions.csv", "")
    name = re.sub(r"\.(pdb|cif)$", "", name)
    match = re.fullmatch(
        r"(?:AF-|ESMFold-)?([A-Z0-9]+)(?:-F\d+-model_v\d+|-v\d+)?", name
    )
    if not match:
        raise ValueError(f"Cannot normalize protein ID: {name}")
    return match.group(1)


def ca_residues(pdb):
    residues = {}
    for line in Path(pdb).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if (
            line.startswith("ATOM  ")
            and line[12:16].strip() == "CA"
            and line[16] in " A"
        ):
            key = f"{line[21]}_{line[22:26].strip()}{line[26].strip()}"
            residues.setdefault(
                key,
                (
                    line[17:20].strip(),
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                ),
            )
    if not residues:
        raise ValueError("PDB contains no CA residues")
    return residues


def run_p2rank(pdb, output, vendor, threads=1, timeout=180, config=None):
    if threads < 1 or timeout <= 0:
        raise ValueError("Invalid process limits")
    pdb, output, vendor = (
        Path(pdb).resolve(),
        Path(output).resolve(),
        Path(vendor).resolve(),
    )
    executable = vendor / "prank"
    if not pdb.is_file() or not executable.is_file():
        raise FileNotFoundError("PDB input or P2Rank distribution not found")
    output.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    if not environment.get("JAVA_HOME"):
        try:
            import jdk4py

            environment["JAVA_HOME"] = str(jdk4py.JAVA_HOME)
        except ImportError:
            pass
    command = [
        str(executable),
        "predict",
        "-f",
        str(pdb),
        "-o",
        str(output),
        "-threads",
        str(threads),
    ]
    if config:
        command.extend(["-c", config])
    completed = subprocess.run(
        command,
        env=environment,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )
    (output / "p2rank.stdout.log").write_text(completed.stdout)
    (output / "p2rank.stderr.log").write_text(completed.stderr)
    if completed.returncode:
        raise RuntimeError(f"P2Rank failed ({completed.returncode}); see {output}")
    expected = output / f"{pdb.name}_predictions.csv"
    if not expected.is_file():
        raise RuntimeError("P2Rank returned no prediction CSV")
    return expected


def pocket_features(predictions, pdb=None):
    predictions = Path(predictions)
    frame = pd.read_csv(predictions, skip_blank_lines=True)
    frame.columns = frame.columns.str.strip()
    if frame.empty:
        raise ValueError("No binding pockets predicted")
    required = {"rank", "score", "residue_ids"}
    if not required <= set(frame):
        raise ValueError(f"Missing P2Rank columns: {required - set(frame)}")
    top = frame.sort_values("rank").iloc[0]
    keys = str(top.residue_ids).split() if pd.notna(top.residue_ids) else []
    if not keys:
        raise ValueError("Top pocket has no residue IDs")
    values = {"score": float(top.score), "residue_count": len(set(keys))}
    units = {"score": "P2Rank_score", "residue_count": "count"}
    inputs = {
        "predictions": str(predictions),
        "predictions_sha256": digest_file(predictions),
        "protein_id": protein_id(predictions),
    }
    names = []
    matched = 0
    limitations = [
        "Predicted binding pocket is not automatically a verified catalytic site",
        "Charge proxy counts residue types and is not an electrostatics calculation",
    ]
    if pdb is not None:
        pdb = Path(pdb)
        residues = ca_residues(pdb)
        selected = [residues[k] for k in dict.fromkeys(keys) if k in residues]
        names = [name for name, _ in selected]
        matched = len(selected)
        if selected:
            positions = np.asarray([position for _, position in selected])
            values["radius_gyration_A"] = float(
                np.sqrt(
                    np.mean(np.sum((positions - positions.mean(axis=0)) ** 2, axis=1))
                )
            )
            units["radius_gyration_A"] = "Angstrom"
        inputs.update(pdb=str(pdb), pdb_sha256=digest_file(pdb))
    counts = Counter(names)
    for name in CATALYTIC:
        values[f"freq_{name}"] = counts[name] / len(names) if names else None
        units[f"freq_{name}"] = "fraction"
    values["charge_proxy"] = (
        (
            sum(counts[x] for x in ("LYS", "ARG", "HIS"))
            - sum(counts[x] for x in ("ASP", "GLU"))
        )
        / len(names)
        if names
        else None
    )
    values["hydrophobic_fraction"] = (
        sum(counts[x] for x in ("ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"))
        / len(names)
        if names
        else None
    )
    values["residue_mapping_fraction"] = matched / len(set(keys))
    units.update(
        charge_proxy="residue_type_proxy",
        hydrophobic_fraction="fraction",
        residue_mapping_fraction="fraction",
    )
    if matched != len(set(keys)):
        limitations.append("Some or all pocket residues lack a matching PDB CA record")
    if all(f"center_{axis}" in frame for axis in "xyz"):
        values["center_A"] = [float(top[f"center_{axis}"]) for axis in "xyz"]
        units["center_A"] = "Angstrom"
    return artifact(
        "P2Rank",
        inputs,
        {"selected_rank": int(top["rank"])},
        values,
        units,
        limitations=limitations,
    )
