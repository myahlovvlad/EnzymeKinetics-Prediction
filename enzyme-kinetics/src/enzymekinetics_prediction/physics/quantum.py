"""GFN2-xTB molecular descriptors via tblite, with explicit geometry and units."""

from collections import Counter
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from scipy.constants import Boltzmann, physical_constants
from threadpoolctl import threadpool_limits

from .artifacts import artifact, content_key, package_version, write_json

HARTREE_J = physical_constants["Hartree energy"][0]
HARTREE_EV = physical_constants["Hartree energy in eV"][0]
BOHR_A = physical_constants["Bohr radius"][0] * 1e10
DEBYE_C_M = 3.33564e-30
E_CHARGE = physical_constants["elementary charge"][0]
AVOGADRO = physical_constants["Avogadro constant"][0]


@dataclass(frozen=True)
class QuantumSettings:
    method: str = "GFN2-xTB"
    seed: int = 42
    max_iterations: int = 250
    electronic_temperature_K: float = 300.0
    solvent: str | None = None
    fukui: bool = True
    threads: int = 1

    def __post_init__(self):
        if self.method not in {"GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"}:
            raise ValueError("Unsupported xTB Hamiltonian")
        if (
            self.threads < 1
            or self.max_iterations < 1
            or not np.isfinite(self.electronic_temperature_K)
            or self.electronic_temperature_K < 0
        ):
            raise ValueError("Invalid QM iteration count or electronic temperature")


def prepare_molecule(smiles, seed=42):
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles)
    if mol is None or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError("QM requires a valid non-generic molecular SMILES")
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError(
            "Compute separate molecular species explicitly; disconnected reactants are not a geometry"
        )
    canonical = Chem.MolToSmiles(mol)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError("Conformer generation failed")
    geometry_method = "ETKDGv3"
    if AllChem.MMFFHasAllMoleculeParams(mol):
        status = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        if status != 0:
            raise ValueError("MMFF geometry optimization did not converge")
        geometry_method += "+MMFF94"
    numbers = np.asarray([a.GetAtomicNum() for a in mol.GetAtoms()], dtype=np.int32)
    positions = np.asarray(mol.GetConformer().GetPositions(), dtype=np.float64)
    charge = int(Chem.GetFormalCharge(mol))
    unpaired = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    electrons = int(numbers.sum()) - charge
    if unpaired > electrons or (electrons - unpaired) % 2:
        raise ValueError("Charge and spin parity are inconsistent")
    return canonical, numbers, positions, charge, unpaired, geometry_method


def singlepoint(numbers, positions_A, charge, unpaired, settings):
    try:
        from tblite.interface import Calculator
    except ImportError as exc:
        raise RuntimeError("Install the physics extra for tblite calculations") from exc
    calculator = Calculator(
        settings.method,
        numbers,
        positions_A / BOHR_A,
        charge=float(charge),
        uhf=int(unpaired),
    )
    calculator.set("verbosity", 0)
    calculator.set("max-iter", settings.max_iterations)
    calculator.set(
        "temperature", Boltzmann * settings.electronic_temperature_K / HARTREE_J
    )
    if settings.solvent:
        calculator.add("alpb-solvation", settings.solvent)
    with threadpool_limits(limits=settings.threads):
        result = calculator.singlepoint()
    values = {
        name: np.asarray(result.get(name))
        for name in [
            "energy",
            "charges",
            "dipole",
            "orbital-energies",
            "orbital-occupations",
        ]
    }
    if not all(np.isfinite(value).all() for value in values.values()):
        raise ValueError("Nonfinite QM result")
    if not np.isclose(values["charges"].sum(), charge, atol=1e-4):
        raise ValueError("Atomic charges do not sum to molecular charge")
    return values


def molecular_descriptors(smiles, settings=QuantumSettings(), cache_dir=None):
    canonical, numbers, positions, charge, unpaired, geometry = prepare_molecule(
        smiles, settings.seed
    )
    inputs = {
        "smiles": canonical,
        "charge_e": charge,
        "unpaired_electrons": unpaired,
        "atomic_numbers": numbers.tolist(),
        "positions_A": positions.tolist(),
    }
    settings_dict = {
        **asdict(settings),
        "geometry_method": geometry,
        "tblite_version": package_version("tblite"),
        "rdkit_version": package_version("rdkit"),
    }
    key = content_key(
        {"schema": 1, "backend": "tblite", "inputs": inputs, "settings": settings_dict}
    )
    path = Path(cache_dir) / f"{key}.json" if cache_dir else None
    if path and path.is_file():
        cached = json.loads(path.read_text())
        if cached.get("cache_key") == key and cached.get("status") == "complete":
            return cached
    state = singlepoint(numbers, positions, charge, unpaired, settings)
    energy = state["orbital-energies"]
    occupations = state["orbital-occupations"]
    homo, lumo = None, None
    # Spin-summed molecular frontier orbitals are only reported for closed shells.
    if unpaired == 0 and energy.ndim == 1:
        occupied, vacant = energy[occupations > 1.0], energy[occupations < 1.0]
        if len(occupied) and len(vacant):
            homo, lumo = (
                float(occupied.max() * HARTREE_EV),
                float(vacant.min() * HARTREE_EV),
            )
    dipole = state["dipole"] * E_CHARGE * BOHR_A * 1e-10 / DEBYE_C_M
    values = {
        "energy_hartree": float(state["energy"]),
        "homo_eV": homo,
        "lumo_eV": lumo,
        "gap_eV": lumo - homo if homo is not None else None,
        "atomic_charges_e": state["charges"].tolist(),
        "dipole_D": dipole.tolist(),
        "dipole_norm_D": float(np.linalg.norm(dipole)),
    }
    units = {
        "energy_hartree": "Hartree",
        "homo_eV": "eV",
        "lumo_eV": "eV",
        "gap_eV": "eV",
        "atomic_charges_e": "elementary_charge",
        "dipole_D": "Debye",
        "dipole_norm_D": "Debye",
    }
    limitations = [
        "Single-point electronic descriptors at an ETKDG/MMFF geometry; not a TS or activation free energy",
        "Charged-molecule dipole depends on origin; coordinates record the chosen origin",
        "Open-shell frontier gap is intentionally not reported",
    ]
    if settings.fukui:
        plus_charge = charge - 1  # N+1 electrons
        minus_charge = charge + 1  # N-1 electrons
        if int(numbers.sum()) - minus_charge < 1:
            raise ValueError(
                "Fukui calculation requires at least one electron in all charge states"
            )
        added = singlepoint(
            numbers,
            positions,
            plus_charge,
            (int(numbers.sum()) - plus_charge) % 2,
            settings,
        )
        removed = singlepoint(
            numbers,
            positions,
            minus_charge,
            (int(numbers.sum()) - minus_charge) % 2,
            settings,
        )
        values["fukui_plus"] = (state["charges"] - added["charges"]).tolist()
        values["fukui_minus"] = (removed["charges"] - state["charges"]).tolist()
        units.update(fukui_plus="dimensionless", fukui_minus="dimensionless")
        limitations.append(
            "Vertical finite-difference Fukui; adjacent states use minimum spin consistent with electron parity"
        )
    output = artifact(
        "tblite", inputs, settings_dict, values, units, limitations=limitations
    )
    if path:
        write_json(path, output)
    return output


def electronic_dissociation(
    parent_smiles, fragment_smiles, settings=QuantumSettings(fukui=False)
):
    if len(fragment_smiles) != 2:
        raise ValueError("Specify exactly two explicit dissociation fragments")
    molecules = [Chem.MolFromSmiles(s) for s in [parent_smiles, *fragment_smiles]]
    if any(m is None for m in molecules):
        raise ValueError("Invalid dissociation species")
    counts = [
        Counter(a.GetAtomicNum() for a in Chem.AddHs(m).GetAtoms()) for m in molecules
    ]
    if counts[0] != counts[1] + counts[2] or Chem.GetFormalCharge(molecules[0]) != sum(
        Chem.GetFormalCharge(m) for m in molecules[1:]
    ):
        raise ValueError("Dissociation must conserve all atoms and charge")
    states = [
        molecular_descriptors(s, settings) for s in [parent_smiles, *fragment_smiles]
    ]
    delta = (
        sum(s["values"]["energy_hartree"] for s in states[1:])
        - states[0]["values"]["energy_hartree"]
    )
    return {
        "electronic_dissociation_kj_mol": delta * HARTREE_J * AVOGADRO / 1000,
        "species": states,
        "limitation": "Electronic energy difference of prepared species; not thermal BDE or activation barrier",
    }
