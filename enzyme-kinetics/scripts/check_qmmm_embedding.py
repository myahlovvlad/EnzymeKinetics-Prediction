"""PySCF QM/point-charge electrostatic embedding smoke test; not enzyme QM/MM."""

import argparse
import json
from pathlib import Path

from pyscf import gto, qmmm, scf, __version__


def run():
    mol = gto.M(
        atom="O 0 0 0; H 0 .757 .587; H 0 -.757 .587",
        unit="Angstrom",
        basis="sto-3g",
        verbose=0,
    )
    vac = scf.RHF(mol)
    vacuum = float(vac.kernel())
    mm = qmmm.mm_charge(scf.RHF(mol), [[3.0, 0.0, 0.0]], [0.5], unit="Angstrom")
    embedded = float(mm.kernel())
    if not vac.converged or not mm.converged:
        raise RuntimeError("SCF did not converge")
    return {
        "pyscf_version": __version__,
        "qm_system": "water RHF/STO-3G",
        "mm_environment": "one fixed +0.5e point charge at (3,0,0) Å",
        "vacuum_energy_hartree": vacuum,
        "embedded_energy_hartree": embedded,
        "energy_difference_hartree": embedded - vacuum,
        "scf_converged": True,
        "limitations": [
            "electrostatic embedding smoke test only",
            "not a coupled PySCF/OpenMM enzyme dynamics driver",
            "no QM/MM boundary, PME, MM relaxation, reaction path or free energy",
        ],
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError("Choose a new output")
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
