"""Generate a GAFF2/AM1-BCC ligand system from a specific screened 5-AND pose.

This creates a ligand-only model, not a parameterized protein-ligand complex.
Requires an external conda environment with OpenFF, AmberTools and OpenMM.
"""

import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

SMILES = "C[C@]12CCC(=O)CC1=CC[C@@H]1[C@@H]2CC[C@]2(C)C(=O)CC[C@@H]12"


def read_xyz(path, expected_symbols):
    lines = path.read_text().splitlines()
    entries = [line.split() for line in lines[2:]]
    if int(lines[0]) != len(entries) or len(entries) != len(expected_symbols):
        raise ValueError("XYZ atom count differs from ligand topology")
    if [r[0] for r in entries] != expected_symbols:
        raise ValueError("XYZ element order differs from ligand topology")
    return [[float(v) for v in row[1:]] for row in entries]


def main():
    from openff.toolkit import Molecule
    from openff.toolkit.utils import AmberToolsToolkitWrapper
    from openff.units import unit as offunit
    import openmm
    from openmm import app, unit
    from openmmforcefields.generators import GAFFTemplateGenerator

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--template", type=Path, required=True)
    p.add_argument("--xyz", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError("Choose a new output directory")
    lines = args.template.read_text().splitlines()
    block = (
        "\n".join(
            x
            for x in lines
            if x.startswith("HETATM") and x[17:20] == "ASD" and x[21] == "A"
        )
        + "\nEND\n"
    )
    mol = AllChem.AssignBondOrdersFromTemplate(
        Chem.MolFromSmiles(SMILES), Chem.MolFromPDBBlock(block, removeHs=True)
    )
    Chem.AssignStereochemistryFrom3D(mol)
    if Chem.MolToSmiles(mol) != Chem.MolToSmiles(Chem.MolFromSmiles(SMILES)):
        raise ValueError("Template ligand stereo mismatch")
    mol = Chem.AddHs(mol, addCoords=True)
    positions = read_xyz(args.xyz, [a.GetSymbol() for a in mol.GetAtoms()])
    for i, pos in enumerate(positions):
        mol.GetConformer().SetAtomPosition(i, pos)
    offmol = Molecule.from_rdkit(
        mol, allow_undefined_stereo=False, hydrogens_are_explicit=True
    )
    offmol.assign_partial_charges(
        "am1bcc",
        use_conformers=[offmol.conformers[0]],
        toolkit_registry=AmberToolsToolkitWrapper(),
    )
    ff = app.ForceField()
    generator = GAFFTemplateGenerator(molecules=[offmol], forcefield="gaff-2.11")
    ff.registerTemplateGenerator(generator.generator)
    topology = offmol.to_topology().to_openmm()
    system = ff.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=None)
    if system.getNumParticles() != mol.GetNumAtoms():
        raise ValueError("Force field particle count mismatch")
    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    context = openmm.Context(
        system, integrator, openmm.Platform.getPlatformByName("CPU")
    )
    context.setPositions(positions * unit.angstrom)
    before = (
        context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )
    openmm.LocalEnergyMinimizer.minimize(context, maxIterations=200)
    state = context.getState(getEnergy=True, getPositions=True)
    after = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "ligand.sdf").open("w") as f:
        f.write(Chem.MolToMolBlock(mol) + "\n$$$$\n")
    with (args.output / "ligand-minimized.pdb").open("w") as f:
        app.PDBFile.writeFile(topology, state.getPositions(), f, keepIds=False)
    (args.output / "ligand-system.xml").write_text(
        openmm.XmlSerializer.serialize(system)
    )
    report = {
        "input_sha256": {
            str(x): hashlib.sha256(x.read_bytes()).hexdigest()
            for x in (args.template, args.xyz)
        },
        "rdkit": rdBase.rdkitVersion,
        "openmm": openmm.__version__,
        "forcefield": "gaff-2.11",
        "charge_method": "AM1-BCC via AmberTools",
        "total_partial_charge_e": float(
            offmol.partial_charges.m_as(offunit.elementary_charge).sum()
        ),
        "atoms": mol.GetNumAtoms(),
        "potential_kj_mol_before": before,
        "potential_kj_mol_after": after,
        "scope": "isolated ligand, gas-phase NoCutoff; not binding energy or complex minimization",
        "limitations": [
            "Conformation depends on unvalidated mutant-product template",
            "No protein, crystal waters, solvent or intermolecular interaction included.",
        ],
    }
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
