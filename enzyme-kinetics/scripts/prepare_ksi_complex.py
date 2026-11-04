"""Exploratory explicit-water WT/5-AND OpenMM minimization, not binding affinity."""
import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem


def select_ash_variants(topology):
    """Specify protonated Asp99 in both WT chains; other residues use OpenMM defaults."""
    residues = list(topology.residues())
    selected = [i for i, r in enumerate(residues) if r.name == "ASP" and r.id == "99" and r.chain.id in ("A", "B")]
    if len(selected) != 2:
        raise ValueError("Expected WT Asp99 in chains A and B")
    return ["ASH" if i in selected else None for i in range(len(residues))]


def main():
    import openmm
    from openmm import app, unit
    from pdbfixer import PDBFixer
    from openff.toolkit import Molecule
    from openff.toolkit.utils import AmberToolsToolkitWrapper
    from openmmforcefields.generators import GAFFTemplateGenerator

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--protein", type=Path, required=True)
    p.add_argument("--ligand-sdf", type=Path, required=True)
    p.add_argument("--crystal-assembly", type=Path, help="original 8CHO assembly, same coordinate frame")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--max-iterations", type=int, default=100)
    p.add_argument("--padding-nm", type=float, default=0.8)
    args = p.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError("Choose a new output directory")
    args.output.mkdir(parents=True, exist_ok=True)
    mol = next((m for m in Chem.SDMolSupplier(str(args.ligand_sdf), removeHs=False) if m), None)
    if mol is None:
        raise ValueError("Cannot parse ligand SDF")
    offmol = Molecule.from_rdkit(mol, hydrogens_are_explicit=True, allow_undefined_stereo=False)
    offmol.assign_partial_charges("am1bcc", use_conformers=[offmol.conformers[0]], toolkit_registry=AmberToolsToolkitWrapper())
    fixer = PDBFixer(filename=str(args.protein))
    fixer.findMissingResidues()
    missing_residues = {str(k): list(v) for k, v in fixer.missingResidues.items()}
    fixer.missingResidues = {}  # do not invent unresolved residues
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        raise ValueError("Nonstandard protein residues require manual review")
    fixer.findMissingAtoms()
    missing_atoms = sum(len(x) for x in fixer.missingAtoms.values())
    fixer.addMissingAtoms()
    ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
    generator = GAFFTemplateGenerator(molecules=[offmol], forcefield="gaff-2.11")
    ff.registerTemplateGenerator(generator.generator)
    modeller = app.Modeller(fixer.topology, fixer.positions)
    variants = select_ash_variants(modeller.topology)
    chosen = modeller.addHydrogens(ff, pH=7.0, variants=variants)
    modeller.add(offmol.to_topology().to_openmm(), mol.GetConformer().GetPositions() * unit.angstrom)
    crystal_water_count = 0
    if args.crystal_assembly:
        from Bio.PDB import PDBParser
        import numpy as np
        parsed = PDBParser(QUIET=True).get_structure("crystal", args.crystal_assembly)
        protein = PDBParser(QUIET=True).get_structure("protein", args.protein)[0]
        if len(parsed) != 2:
            raise ValueError("Expected two models in crystal assembly")
        for model_id, chain_id in ((0, "A"), (1, "B")):
            for resid in range(1, 126):
                if np.linalg.norm(parsed[model_id]["A"][resid]["CA"].coord - protein[chain_id][resid]["CA"].coord) > 0.01:
                    raise ValueError("Crystal and WT coordinate frames differ")
        oxygens = [a.coord.astype(float) for model in parsed for a in model.get_atoms()
                   if a.get_parent().id[0] == "W" and a.element == "O"]
        if not oxygens:
            raise ValueError("Crystal assembly has no parseable waters")
        ligand_xyz = mol.GetConformer().GetPositions()
        if min(np.linalg.norm(ligand_xyz[:21] - o, axis=1).min() for o in oxygens) < 2.4:
            raise ValueError("Crystal water overlaps ligand within 2.4 Å")
        water_top = app.Topology()
        chain = water_top.addChain("D")
        water_pos = []
        for oxygen in oxygens:
            residue = water_top.addResidue("HOH", chain)
            o = water_top.addAtom("O", app.element.oxygen, residue)
            h1 = water_top.addAtom("H1", app.element.hydrogen, residue)
            h2 = water_top.addAtom("H2", app.element.hydrogen, residue)
            water_top.addBond(o, h1)
            water_top.addBond(o, h2)
            # Initial water orientation is arbitrary, not measured in the crystal.
            water_pos.extend([oxygen, oxygen + [0.9572, 0, 0], oxygen + [-0.239987, 0.92663, 0]])
        modeller.add(water_top, water_pos * unit.angstrom)
        crystal_water_count = len(oxygens)
    modeller.addSolvent(ff, model="tip3p", padding=args.padding_nm * unit.nanometer, neutralize=True)
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=0.8 * unit.nanometer, constraints=app.HBonds)
    integ = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    context = openmm.Context(system, integ, openmm.Platform.getPlatformByName("CPU"))
    context.setPositions(modeller.positions)
    before = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    openmm.LocalEnergyMinimizer.minimize(context, maxIterations=args.max_iterations)
    state = context.getState(getPositions=True, getEnergy=True, getForces=True)
    after = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    with (args.output / "minimized.pdb").open("w") as f:
        app.PDBFile.writeFile(modeller.topology, state.getPositions(), f, keepIds=False)
    (args.output / "system.xml").write_text(openmm.XmlSerializer.serialize(system))
    inputs = [args.protein, args.ligand_sdf] + ([args.crystal_assembly] if args.crystal_assembly else [])
    report = {"input_sha256": {str(x): hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs},
              "protein_forcefield": "amber14-all.xml", "ligand_forcefield": "gaff-2.11", "charges": "AM1-BCC AmberTools",
              "water": "TIP3P with initial crystal oxygen positions" if crystal_water_count else "TIP3P added by OpenMM, not preserved crystallographic waters", "crystal_water_oxygen_count": crystal_water_count, "pH": 7.0,
              "ash99_forced_both_chains": True, "chosen_protonation_variants": chosen,
              "missing_residues_not_added": missing_residues, "missing_heavy_atoms_added": missing_atoms,
              "box_padding_nm": args.padding_nm, "atoms": system.getNumParticles(), "max_iterations": args.max_iterations,
              "potential_kj_mol_before": before, "potential_kj_mol_after": after,
              "scope": "exploratory minimized solvated complex; potential energy is not binding energy",
              "limitations": ["Product-templated pose unvalidated", "Water hydrogen orientations initially arbitrary" if crystal_water_count else "Crystal waters not retained", "Only a single protonation hypothesis", "No equilibration or binding-free-energy calculation"]}
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
