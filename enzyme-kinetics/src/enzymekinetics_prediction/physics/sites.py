"""Evidence-labelled active-site residues, conformations, contacts and assembly handling."""

from pathlib import Path

from Bio.PDB import PDBIO, PDBParser, Structure, Model, Chain
from Bio.PDB.vectors import calc_dihedral
import numpy as np

from .artifacts import artifact, digest_file

KSI_ROLES = {
    14: ("TYR", "oxyanion_hole_donor"),
    38: ("ASP", "general_base_acid"),
    99: ("ASP", "oxyanion_hole_donor"),
    54: ("PHE", "base_positioning_and_hydrophobic_binding"),
    116: ("PHE", "base_positioning_and_hydrophobic_binding"),
    39: ("PRO", "backbone_positioning_of_Asp38"),
}
KSI_BINDING_CANDIDATES = {18, 55, 58, 61, 63, 82, 84, 86, 95, 97, 112, 114, 116, 121}
CHI1 = {
    "TYR": ["N", "CA", "CB", "CG"],
    "ASP": ["N", "CA", "CB", "CG"],
    "ASN": ["N", "CA", "CB", "CG"],
    "PHE": ["N", "CA", "CB", "CG"],
    "PRO": ["N", "CA", "CB", "CG"],
}
CHI2 = {
    "TYR": ["CA", "CB", "CG", "CD1"],
    "PHE": ["CA", "CB", "CG", "CD1"],
    "ASP": ["CA", "CB", "CG", "OD1"],
    "ASN": ["CA", "CB", "CG", "OD1"],
}
CHI1.update(
    {
        name: ["N", "CA", "CB", last]
        for name, last in {
            "ARG": "CG",
            "CYS": "SG",
            "GLU": "CG",
            "GLN": "CG",
            "HIS": "CG",
            "ILE": "CG1",
            "LEU": "CG",
            "LYS": "CG",
            "MET": "CG",
            "SER": "OG",
            "THR": "OG1",
            "TRP": "CG",
            "VAL": "CG1",
        }.items()
    }
)
CHI2.update(
    {
        name: ["CA", "CB", "CG", last]
        for name, last in {
            "ARG": "CD",
            "GLU": "CD",
            "GLN": "CD",
            "HIS": "ND1",
            "LEU": "CD1",
            "LYS": "CD",
            "MET": "SD",
            "TRP": "CD1",
        }.items()
    }
)
CHI2["ILE"] = ["CA", "CB", "CG1", "CD1"]


def heavy_atoms(residue):
    return [a for a in residue.get_atoms() if a.element not in {"H", "D"}]


def dihedral(residue, names):
    if not all(n in residue for n in names):
        return None
    return float(np.degrees(calc_dihedral(*[residue[n].get_vector() for n in names])))


def flatten_assembly(source, output):
    """RCSB PDB biological assemblies can place symmetry mates in separate MODELs."""
    structure = PDBParser(QUIET=True).get_structure("assembly", source)
    combined = Structure.Structure("normalized")
    model = Model.Model(0)
    combined.add(model)
    mapping = []
    for old_model in structure:
        for old_chain in old_model:
            residues = [r for r in old_chain if r.id[0] == " "]
            if not residues:
                continue
            if len(mapping) >= 26:
                raise ValueError("Too many chains for normalized PDB")
            chain = Chain.Chain(chr(65 + len(mapping)))
            model.add(chain)
            for residue in residues:
                cleaned = residue.copy()
                for atom in list(cleaned):
                    if atom.element in {"H", "D"}:
                        cleaned.detach_child(atom.id)
                chain.add(cleaned)
            mapping.append(
                {
                    "model": old_model.id,
                    "chain": old_chain.id,
                    "new_chain": chain.id,
                    "residues": len(residues),
                }
            )
    if not mapping:
        raise ValueError("No protein chains in assembly")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(combined)
    io.save(str(output))
    return {
        "source_sha256": digest_file(source),
        "normalized_sha256": digest_file(output),
        "chain_mapping": mapping,
        "selection": "protein heavy atoms only; hydrogens reassigned for stated pH; crystal waters/additives excluded explicitly",
    }


def site_report(
    pdb, chain_id="A", ligand_name=None, contact_cutoff_A=4.5, roles=KSI_ROLES
):
    if not 0 < contact_cutoff_A <= 10:
        raise ValueError("Invalid contact cutoff")
    model = PDBParser(QUIET=True).get_structure("site", pdb)[0]
    chain = model[chain_id]
    ligand = (
        [r for r in chain if r.resname == ligand_name and r.id[0] != " "]
        if ligand_name
        else []
    )
    if ligand_name and not ligand:
        raise ValueError("Requested ligand is absent from selected chain")
    ligand_xyz = np.array([a.coord for r in ligand for a in heavy_atoms(r)])
    protein = [r for r in chain if r.id[0] == " "]
    neighbors = {
        r.id: (
            protein[i - 1] if i else None,
            protein[i + 1] if i + 1 < len(protein) else None,
        )
        for i, r in enumerate(protein)
    }

    def backbone_angle(specification):
        if any(r is None or name not in r for r, name in specification):
            return None
        return float(
            np.degrees(
                calc_dihedral(*[r[name].get_vector() for r, name in specification])
            )
        )

    residues = []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        xyz = np.array([a.coord for a in heavy_atoms(residue)])
        distance = (
            float(np.linalg.norm(xyz[:, None] - ligand_xyz[None, :], axis=-1).min())
            if len(ligand_xyz)
            else None
        )
        role = roles.get(residue.id[1])
        if role is None and (distance is None or distance > contact_cutoff_A):
            continue
        record = {
            "chain": chain.id,
            "author_residue_number": residue.id[1],
            "insertion_code": residue.id[2].strip(),
            "resname": residue.resname,
            "literature_wt_resname": role[0] if role else None,
            "literature_role": role[1] if role else "geometric_ligand_contact",
            "matches_expected_wt": residue.resname == role[0] if role else None,
            "min_ligand_distance_A": distance,
            "chi1_deg": dihedral(residue, CHI1[residue.resname])
            if residue.resname in CHI1
            else None,
            "chi2_deg": dihedral(residue, CHI2[residue.resname])
            if residue.resname in CHI2
            else None,
            "heavy_atom_coordinates_A": {
                a.name: a.coord.astype(float).tolist() for a in heavy_atoms(residue)
            },
        }
        previous, following = neighbors[residue.id]
        record["phi_deg"] = backbone_angle(
            [(previous, "C"), (residue, "N"), (residue, "CA"), (residue, "C")]
        )
        record["psi_deg"] = backbone_angle(
            [(residue, "N"), (residue, "CA"), (residue, "C"), (following, "N")]
        )
        record["peptide_omega_deg"] = backbone_angle(
            [(previous, "CA"), (previous, "C"), (residue, "N"), (residue, "CA")]
        )
        residues.append(record)
    pairs = []
    selected = [r for r in chain if r.id[0] == " " and r.id[1] in roles]
    for i, a in enumerate(selected):
        for b in selected[i + 1 :]:
            xa = np.array([x.coord for x in heavy_atoms(a)])
            xb = np.array([x.coord for x in heavy_atoms(b)])
            pairs.append(
                {
                    "residue_a": a.id[1],
                    "residue_b": b.id[1],
                    "min_heavy_distance_A": float(
                        np.linalg.norm(xa[:, None] - xb[None, :], axis=-1).min()
                    ),
                }
            )
    return artifact(
        "structure-active-site",
        {"pdb_sha256": digest_file(pdb)},
        {
            "chain": chain_id,
            "ligand": ligand_name,
            "contact_cutoff_A": contact_cutoff_A,
        },
        {"residues": residues, "pair_contacts": pairs},
        {
            "residues": "mixed explicit Angstrom/degree records",
            "pair_contacts": "Angstrom",
        },
        limitations=[
            "Residue roles are literature annotations, not inferred from distance alone",
            "Contacts with EQU/5SD are analog binding; D38N/Y14F structures are not WT ES",
            "Crystallographic B factors and AlphaFold confidence must not be interpreted as equivalent flexibility",
        ],
    )
