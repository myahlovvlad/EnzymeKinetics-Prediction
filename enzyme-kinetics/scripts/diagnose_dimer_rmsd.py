"""Compare raw dimer RMSD against a minimum-image chain reconstruction.

Diagnostic assumes an orthorhombic box and two already-whole protein chains;
it does not prove that the resulting dimer is physically stable.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import mdtraj as md
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from enzymekinetics_prediction.physics.trajectory import load_dcd, protein_structure  # noqa: E402
from enzymekinetics_prediction.physics.geometry import (  # noqa: E402
    reconstruct_assembly,
    structural_metrics,
    STRUCTURAL_METHOD,
)


def reconstruct_pair(xyz, box, chain_a, chain_b):
    """Compatibility wrapper: input chains must already be whole."""
    fixed, applied = reconstruct_assembly(xyz, box, [chain_a, chain_b])
    return fixed, -applied[:, 1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    traj = load_dcd(args.trajectory, args.topology)
    ca = traj.atom_slice(traj.topology.select("protein and name CA"))
    chains = sorted(set(a.residue.chain.index for a in ca.topology.atoms))
    if len(chains) != 2 or ca.unitcell_lengths is None:
        raise ValueError("Expected exactly two protein chains and periodic box")
    if not np.allclose(ca.unitcell_angles, 90, atol=1e-3):
        raise ValueError("Non-orthorhombic box unsupported")
    a = np.array(
        [x.index for x in ca.topology.atoms if x.residue.chain.index == chains[0]]
    )
    b = np.array(
        [x.index for x in ca.topology.atoms if x.residue.chain.index == chains[1]]
    )
    raw = ca.xyz.copy()
    protein, canonical = protein_structure(traj)
    fixed = protein.xyz[:, protein.topology.select("name CA")]
    # Retain legacy diagnostic shift convention: subtracted from chain B.
    _, shifts = reconstruct_pair(raw, ca.unitcell_lengths, a, b)

    def rmsd(xyz, select=None):
        selected = xyz if select is None else xyz[:, select]
        return structural_metrics(selected)["rmsd_nm"]

    original, corrected = rmsd(raw), canonical["rmsd_nm"]
    separate_a, separate_b = rmsd(raw, a), rmsd(raw, b)
    mid = len(original) // 2
    result = {
        "inputs": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (args.trajectory, args.topology)
        },
        "mdtraj_version": md.__version__,
        "structural_method": STRUCTURAL_METHOD,
        "frames": len(raw),
        "ca_per_chain": [len(a), len(b)],
        "box_lengths_first_nm": ca.unitcell_lengths[0].tolist(),
        "shifted_frames": int(np.count_nonzero(np.linalg.norm(shifts, axis=1) > 0)),
        "raw_rmsd_half_means_nm": [
            float(original[:mid].mean()),
            float(original[mid:].mean()),
        ],
        "reconstructed_rmsd_half_means_nm": [
            float(corrected[:mid].mean()),
            float(corrected[mid:].mean()),
        ],
        "chain_a_rmsd_half_means_nm": [
            float(separate_a[:mid].mean()),
            float(separate_a[mid:].mean()),
        ],
        "chain_b_rmsd_half_means_nm": [
            float(separate_b[:mid].mean()),
            float(separate_b[mid:].mean()),
        ],
        "raw_chain_center_distances_nm": np.linalg.norm(
            raw[:, b].mean(1) - raw[:, a].mean(1), axis=1
        ).tolist(),
        "reconstructed_chain_center_distances_nm": np.linalg.norm(
            fixed[:, b].mean(1) - fixed[:, a].mean(1), axis=1
        ).tolist(),
        "shifts_nm": shifts.tolist(),
        "limitations": [
            "Bonded chains are reconstructed before compact minimum-image assembly; no stability claim.",
            "Only CA-based alignment and displacement examined.",
        ],
    }
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
