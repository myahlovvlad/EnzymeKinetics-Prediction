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
from enzymekinetics_prediction.physics.trajectory import load_dcd  # noqa: E402


def reconstruct_pair(xyz, box, chain_a, chain_b):
    out = np.array(xyz, dtype=float, copy=True)
    box = np.asarray(box, dtype=float)
    if box.shape != (len(out), 3) or np.any(box <= 0):
        raise ValueError("Expected positive orthorhombic box lengths for every frame")
    a = out[:, chain_a, :].mean(axis=1)
    b = out[:, chain_b, :].mean(axis=1)
    shifts = np.round((b - a) / box) * box
    out[:, chain_b, :] -= shifts[:, None, :]
    return out, shifts


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
    fixed, shifts = reconstruct_pair(raw, ca.unitcell_lengths, a, b)

    def rmsd(xyz, select=None):
        tmp = md.Trajectory(xyz.astype(np.float32), ca.topology)
        return md.rmsd(tmp, tmp, 0, atom_indices=select, parallel=False).astype(float)

    original, corrected = rmsd(raw), rmsd(fixed)
    separate_a, separate_b = rmsd(raw, a), rmsd(raw, b)
    mid = len(original) // 2
    result = {
        "inputs": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (args.trajectory, args.topology)
        },
        "mdtraj_version": md.__version__,
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
            "Assumes each chain whole and minimum-image dimer; no physical stability claim.",
            "Only CA-based alignment and displacement examined.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
