"""Static molecular coordinates with an orbiting camera; NOT an MD animation."""
import argparse
import hashlib
import json
from io import StringIO
from pathlib import Path

import imageio.v2 as imageio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBParser


def collect(structure):
    model = structure[0]
    backbone = {}
    for chain in ("A", "B"):
        backbone[chain] = np.array([r["CA"].coord for r in model[chain] if r.id[0] == " " and "CA" in r])
    lig = [a for a in model.get_atoms() if a.get_parent().resname == "UNK" and a.element in ("C", "O")]
    if any(len(x) != 125 for x in backbone.values()) or len(lig) != 21:
        raise ValueError("Unexpected WT dimer or ligand heavy-atom count")
    return backbone, np.array([a.coord for a in lig]), [a.element for a in lig]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pdb", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--frames", type=int, default=24)
    args = p.parse_args()
    if args.frames < 2 or (args.output.exists() and any(args.output.iterdir())):
        raise ValueError("Use >=2 frames and an empty output directory")
    lines = args.pdb.read_text().splitlines()
    selected = [line for line in lines if line.startswith("ATOM  ") or (line.startswith("HETATM") and line[17:20] == "UNK")]
    backbone, ligand, elements = collect(PDBParser(QUIET=True).get_structure("p", StringIO("\n".join(selected) + "\nEND\n")))
    args.output.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(6, 6), dpi=100)
    ax = fig.add_subplot(projection="3d")
    for chain, xyz in backbone.items():
        ax.plot(*xyz.T, lw=1.5, label=f"chain {chain}")
    ax.scatter(*ligand.T, c=["red" if x == "O" else "orange" for x in elements], s=25, depthshade=False, label="5-AND")
    allxyz = np.vstack(list(backbone.values()) + [ligand])
    center = (allxyz.min(0) + allxyz.max(0)) / 2
    radius = float((allxyz.max(0) - allxyz.min(0)).max() / 2)
    ax.set_xlim(center[0]-radius, center[0]+radius)
    ax.set_ylim(center[1]-radius, center[1]+radius)
    ax.set_zlim(center[2]-radius, center[2]+radius)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()
    ax.set_title("WT–5-AND model: static coordinates, rotating camera")
    ax.legend(loc="upper left")
    frames = []
    for i in range(args.frames):
        ax.view_init(elev=20, azim=360 * i / args.frames)
        fig.canvas.draw()
        frames.append(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy())
    fig.savefig(args.output / "preview.png")
    plt.close(fig)
    imageio.mimwrite(args.output / "camera-orbit.mp4", frames, fps=8, codec="libx264", macro_block_size=16)
    report = {"pdb": str(args.pdb), "pdb_sha256": hashlib.sha256(args.pdb.read_bytes()).hexdigest(),
              "frames": args.frames, "fps": 8, "chain_CA_counts": {k: len(v) for k, v in backbone.items()},
              "ligand_heavy_atoms": len(ligand), "scope": "static model with rotating camera, not reaction dynamics or MD trajectory"}
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
