"""Synthetic DCD integration and analysis output isolation, without simulation."""

import json

import numpy as np
import pandas as pd
import pytest

from enzymekinetics_prediction.physics.geometry import STRUCTURAL_METHOD
from enzymekinetics_prediction.physics.trajectory import protein_structure, load_dcd
from enzymekinetics_prediction.research.ensemble import analyze_replica

md = pytest.importorskip("mdtraj")


def trajectory():
    topology = md.Topology()
    for _ in range(2):
        chain = topology.add_chain()
        previous = None
        for resid in (14, 38, 39, 2):
            residue = topology.add_residue("ALA", chain, resSeq=resid)
            atom = topology.add_atom("CA", md.element.carbon, residue)
            bridge = topology.add_atom("C", md.element.carbon, residue)
            topology.add_bond(atom, bridge)
            if previous is not None:
                topology.add_bond(previous, atom)
            previous = bridge
    tetra = np.array([[0, 0, 0], [0.25, 0, 0], [0, 0.5, 0], [0, 0, 0.75]])
    ca_positions = np.vstack((tetra, tetra + [1.5, 0.25, 0.125])) + 2
    positions = np.array(
        [p + offset for p in ca_positions for offset in ([0, 0, 0], [0.03125, 0, 0])]
    )
    xyz = np.repeat(positions[None], 6, axis=0)
    xyz[:, 6, 2] += np.arange(6) / 64
    return md.Trajectory(
        xyz,
        topology,
        unitcell_lengths=np.full((6, 3), 10.0),
        unitcell_angles=np.full((6, 3), 90.0),
    )


def test_full_bonded_reconstruction_precedes_ca_alignment():
    clean = trajectory()
    moved = clean[:]
    images = np.random.default_rng(202).integers(-2, 3, size=moved.xyz.shape)
    moved.xyz += 10 * images
    _, expected = protein_structure(clean)
    _, actual = protein_structure(moved)
    for key in ("rmsd_nm", "rmsf_nm", "dccm"):
        np.testing.assert_allclose(actual[key], expected[key], atol=1e-12)


def test_missing_connectivity_and_triclinic_metadata_are_not_silently_ignored():
    traj = trajectory()
    traj.unitcell_angles = np.full((6, 3), [90, 90, 60])
    with pytest.raises(ValueError, match="orthorhombic"):
        protein_structure(traj)
    topo = md.Topology()
    chain = topo.add_chain()
    for i in range(2):
        topo.add_atom("CA", md.element.carbon, topo.add_residue("ALA", chain))
    incomplete = md.Trajectory(
        np.zeros((2, 2, 3)),
        topo,
        unitcell_lengths=np.full((2, 3), 10.0),
        unitcell_angles=np.full((2, 3), 90.0),
    )
    with pytest.raises(ValueError, match="disconnected"):
        protein_structure(incomplete)
    incomplete.unitcell_vectors = None
    with pytest.raises(ValueError, match="box vectors"):
        protein_structure(incomplete, periodic=True)


def test_replica_outputs_are_new_and_share_canonical_metrics(tmp_path):
    source = tmp_path / "replica-101"
    source.mkdir()
    traj = trajectory()
    traj.save_dcd(str(source / "trajectory.dcd"))
    traj[0].save_pdb(str(source / "topology.pdb"))
    # PDB cannot encode the arbitrary coarse-grained bonds: record CONECT.
    serials = [
        int(line[6:11])
        for line in (source / "topology.pdb").read_text().splitlines()
        if line.startswith("ATOM")
    ]
    with (source / "topology.pdb").open("a") as handle:
        for a, b in traj.topology.bonds:
            handle.write(f"CONECT{serials[a.index]:5d}{serials[b.index]:5d}\n")
    pd.DataFrame({"Temperature (K)": [298.15] * 6}).to_csv(
        source / "state.csv", index=False
    )
    historical = source / "analysis.json"
    historical.write_text('{"historical": true}\n')
    result = analyze_replica(source)
    assert historical.read_text() == '{"historical": true}\n'
    output = source / "analysis-c1"
    assert result["structural_method"] == STRUCTURAL_METHOD
    assert json.loads((output / "analysis.json").read_text())["frames"] == 6
    _, expected = protein_structure(
        load_dcd(source / "trajectory.dcd", source / "topology.pdb")
    )
    with np.load(output / "structure.npz") as actual:
        for key in ("rmsd_nm", "rmsf_nm", "dccm"):
            np.testing.assert_allclose(actual[key], expected[key], atol=1e-12)
    with pytest.raises(FileExistsError):
        analyze_replica(source)


def trace_system(path, periodic, atoms):
    mm = pytest.importorskip("openmm")
    system = mm.System()
    force = mm.NonbondedForce()
    if periodic:
        force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
    for i in range(atoms):
        system.addParticle(12)
        force.addParticle(1.0 if i == atoms - 1 else 0.0, 0.1, 0)
    system.addForce(force)
    path.write_text(mm.XmlSerializer.serialize(system))


def test_periodic_trace_requires_explicit_nonperiodic_preparation(
    monkeypatch, tmp_path
):
    pytest.importorskip("openmm")
    from enzymekinetics_prediction.physics import trajectory as module

    traj = trajectory()
    system = tmp_path / "system.xml"
    trace_system(system, True, traj.n_atoms)
    monkeypatch.setattr(module, "load_dcd", lambda *args: traj)
    with pytest.raises(ValueError, match="prepared nonperiodic"):
        module.analyze_trajectory(
            "unused.dcd",
            "unused.pdb",
            system,
            [0, 1],
            "protein",
            "test",
            tmp_path / "trace",
        )


def test_nonperiodic_trace_accepts_default_dcd_box_without_new_flag(tmp_path):
    pytest.importorskip("openmm")
    from enzymekinetics_prediction.physics.trajectory import analyze_trajectory

    traj = trajectory()
    dcd, pdb, system = (
        tmp_path / name for name in ("trajectory.dcd", "topology.pdb", "system.xml")
    )
    traj.save_dcd(str(dcd))
    traj[0].save_pdb(str(pdb))
    trace_system(system, False, traj.n_atoms)
    result = analyze_trajectory(
        dcd,
        pdb,
        system,
        [0, 1],
        "protein",
        "test",
        tmp_path / "trace",
        ionic_strength_M=0,
        dielectric=1,
    )
    assert result["values"]["frames"] == 6
    assert (tmp_path / "trace/trace.csv").is_file()
