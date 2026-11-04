import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

openmm = pytest.importorskip("openmm")
from openmm import unit

scripts = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(scripts))
try:
    spec = importlib.util.spec_from_file_location("lowT_quench", scripts / "lowT_quench.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
finally:
    sys.path.remove(str(scripts))


def test_checks_reject_bad_frame():
    x = np.array([[0., 0., 0.], [.1, 0., 0.]])
    f = np.zeros_like(x)
    box = np.eye(3)
    m.assert_valid_frame(x, f, 1., 2., 0., 1e-5, 0., 0., box, box)
    with pytest.raises(ValueError, match="constraint"):
        m.assert_valid_frame(x, f, 1., 2., 1e-3, 1e-5, 0., 0., box, box)
    with pytest.raises(ValueError, match="RMSD"):
        m.assert_valid_frame(x, f, 1., 2., 0., 1e-5, .5, 0., box, box)
    with pytest.raises(ValueError, match="nonfinite"):
        m.assert_valid_frame(x, f, np.nan, 2., 0., 1e-5, 0., 0., box, box)


def test_end_checks_reject_drift():
    rows = [{"potential_energy_kJ_mol": 1000-i*200, "temperature_K": 298.15} for i in range(5)]
    assert not m.end_checks(rows, 298.15)["energy_stable"]
    rows = [{"potential_energy_kJ_mol": 1000+i, "temperature_K": 5000} for i in range(5)]
    assert not m.end_checks(rows, 298.15)["temperature_stable"]


def test_periodic_protein_ligand_hbond_candidate():
    x = np.array([[1.90, 0, 0], [1.99, 0, 0], [.15, 0, 0]])
    assert m.protein_ligand_hbonds(x, [2., 2., 2.], [(0, 1)], [2], ["protein", "protein", "ligand"]) == 1
    assert m.protein_ligand_hbonds(x, [2., 2., 2.], [(0, 1)], [2], ["protein", "protein", "protein"]) == 0


def test_restraint_energy_force_and_periodic_image():
    system = openmm.System()
    system.addParticle(12.)
    system.setDefaultPeriodicBoxVectors(openmm.Vec3(2,0,0), openmm.Vec3(0,2,0), openmm.Vec3(0,0,2))
    _, info = m.add_backbone_restraint(system, [0], [[0., 0., 0.]], 5.)
    assert info["atom_indices"] == [0]
    integ = openmm.VerletIntegrator(.001)
    context = openmm.Context(system, integ, openmm.Platform.getPlatformByName("CPU"), {"Threads": "4"})
    context.setPositions([[.1, 0, 0]] * unit.nanometer)
    state = context.getState(getEnergy=True, getForces=True)
    assert state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) == pytest.approx(20.92)
    assert state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)[0, 0] == pytest.approx(-418.4)
    context.setPositions([[2.1, 0, 0]] * unit.nanometer)
    assert context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) == pytest.approx(20.92)
    with pytest.raises(ValueError):
        m.add_backbone_restraint(system, [0, 0], [[0., 0., 0.]], 5.)


def test_first_failing_rolling_energy_window():
    energies = [-555212.26, -551989.37, -552291.39, -552453.93, -552763.59, -552961.10]
    rows = [{"potential_energy_kJ_mol": e, "temperature_K": 298.15} for e in energies]
    assert m.end_checks(rows[:5], 298.15)["energy_stable"]
    assert not m.end_checks(rows[:6], 298.15)["energy_stable"]


def test_shared_md_loop_stops_on_first_failing_window(tmp_path):
    # Minimal peptide-like two-particle fragment and a rigid three-atom water.
    # The time-dependent offset is deliberate: it creates a deterministic
    # monotone potential without changing forces or destabilizing the system.
    system = openmm.System()
    for mass in (12., 14., 16., 1., 1.):
        system.addParticle(mass)
    bonds = openmm.HarmonicBondForce()
    bonds.addBond(0, 1, .14, 1000.)
    system.addForce(bonds)
    for i, j, distance in ((2, 3, .09572), (2, 4, .09572), (3, 4, .15139)):
        system.addConstraint(i, j, distance)
    bias = openmm.CustomExternalForce("offset")
    bias.addGlobalParameter("offset", 0.)
    bias.addParticle(0, [])
    system.addForce(bias)
    integ = openmm.VerletIntegrator(.001 * unit.picoseconds)
    context = openmm.Context(system, integ, openmm.Platform.getPlatformByName("CPU"), {"Threads": "4"})
    context.setPositions([[0, 0, 0], [.14, 0, 0], [1, 1, 1], [1.09572, 1, 1], [.976, 1.09266, 1]] * unit.nanometer)
    context.applyConstraints(1e-8)
    context.setVelocitiesToTemperature(10 * unit.kelvin, 42)
    metrics = []
    saved = []
    with pytest.raises(ValueError, match="rolling energy"):
        for block in m.checked_blocks(20, metrics, 298.15):
            if block:
                integ.step(1)
                context.setParameter("offset", -200. * block)
            state = context.getState(getPositions=True, getEnergy=True)
            (tmp_path / f"state-{block}.xml").write_text(openmm.XmlSerializer.serialize(state))
            saved.append(block)
            metrics.append({"potential_energy_kJ_mol": state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole), "temperature_K": 298.15})
    assert saved == list(range(5))
    assert not (tmp_path / "state-5.xml").exists()
    assert len(metrics) == 5


def test_state_roundtrip(tmp_path):
    system = openmm.System()
    system.addParticle(12.)
    integ = openmm.VerletIntegrator(.001)
    context = openmm.Context(system, integ, openmm.Platform.getPlatformByName("CPU"), {"Threads": "4"})
    context.setPositions([[.1, .2, .3]] * unit.nanometer)
    state = context.getState(getPositions=True)
    (tmp_path / "state.xml").write_text(openmm.XmlSerializer.serialize(state))
    restored = openmm.XmlSerializer.deserialize((tmp_path / "state.xml").read_text())
    coords = np.asarray(restored.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    np.save(tmp_path / "positions.npy", coords)
    np.save(tmp_path / "box.npy", np.asarray(restored.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)))
    np.testing.assert_allclose(np.load(tmp_path / "positions.npy"), [[.1, .2, .3]])
    assert np.load(tmp_path / "box.npy").shape == (3, 3)
    assert "temperature_K" in m.METRICS and "projected_force_max_kJ_mol_nm" in m.METRICS
