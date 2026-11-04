"""Analytic reference cases and native-backend contracts for the physical workflow."""

import json

import numpy as np
import pandas as pd
import pytest
from scipy.constants import Boltzmann, elementary_charge

from enzymekinetics_prediction.physics.artifacts import (
    artifact,
    content_key,
    write_json,
)
from enzymekinetics_prediction.physics.dynamics import MDSettings, ensemble_observables
from enzymekinetics_prediction.physics.electrostatics import screened_coulomb, sample_dx
from enzymekinetics_prediction.physics.kinetics import CycleRates, eyring, tunneling
from enzymekinetics_prediction.physics.pipeline import load_config, run_pipeline
from enzymekinetics_prediction.physics.tables import key_for, load_tables, write_table


def test_cycle_flux_matches_michaelis_menten_and_release_limit():
    rates = CycleRates(1e7, 200, 100, 30, 10)
    parameters = rates.parameters()
    assert parameters["kcat_s_inv"] < 10
    for substrate in [0, 1e-9, 1e-6, 1, 100]:
        stationary = rates.stationary(substrate, 2e-6)
        fractions = stationary["E_ES_EP_fractions"]
        assert sum(fractions) == pytest.approx(1)
        expected = (
            2e-6
            * parameters["kcat_s_inv"]
            * substrate
            / (parameters["Km_M"] + substrate)
        )
        assert stationary["velocity_M_per_s"] == pytest.approx(expected, abs=1e-15)
    fast_release = CycleRates(1e7, 200, 100, 0, 1e12).parameters()
    assert fast_release["kcat_s_inv"] == pytest.approx(100)
    assert fast_release["Km_M"] == pytest.approx(300 / 1e7)


def test_tst_dimensions_and_tunneling_domain():
    assert eyring(60, 298.15) == pytest.approx(191.3008269, rel=1e-8)
    assert eyring(65, 298.15) < eyring(60, 298.15)
    assert tunneling(10, 298.15, "bell-parabolic") == pytest.approx(
        tunneling(10, 298.15), rel=1e-7
    )
    with pytest.raises(ValueError):
        tunneling(3000, 100, "bell-parabolic")
    with pytest.raises(ValueError):
        CycleRates(0, 1, 1, 0, 1)


def test_point_charge_units_symmetry_screening_and_self_exclusion():
    phi, field = screened_coulomb(
        [[0, 0, 0]], [1], [[1, 0, 0], [-1, 0, 0]], dielectric=1, ionic_strength_M=0
    )
    np.testing.assert_allclose(phi, 14.39964548, rtol=1e-7)
    np.testing.assert_allclose(field[:, 0], [14.39964548, -14.39964548], rtol=1e-7)
    np.testing.assert_allclose(field[:, 1:], 0)
    screened, _ = screened_coulomb(
        [[0, 0, 0]], [1], [[1, 0, 0]], dielectric=1, ionic_strength_M=0.15
    )
    assert 0 < screened[0] < phi[0]
    with pytest.raises(ValueError):
        screened_coulomb([[0, 0, 0]], [1], [[0, 0, 0]])


def test_dx_linear_potential_axis_order_and_thermal_units(tmp_path):
    # z is fastest; potential in kBT/e is x + 2y + 3z.
    values = [x + 2 * y + 3 * z for x in range(2) for y in range(2) for z in range(2)]
    dx = tmp_path / "linear.dx"
    dx.write_text(
        "object 1 class gridpositions counts 2 2 2\norigin 0 0 0\ndelta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\nobject 3 class array type double rank 0 items 8 data follows\n"
        + " ".join(map(str, values))
        + "\n"
    )
    phi, field = sample_dx(dx, [[0.5, 0.5, 0.5]], 300)
    thermal = Boltzmann * 300 / elementary_charge
    np.testing.assert_allclose(phi, [3 * thermal])
    np.testing.assert_allclose(field, [[-thermal, -2 * thermal, -3 * thermal]])
    with pytest.raises(ValueError):
        sample_dx(dx, [[2, 0, 0]], 300)


def test_ensemble_rigid_motion_and_distance_population():
    xyz = np.random.default_rng(42).normal(size=(5, 8, 3))
    values, _ = ensemble_observables(xyz, [(0, 1)], [2])
    rotation = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    other, _ = ensemble_observables(xyz @ rotation + 10, [(0, 1)], [2])
    assert other == pytest.approx(values, abs=1e-10)
    expected = (np.linalg.norm(xyz[:, 0] - xyz[:, 1], axis=1) <= 2).mean()
    assert values["productive_fraction"] == expected
    with pytest.raises(ValueError):
        MDSettings(production_steps=1)


def test_table_provenance_integrity_and_condition_join(tmp_path):
    table = tmp_path / "electro.csv"
    frame = pd.DataFrame(
        {
            "uniprot": ["P1"],
            "ph": [7.0],
            "ionic_strength_M": [0.15],
            "electro_potential": [0.1],
        }
    )
    write_table(
        table,
        frame,
        "electro",
        ["uniprot", "ph", "ionic_strength_M"],
        {"electro_potential": "Volt"},
        {"test": True},
    )
    loaded = load_tables([table])[0]
    assert loaded["lookup"][key_for(frame.iloc[0], loaded["metadata"]["join_on"])][
        0
    ] == pytest.approx(0.1)
    missing_ionic = key_for({"uniprot": "P1", "ph": 7}, loaded["metadata"]["join_on"])
    assert missing_ionic not in loaded["lookup"]
    table.write_text(table.read_text().replace("0.1", "0.2"))
    with pytest.raises(ValueError, match="checksum"):
        load_tables([table])


def test_artifacts_units_atomic_json_and_cache_identity(tmp_path):
    assert content_key({"a": 1, "b": 2}) == content_key({"b": 2, "a": 1})
    assert content_key({"T": 300}) != content_key({"T": 301})
    with pytest.raises(ValueError):
        artifact("test", {}, {}, {"energy": 1}, {})
    target = tmp_path / "a.json"
    write_json(target, {"valid": True})
    with pytest.raises(ValueError):
        write_json(target, {"invalid": float("nan")})
    assert json.loads(target.read_text()) == {"valid": True}


def test_pipeline_rejects_unknown_fields_and_records_stage_failure(tmp_path):
    pdb = tmp_path / "fake.pdb"
    pdb.write_text("END\n")
    config = {
        "schema_version": 1,
        "protein_id": "P0",
        "protein_pdb": str(pdb),
        "substrate_smiles": "CCO",
        "temperature_K": 298.15,
        "ph": 7,
        "quantum": {"enabled": False},
        "pocket": {"enabled": False},
        "electrostatics": {"enabled": True},
        "dynamics": {"enabled": False},
        "kinetics": {"enabled": False},
    }
    source = tmp_path / "config.json"
    source.write_text(json.dumps(config | {"invented": True}))
    with pytest.raises(ValueError, match="Unknown"):
        load_config(source)
    source.write_text(json.dumps(config))
    output = tmp_path / "run"
    run_pipeline(load_config(source), output)
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["stages"]["electrostatics"]["status"] == "failed"
    with pytest.raises(ValueError, match="immutable"):
        run_pipeline(load_config(source), output)


def test_native_tblite_charge_fukui_and_cache(tmp_path):
    pytest.importorskip("tblite")
    from enzymekinetics_prediction.physics.quantum import (
        QuantumSettings,
        molecular_descriptors,
    )

    settings = QuantumSettings(fukui=True, threads=1)
    first = molecular_descriptors("O", settings, tmp_path)
    second = molecular_descriptors("O", settings, tmp_path)
    assert first == second
    assert first["backend"] == "tblite"
    assert sum(first["values"]["atomic_charges_e"]) == pytest.approx(0, abs=1e-5)
    assert sum(first["values"]["fukui_plus"]) == pytest.approx(1, abs=1e-5)
    assert sum(first["values"]["fukui_minus"]) == pytest.approx(1, abs=1e-5)


def test_reversible_reaction_conservation_equilibrium_and_haldane():
    from enzymekinetics_prediction.physics.reaction import ReversibleCycle

    rates = ReversibleCycle(1e7, 200, 100, 30, 1e6, 10)
    e, s = 1e-6, 1e-5
    es = rates.kon_substrate_M_inv_s * e * s / rates.koff_substrate_s
    ep = rates.chemistry_forward_s * es / rates.chemistry_reverse_s
    p = s * rates.equilibrium_constant()
    np.testing.assert_allclose(rates.fluxes([e, es, ep, s, p]), 0, atol=1e-15)
    frame, evidence = rates.simulate([1e-6, 0, 0, 1e-4, 0], np.linspace(0, 30, 100))
    np.testing.assert_allclose(
        frame[["E_M", "ES_M", "EP_M"]].sum(axis=1), 1e-6, atol=1e-14
    )
    np.testing.assert_allclose(
        frame[["S_M", "ES_M", "EP_M", "P_M"]].sum(axis=1), 1e-4, atol=1e-14
    )
    assert evidence["values"]["substrate_conservation_max_error_M"] < 1e-14
    with pytest.raises(ValueError, match="Haldane"):
        rates.check_thermodynamics(100, 300)


def test_projected_trajectory_field_orientation_and_distances():
    from enzymekinetics_prediction.physics.trajectory import projected_fields

    # Probe midpoint x=1 A; positive point charge at origin, vacuum.
    xyz = np.array([[[0.05, 0, 0], [0.15, 0, 0], [0, 0, 0]]] * 2)
    forward = projected_fields(
        xyz, [0, 0, 1], [0, 1], [2], dielectric=1, ionic_strength_M=0
    )
    backward = projected_fields(
        xyz, [0, 0, 1], [1, 0], [2], dielectric=1, ionic_strength_M=0
    )
    np.testing.assert_allclose(forward.distance_nm, 0.1)
    np.testing.assert_allclose(forward.field_projection_V_per_A, 14.39964548, rtol=1e-7)
    np.testing.assert_allclose(
        forward.field_projection_V_per_A, -backward.field_projection_V_per_A
    )
    with pytest.raises(ValueError, match="exclude"):
        projected_fields(xyz, [0, 0, 1], [0, 1], [0, 2])


def test_native_openmm_prepared_system_in_unicode_workspace(tmp_path):
    mm = pytest.importorskip("openmm")
    pytest.importorskip("mdtraj")
    from enzymekinetics_prediction.physics.dynamics import run_md

    folder = tmp_path / "траектории"
    folder.mkdir()
    pdb = folder / "protein.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00  0.00           C\nATOM      3  CA  ALA A   3       7.600   0.000   0.000  1.00  0.00           C\nEND\n"
    )
    system = mm.System()
    for _ in range(3):
        system.addParticle(12)
    bond = mm.HarmonicBondForce()
    bond.addBond(0, 1, 0.38, 1000)
    bond.addBond(1, 2, 0.38, 1000)
    system.addForce(bond)
    xml = folder / "system.xml"
    xml.write_text(mm.XmlSerializer.serialize(system))
    settings = MDSettings(
        equilibration_steps=0, production_steps=6, report_interval=2, threads=1
    )
    result = run_md(pdb, folder / "md", settings, xml)
    assert result["values"]["frames"] == 3
    assert result["values"]["production_time_ps"] == pytest.approx(0.012)
    assert (folder / "md/checkpoint.chk").is_file()


def test_pqr_accepts_amber_zero_radius_polar_hydrogen(tmp_path):
    from enzymekinetics_prediction.physics.electrostatics import read_pqr

    pqr = tmp_path / "protein.pqr"
    pqr.write_text("ATOM 1 HG1 THR A 1 1.0 2.0 3.0 0.4102 0.0000\n")
    xyz, charge, radius = read_pqr(pqr)
    assert radius[0] == 0
    assert charge[0] == pytest.approx(0.4102)
    np.testing.assert_allclose(xyz, [[1, 2, 3]])


def test_unavailable_apbs_fails_without_silent_backend_substitution(
    tmp_path, monkeypatch
):
    from enzymekinetics_prediction.physics import electrostatics

    monkeypatch.delenv("APBS_EXECUTABLE", raising=False)
    monkeypatch.setattr(electrostatics, "VENDOR_DIR", tmp_path)
    monkeypatch.setattr(electrostatics.shutil, "which", lambda name: None)
    with pytest.raises(RuntimeError, match="APBS executable unavailable"):
        electrostatics.run_apbs(tmp_path / "missing.pqr", tmp_path / "out", [[0, 0, 0]])
