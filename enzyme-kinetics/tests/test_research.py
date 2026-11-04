"""Scientific contracts: independent units, geometry, protonation and provenance."""

import json

import numpy as np
import pandas as pd
import pytest

from enzymekinetics_prediction.physics.conditions import (
    sequence_titration,
    water_dielectric,
    ionized_fraction,
    propka_summary,
)
from enzymekinetics_prediction.physics.mechanism import (
    classify_stationary,
    numerical_hessian,
    KSI_MECHANISM,
)
from enzymekinetics_prediction.physics.sites import flatten_assembly
from enzymekinetics_prediction.research.journal import Journal
from enzymekinetics_prediction.research.statistics import (
    paired_sign_flip,
    benjamini_hochberg,
    autocorrelation_diagnostics,
)
from enzymekinetics_prediction.research.torsions import rotate


def test_sequence_pi_is_actual_charge_root_and_charge_decreases_with_ph():
    for sequence in ["DDDD", "KKKK", "ACDEFGHIKLMNPQRSTVWY"]:
        pi, _, _ = sequence_titration(sequence)
        _, charge, _ = sequence_titration(sequence, pi)
        assert abs(charge) < 1e-8
        assert sequence_titration(sequence, 5)[1] > sequence_titration(sequence, 9)[1]
    assert np.isnan(sequence_titration("AX")[0])


def test_environment_temperature_and_titration_units():
    assert water_dielectric(298.15) == pytest.approx(78.30334, abs=1e-4)
    assert water_dielectric(310.15) < water_dielectric(288.15)
    assert ionized_fraction(7, 7) == 0.5
    assert ionized_fraction(9, 7) > ionized_fraction(5, 7)
    with pytest.raises(ValueError):
        water_dielectric(200)


def test_propka_parser_uses_summary_not_determinants(tmp_path):
    p = tmp_path / "pka.log"
    p.write_text(
        "ASP 38 A 999.0 22.0\nSUMMARY OF THIS PREDICTION\n ASP  38 A  6.31  3.80\n ASP  99 A 9.03 3.80\n TYR 14 A 18.39 10.00\nThe pI is 5.83 (folded) and 5.38 (unfolded)\n"
    )
    result = propka_summary(p)["values"]
    assert result["residue_pka"][0]["pka_predicted"] == 6.31
    assert result["folded_pi"] == 5.83


def test_two_step_mechanism_and_stationary_hessian_criteria():
    assert [s["id"] for s in KSI_MECHANISM["states"]] == [
        "E_S",
        "ES",
        "TS1",
        "EI",
        "TS2",
        "EP",
        "E_P",
    ]
    h = numerical_hessian(
        lambda x: np.array([2 * x[0], -4 * x[1]]), np.array([0.0, 0.0])
    )
    np.testing.assert_allclose(h, [[2, 0], [0, -4]])
    assert classify_stationary([0, 0], h)["classification"] == "index_1_saddle"
    assert classify_stationary([1, 0], h)["classification"] == "nonstationary"
    assert classify_stationary([0, 0], np.eye(2))["classification"] == "minimum"
    assert not KSI_MECHANISM["computed_full_TS"]


def test_small_replication_cannot_invent_significance_and_bh():
    test = paired_sign_flip([-1, -2, -3])
    assert test["p_two_sided_exact"] == 0.25
    np.testing.assert_allclose(
        benjamini_hochberg([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04]
    )
    with pytest.raises(ValueError):
        paired_sign_flip([1])


def test_temporal_samples_are_correlated_not_independent():
    rng = np.random.default_rng(42)
    x = np.zeros(2000)
    for i in range(1, len(x)):
        x[i] = 0.95 * x[i - 1] + rng.normal()
    result = autocorrelation_diagnostics(x)
    assert result["effective_samples"] < len(x) / 5
    assert autocorrelation_diagnostics(np.ones(10))["effective_samples"] is None


def test_journal_registration_immutable_and_failures_persist(tmp_path):
    j = Journal(tmp_path)
    j.register({"H0": "x=0"})
    j.register({"H0": "x=0"})
    with pytest.raises(ValueError):
        j.register({"H0": "x=1"})
    with pytest.raises(RuntimeError):
        with j.run("example", {"seed": 42}):
            raise RuntimeError("recorded failure")
    j.result("example", {"value": 1})
    j.result("example", {"value": 2})
    assert len(list((tmp_path / "results").glob("*.json"))) == 2
    events = [
        json.loads(line)
        for line in (tmp_path / "events.jsonl").read_text().splitlines()
    ]
    assert any(r["event"] == "failed" for r in events)


def test_flatten_assembly_preserves_both_models_and_strips_unknown_hydrogens(tmp_path):
    text = "MODEL        1\nATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nATOM      2  H   ALA A   1       0.000   0.000   1.000  1.00  0.00           H\nENDMDL\nMODEL        2\nATOM      1  CA  ALA A   1       2.000   0.000   0.000  1.00  0.00           C\nENDMDL\nEND\n"
    source = tmp_path / "source.pdb"
    source.write_text(text)
    output = tmp_path / "normalized.pdb"
    result = flatten_assembly(source, output)
    assert [r["new_chain"] for r in result["chain_mapping"]] == ["A", "B"]
    assert output.read_text().count("ATOM") == 2


def test_conditional_torsion_rotation_preserves_axis_and_distances():
    xyz = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1]], float)
    rotated = rotate(xyz, (0, 1, [1, 2, 3]), 90)
    np.testing.assert_allclose(rotated[:2], xyz[:2], atol=1e-10)
    np.testing.assert_allclose(
        np.linalg.norm(rotated[2:] - rotated[1], axis=1),
        np.linalg.norm(xyz[2:] - xyz[1], axis=1),
    )


def test_feature_titration_block_is_present_and_legacy_shape_stable(tmp_path):
    from enzymekinetics_prediction.learning.features import FeatureBuilder

    frame = pd.DataFrame(
        {
            "sequence": ["ACDEFGHIK"],
            "uniprot": ["P0"],
            "ec": ["1.1.1.1"],
            "reactant_smiles": ["CCO"],
            "ph": [7.0],
            "temperature": [25.0],
        }
    )
    builder = FeatureBuilder(tmp_path).fit(frame)
    x = builder.transform(frame)
    start, stop = builder.feature_blocks_["titration"]
    np.testing.assert_allclose(
        x[0, start:stop], sequence_titration("ACDEFGHIK", 7), rtol=1e-5
    )
    builder.environment_features_ = False
    assert builder.transform(frame).shape[1] == x.shape[1] - 3


def test_prepared_checkpoint_continues_without_new_minimization(tmp_path):
    mm = pytest.importorskip("openmm")
    pytest.importorskip("mdtraj")
    from enzymekinetics_prediction.physics.dynamics import run_md, MDSettings

    pdb = tmp_path / "a.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00  0.00           C\nATOM      3  CA  ALA B  99       7.600   0.000   0.000  1.00  0.00           C\nEND\n"
    )
    system = mm.System()
    for _ in range(3):
        system.addParticle(12)
    xml = tmp_path / "system.xml"
    xml.write_text(mm.XmlSerializer.serialize(system))
    settings = MDSettings(
        threads=1,
        equilibration_steps=0,
        production_steps=6,
        report_interval=2,
        minimization_iterations=0,
    )
    run_md(pdb, tmp_path / "first", settings, xml)
    result = run_md(
        tmp_path / "first/topology.pdb",
        tmp_path / "second",
        settings,
        xml,
        restart_checkpoint=tmp_path / "first/checkpoint.chk",
    )
    assert result["values"]["frames"] == 3
    assert "restart_checkpoint_sha256" in result["inputs"]
    state = pd.read_csv(tmp_path / "second/state.csv")
    assert state.iloc[-1, 0] == 12
    from Bio.PDB import PDBParser

    structure = PDBParser(QUIET=True).get_structure(
        "s", tmp_path / "second/topology.pdb"
    )
    assert structure[0]["B"][99].resname == "ALA"


def test_focused_apbs_grid_contract_and_resource_limits(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from enzymekinetics_prediction.physics import electrostatics

    executable = tmp_path / "fake-apbs"
    executable.write_text("contract stub, not a solver")
    monkeypatch.setenv("APBS_EXECUTABLE", str(executable))
    pqr = tmp_path / "a.pqr"
    pqr.write_text("ATOM 1 CA ALA A 1 0 0 0 1.0 1.5\n")

    def stub(command, **kwargs):
        assert kwargs["env"]["OMP_NUM_THREADS"] == "2"
        out = Path(kwargs["cwd"])
        (out / "potential.dx").write_text(
            "object 1 class gridpositions counts 2 2 2\norigin 0 0 0\ndelta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\nobject 3 class array type double rank 0 items 8 data follows\n0 1 1 2 1 2 2 3\n"
        )
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    from pathlib import Path

    monkeypatch.setattr(electrostatics.subprocess, "run", stub)
    result = electrostatics.run_apbs(
        pqr,
        tmp_path / "run",
        [[0.5, 0.5, 0.5]],
        threads=2,
        grid_points=129,
        focus_length_A=20,
    )
    assert result["settings"]["fine_lengths_A"] == [20, 20, 20]
    assert result["settings"]["fine_center_A"] == [0.5, 0.5, 0.5]
    assert result["settings"]["grid"] == [129] * 3
    assert "dime 129 129 129" in (tmp_path / "run/apbs.in").read_text()
