"""Independent structural/PBC references; coordinates only, no MD execution."""

import numpy as np
import pytest

from enzymekinetics_prediction.physics.geometry import (
    make_whole,
    minimum_image,
    orthorhombic_lengths,
    reconstruct_assembly,
    structural_metrics,
)
from enzymekinetics_prediction.physics.dynamics import ensemble_observables


def sample_assembly():
    tetra = np.array([[0, 0, 0], [0.2, 0, 0], [0, 0.3, 0], [0, 0, 0.4]])
    initial = np.vstack((tetra, tetra + [1.4, 0.2, 0.1])) + [2, 2, 2]
    xyz = np.repeat(initial[None], 6, axis=0)
    xyz[:, 3, 2] += np.linspace(0, 0.05, 6)
    xyz[:, 5, 0] -= np.linspace(0, 0.03, 6)
    bonds = [(root, root + i) for root in (0, 4) for i in (1, 2, 3)]
    return xyz, np.repeat((10 * np.eye(3))[None], 6, axis=0), bonds


def test_independent_atom_and_chain_images_preserve_geometry_and_metrics():
    xyz, boxes, bonds = sample_assembly()
    images = np.random.default_rng(101).integers(-4, 5, size=xyz.shape)
    wrapped = xyz + 10 * images
    untouched = wrapped.copy()
    fixed, _, _ = make_whole(wrapped, boxes, bonds)
    np.testing.assert_array_equal(wrapped, untouched)
    for group in (slice(0, 4), slice(4, 8)):
        expected = xyz[:, group, None] - xyz[:, None, group]
        actual = fixed[:, group, None] - fixed[:, None, group]
        np.testing.assert_allclose(actual, expected, atol=3e-14)
    reference = structural_metrics(xyz)
    observed = structural_metrics(fixed)
    for key in ("rmsd_nm", "rmsf_nm", "radius_gyration_nm", "dccm"):
        np.testing.assert_allclose(observed[key], reference[key], atol=1e-10)
    summary, _ = ensemble_observables(fixed)
    assert summary["rmsd_mean_nm"] == pytest.approx(
        reference["rmsd_nm"].mean(), abs=1e-12
    )


def test_kabsch_removes_framewise_rigid_motion_without_reflections():
    xyz, _, _ = sample_assembly()
    moved = xyz.copy()
    for i, theta in enumerate(np.linspace(0, 2, len(xyz))):
        c, s = np.cos(theta), np.sin(theta)
        moved[i] = xyz[i] @ np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]]) + i * 2
    for key in ("rmsd_nm", "rmsf_nm", "dccm"):
        np.testing.assert_allclose(
            structural_metrics(moved)[key], structural_metrics(xyz)[key], atol=1e-10
        )
    mirrored = xyz[0] * [-1, 1, 1]
    metrics = structural_metrics(np.stack((xyz[0], mirrored)))
    assert metrics["rmsd_nm"][1] > 0.01


def test_population_rmsf_and_dccm_have_an_analytic_reference():
    tetra = 0.2 * np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])
    xyz = np.array([tetra, 2 * tetra, 3 * tetra])
    metrics = structural_metrics(xyz)
    np.testing.assert_allclose(
        metrics["rmsd_nm"], np.arange(3) * np.sqrt(0.12), atol=1e-14
    )
    np.testing.assert_allclose(metrics["rmsf_nm"], np.sqrt(0.08), atol=1e-14)
    np.testing.assert_allclose(metrics["dccm"], 4 / 3 * np.eye(4) - 1 / 3, atol=1e-14)
    static = structural_metrics(np.repeat(tetra[None], 3, axis=0))
    np.testing.assert_allclose(static["rmsf_nm"], 0, atol=1e-14)
    np.testing.assert_array_equal(static["dccm"], np.zeros((4, 4)))


@pytest.mark.parametrize(
    "box",
    [np.zeros((3, 3)), np.full((3, 3), np.nan), [[3, 0, 0], [1, 3, 0], [0, 0, 3]]],
)
def test_invalid_and_triclinic_cells_fail_explicitly(box):
    with pytest.raises(ValueError):
        orthorhombic_lengths(box)


def test_per_frame_orthorhombic_cells_and_invalid_groups():
    xyz = np.array([[[0, 0, 0], [11, 0, 0]], [[0, 0, 0], [21, 0, 0]]], float)
    fixed, shifts = reconstruct_assembly(xyz, [[10] * 3, [20] * 3], [[0], [1]])
    np.testing.assert_allclose(fixed[:, 1, 0], [1, 1])
    np.testing.assert_allclose(shifts[:, 1, 0], [-10, -20])
    for groups in ([[0], [0, 1]], [[0]], [[0.0], [1.0]]):
        with pytest.raises(ValueError):
            reconstruct_assembly(xyz, [[10] * 3] * 2, groups)
    with pytest.raises(ValueError):
        minimum_image([np.nan, 0, 0], [10] * 3)


def test_inconsistent_periodic_bond_cycle_is_rejected():
    xyz = np.array([[[0, 0, 0], [4, 0, 0], [8, 0, 0]]], float)
    with pytest.raises(ValueError, match="cycle"):
        make_whole(xyz, [10 * np.eye(3)], [(0, 1), (1, 2), (2, 0)])


def test_half_box_reconstruction_is_explicitly_ambiguous():
    xyz = np.array([[[0, 0, 0], [5, 0, 0]]], float)
    with pytest.raises(ValueError, match="half-box"):
        reconstruct_assembly(xyz, [[10] * 3], [[0], [1]])
    with pytest.raises(ValueError, match="half-box"):
        make_whole(xyz, [10 * np.eye(3)], [(0, 1)])
