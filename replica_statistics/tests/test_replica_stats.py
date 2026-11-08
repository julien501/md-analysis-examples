"""Checks of replica_stats against scipy and against known cases. Run with: python -m pytest tests"""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import replica_stats as rs  # noqa: E402


def test_anova_matches_scipy():
    rng = np.random.default_rng(0)
    groups = np.repeat(["WT", "A", "B", "C"], [3, 3, 4, 2])          # unbalanced on purpose
    values = rng.normal(0, 1, len(groups)) + (groups == "B") * 2
    a = rs.variant_anova(values, groups)
    F, p = stats.f_oneway(*[values[groups == g] for g in ["WT", "A", "B", "C"]])
    assert a["F"] == pytest.approx(F)
    assert a["p"] == pytest.approx(p)


def test_icc_limits():
    groups = np.repeat(["WT", "A", "B"], 3)
    no_effect = np.array([1.0, 2.0, 3.0] * 3)                       # identical systems
    assert rs.variant_anova(no_effect, groups)["icc"] < 0.0           # read as 0
    all_effect = np.repeat([1.0, 5.0, 9.0], 3) + np.tile([0.0, 0.01, -0.01], 3)
    assert rs.variant_anova(all_effect, groups)["icc"] > 0.99


def test_constant_descriptor_is_not_significant():
    groups = np.repeat(["WT", "A", "B"], 3)
    a = rs.variant_anova(np.full(9, 0.3) + 1e-17 * np.arange(9), groups)
    assert np.isnan(a["p"])


def test_dunnett_only_after_significant_anova():
    rng = np.random.default_rng(2)
    groups = np.repeat(["WT", "A", "B"], 3)
    _, d = rs.variant_test(rng.normal(0, 1, 9), groups)
    assert d.empty
    _, d = rs.variant_test(rng.normal(0, 0.1, 9) + (groups == "B") * 5, groups)
    assert d.set_index("variant").loc["B", "p"] < 0.05


def test_free_energy_zero_at_maximum():
    g = rs.free_energy([0.1, 0.5, 0.2])
    assert g[1] == 0.0 and (g >= 0).all()


def test_replicas_needed_textbook_value():
    # sd = delta gives 2 * (1.96 + 0.84)^2 = 15.7 -> 16 per group
    assert rs.replicas_needed(1.0, 1.0) == 16


def test_bh_matches_scipy():
    p = np.array([0.001, 0.02, 0.03, 0.2, np.nan, 0.8])
    q = rs.bh(p)
    ok = np.isfinite(p)
    assert np.allclose(q[ok], stats.false_discovery_control(p[ok], method="bh"))
    assert np.isnan(q[~ok]).all()
