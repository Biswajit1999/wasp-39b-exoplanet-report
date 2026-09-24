"""Reproduction checks for the archived published spectrum."""
from pathlib import Path
import numpy as np
import analyze_spectrum as spectrum

def test_spectrum_analysis_is_finite_and_reproducible():
    result = spectrum.main()
    assert result["n"] == 94
    assert result["common_support_n"] == 93
    assert result["excluded_model_bins"] == 1
    assert result["rows"]
    assert all(np.isfinite(row["chi2"]) and row["dof"] > 0 for row in result["rows"])
    for path in (spectrum.STATS_FILE, spectrum.FIGURE_FILE):
        assert Path(path).is_file() and Path(path).stat().st_size > 100


def test_model_comparison_rejects_extrapolation():
    with np.testing.assert_raises(ValueError):
        spectrum.offset_model_test(
            np.asarray([0.0, 1.0]), np.asarray([1.0, 1.0]), np.asarray([0.1, 0.1]),
            np.asarray([0.2, 0.8]), np.asarray([1.0, 1.0]),
        )
