from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import analyze_spectral_robustness as robust


def test_gls_offset_fit_recovers_known_offset():
    values = np.asarray([2.0, 3.0, 4.0])
    model = np.asarray([1.5, 2.5, 3.5])
    covariance = np.diag([0.1, 0.2, 0.3]) ** 2
    fit = robust.gls_offset_fit(values, model, covariance)
    assert np.isclose(fit["offset"], 0.5)
    assert fit["chi2"] < 1e-20


def test_scale_height_is_physical():
    result = robust.atmospheric_scale_height()
    assert 1 < result["gravity_m_s2"] < 30
    assert 100 < result["scale_height_km"] < 5000
    assert 10 < result["one_scale_height_ppm"] < 1000


def test_real_spectrum_robustness_outputs_are_stable():
    result = robust.main()
    delta = result["primary"]["delta_chi2_no_co2_minus_full"]
    assert result["inputs"]["archive_bins"][0] == 94
    assert len(result["inputs"]["wavelength"]) == 93
    assert result["inputs"]["excluded_bins"][0] == 1
    assert delta > 100
    assert result["leave_one_out"].min() > 100
    multiverse = np.asarray([
        float(row["delta_chi2_no_co2_minus_full"])
        for row in result["multiverse_rows"]
    ])
    assert len(multiverse) == 60
    assert np.all(multiverse > 0)
    assert len(result["covariance_rows"]) == 20
    assert min(float(row["delta_chi2_no_co2_minus_full"])
               for row in result["block_rows"]) > 0
    assert result["region_results"]["inside"]["delta_chi2_no_co2_minus_full"] > 0
    for path in (
        robust.SUMMARY_FILE, robust.MULTIVERSE_FILE, robust.COVARIANCE_FILE,
        robust.BLOCK_FILE, robust.CONTRIBUTIONS_FILE, robust.FIGURE_FILE,
    ):
        assert Path(path).is_file() and Path(path).stat().st_size > 200


def test_bin_integral_matches_linear_model_average():
    model_wavelength = np.asarray([0.0, 0.5, 1.0, 1.5, 2.0])
    model_value = 3.0 * model_wavelength + 2.0
    left = np.asarray([0.2, 1.1])
    right = np.asarray([0.8, 1.9])
    expected = 3.0 * (left + right) / 2.0 + 2.0
    assert np.allclose(
        robust._bin_average(left, right, model_wavelength, model_value), expected
    )


def test_rebinning_rejects_invalid_origins():
    with np.testing.assert_raises(ValueError):
        robust.rebin_arrays(robust.load_inputs(), factor=2, origin=2)
