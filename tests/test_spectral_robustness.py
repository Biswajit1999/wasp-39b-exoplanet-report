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
    delta = result["base"]["delta_chi2_no_co2_minus_full"]
    assert len(result["inputs"]["wavelength"]) == 94
    assert delta > 100
    assert result["leave_one_out"].min() > 100
    assert all(row["delta_chi2_no_co2_minus_full"] > 20 for row in result["binning_rows"])
    assert result["amplitude_scale_heights"] > 1
    for path in (robust.SUMMARY_FILE, robust.CONTRIBUTIONS_FILE, robust.FIGURE_FILE):
        assert Path(path).is_file() and Path(path).stat().st_size > 200
