"""Sensitivity audit for the archived WASP-39 b PRISM products.

The output is a fixed-forward-model diagnostic, not a retrieval or molecular
detection significance. Every comparison uses shared support and the same
single fitted vertical offset for the full and remove-CO2 archive products.
"""
from __future__ import annotations

import csv
from pathlib import Path
from astropy import units as u
from astropy.constants import G, M_earth, R_earth, R_sun, k_B, m_p
from astropy.table import Table
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapezoid

ROOT = Path(__file__).resolve().parents[1]
DATA, FIGURES = ROOT / "data" / "spectra", ROOT / "figures"
SUMMARY_FILE = FIGURES / "spectral_robustness_statistics.csv"
MULTIVERSE_FILE = FIGURES / "spectral_multiverse.csv"
COVARIANCE_FILE = FIGURES / "spectral_covariance_grid.csv"
BLOCK_FILE = FIGURES / "spectral_block_deletions.csv"
CONTRIBUTIONS_FILE = FIGURES / "spectral_bin_contributions.csv"
FIGURE_FILE = FIGURES / "wasp39b_spectral_robustness.png"
REBIN_FACTORS = (1, 2, 4, 8)
ERROR_RULES = ("mean", "maximum")
SAMPLING_RULES = ("center", "bin_integrated")
COVARIANCE_FRACTIONS = (0.0, 0.25, 0.5, 1.0, 2.0)
COVARIANCE_LENGTHS = (0.05, 0.10, 0.25, 0.50)
BLOCK_WIDTHS = (0.10, 0.25, 0.50)


def _bin_average(left, right, model_wavelength, model_value):
    """Integrate a piecewise-linear model over each observed bin."""
    result = []
    for lower, upper in zip(left, right):
        inside = model_wavelength[(model_wavelength > lower) & (model_wavelength < upper)]
        grid = np.concatenate(([lower], inside, [upper]))
        result.append(trapezoid(np.interp(grid, model_wavelength, model_value), grid) / (upper - lower))
    return np.asarray(result)


def load_inputs(model_sampling="center", error_rule="mean"):
    if model_sampling not in SAMPLING_RULES or error_rule not in ERROR_RULES:
        raise ValueError("unknown sensitivity-design rule")
    table = Table.read(DATA / "eureka_transmission_spectrum.txt", format="ascii.ecsv")
    wavelength = np.asarray(table["wavelength"], float)
    width = np.asarray(table["bin_width"], float)
    left, right = wavelength - width / 2, wavelength + width / 2
    depth = np.asarray(table["tr_depth"], float) * 1e6
    neg = np.asarray(table["tr_depth_errneg"], float) * 1e6
    pos = np.asarray(table["tr_depth_errpos"], float) * 1e6
    error = (neg + pos) / 2 if error_rule == "mean" else np.maximum(neg, pos)
    full_raw = np.loadtxt(DATA / "scchimeramodel.txt")
    no_raw = np.loadtxt(DATA / "scchimeramodel_no_co2.txt")
    support = np.asarray([max(full_raw[:, 0].min(), no_raw[:, 0].min()),
                          min(full_raw[:, 0].max(), no_raw[:, 0].max())])
    keep = (left >= support[0]) & (right <= support[1])
    wavelength, width, left, right, depth, error = (
        x[keep] for x in (wavelength, width, left, right, depth, error)
    )
    def sample(raw):
        values = raw[:, 1] * 1e6
        if model_sampling == "center":
            return np.interp(wavelength, raw[:, 0], values)
        return _bin_average(left, right, raw[:, 0], values)
    return {"wavelength": wavelength, "bin_width": width, "depth": depth,
            "error": error, "full": sample(full_raw), "no_co2": sample(no_raw),
            "archive_bins": np.asarray([len(table)]),
            "excluded_bins": np.asarray([int((~keep).sum())]), "support": support}


def gls_offset_fit(values, model, covariance):
    precision = np.linalg.pinv(covariance, hermitian=True)
    ones = np.ones(len(values))
    offset = float((ones @ precision @ (values - model)) / (ones @ precision @ ones))
    fitted = model + offset
    residual = values - fitted
    return {"offset": offset, "model": fitted, "residual": residual,
            "chi2": float(residual @ precision @ residual)}


def covariance_matrix(wavelength, error, correlated_fraction, correlation_length_micron):
    covariance = np.diag(error ** 2)
    if correlated_fraction > 0:
        amplitude = correlated_fraction * float(np.median(error))
        separation = np.abs(wavelength[:, None] - wavelength[None, :])
        covariance += amplitude ** 2 * np.exp(-separation / correlation_length_micron)
    return covariance


def compare_models(wavelength, depth, error, full, no_co2,
                   correlated_fraction=0.0, correlation_length_micron=0.1, **_):
    covariance = covariance_matrix(wavelength, error, correlated_fraction,
                                   correlation_length_micron)
    full_fit = gls_offset_fit(depth, full, covariance)
    no_fit = gls_offset_fit(depth, no_co2, covariance)
    eig = np.linalg.eigvalsh(covariance)
    return {"full": full_fit, "no_co2": no_fit,
            "delta_chi2_no_co2_minus_full": no_fit["chi2"] - full_fit["chi2"],
            "covariance_condition_number": float(eig[-1] / eig[0]),
            "covariance_min_eigenvalue": float(eig[0])}


def rebin_arrays(inputs, factor, origin=0):
    if factor < 1 or not 0 <= origin < factor:
        raise ValueError("require factor >= 1 and 0 <= origin < factor")
    keys = ("wavelength", "bin_width", "depth", "error", "full", "no_co2")
    order = np.argsort(inputs["wavelength"])
    source = {key: np.asarray(inputs[key])[order] for key in keys}
    result = {key: [] for key in keys}
    for start in range(origin, len(order) - factor + 1, factor):
        selected = slice(start, start + factor)
        weights = 1 / source["error"][selected] ** 2
        for key in ("wavelength", "depth", "full", "no_co2"):
            result[key].append(float(np.average(source[key][selected], weights=weights)))
        result["error"].append(float(np.sqrt(1 / weights.sum())))
        lower = source["wavelength"][selected] - source["bin_width"][selected] / 2
        upper = source["wavelength"][selected] + source["bin_width"][selected] / 2
        result["bin_width"].append(float(upper.max() - lower.min()))
    return {key: np.asarray(value) for key, value in result.items()}


def atmospheric_scale_height():
    with (ROOT / "data" / "system_parameters.csv").open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle))
    mass, radius = float(row["pl_bmasse"]) * M_earth, float(row["pl_rade"]) * R_earth
    star_radius, temperature = float(row["st_rad"]) * R_sun, float(row["pl_eqt"])
    gravity = G * mass / radius ** 2
    height = k_B * (temperature * u.K) / (2.3 * m_p * gravity)
    return {"gravity_m_s2": gravity.to_value("m/s2"), "scale_height_km": height.to_value("km"),
            "one_scale_height_ppm": float((2 * radius * height / star_radius ** 2 * 1e6).decompose().value)}


def _write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def _subset(inputs, keep):
    n = len(inputs["wavelength"])
    return {key: value[keep] for key, value in inputs.items() if len(value) == n}


def main():
    FIGURES.mkdir(exist_ok=True)
    inputs = load_inputs()
    base = compare_models(**inputs)
    wavelength = inputs["wavelength"]
    full_fit, no_fit = base["full"], base["no_co2"]
    contribution = ((no_fit["residual"] / inputs["error"]) ** 2
                    - (full_fit["residual"] / inputs["error"]) ** 2)
    loo = np.asarray([compare_models(**_subset(inputs, np.arange(len(wavelength)) != i))
                      ["delta_chi2_no_co2_minus_full"] for i in range(len(wavelength))])

    multiverse = []
    for sampling in SAMPLING_RULES:
        for error_rule in ERROR_RULES:
            design = load_inputs(sampling, error_rule)
            for factor in REBIN_FACTORS:
                for origin in range(factor):
                    rebinned = rebin_arrays(design, factor, origin)
                    result = compare_models(**rebinned)
                    multiverse.append({"design_id": len(multiverse) + 1,
                        "model_sampling": sampling, "asymmetric_error_rule": error_rule,
                        "rebin_factor": factor, "rebin_origin": origin,
                        "n_bins": len(rebinned["wavelength"]),
                        "full_chi2": f"{result['full']['chi2']:.12g}",
                        "no_co2_chi2": f"{result['no_co2']['chi2']:.12g}",
                        "delta_chi2_no_co2_minus_full": f"{result['delta_chi2_no_co2_minus_full']:.12g}"})

    covariance = []
    for fraction in COVARIANCE_FRACTIONS:
        for length in COVARIANCE_LENGTHS:
            result = compare_models(**inputs, correlated_fraction=fraction,
                                    correlation_length_micron=length)
            covariance.append({"correlated_amplitude_over_median_error": fraction,
                "correlation_length_micron": length, "n_bins": len(wavelength),
                "covariance_condition_number": f"{result['covariance_condition_number']:.12g}",
                "covariance_min_eigenvalue": f"{result['covariance_min_eigenvalue']:.12g}",
                "full_chi2": f"{result['full']['chi2']:.12g}",
                "no_co2_chi2": f"{result['no_co2']['chi2']:.12g}",
                "delta_chi2_no_co2_minus_full": f"{result['delta_chi2_no_co2_minus_full']:.12g}"})

    blocks = []
    for width in BLOCK_WIDTHS:
        seen = set()
        for center in wavelength:
            deleted = np.abs(wavelength - center) <= width / 2
            signature = tuple(np.flatnonzero(deleted))
            if not signature or signature in seen: continue
            seen.add(signature)
            result = compare_models(**_subset(inputs, ~deleted))
            blocks.append({"requested_width_micron": width, "center_micron": f"{center:.10g}",
                "deleted_min_micron": f"{wavelength[deleted].min():.10g}",
                "deleted_max_micron": f"{wavelength[deleted].max():.10g}",
                "n_deleted": int(deleted.sum()), "n_retained": int((~deleted).sum()),
                "delta_chi2_no_co2_minus_full": f"{result['delta_chi2_no_co2_minus_full']:.12g}"})

    window = (wavelength >= 4.1) & (wavelength <= 4.6)
    regions = {name: compare_models(**_subset(inputs, keep)) for name, keep in
               (("inside", window), ("outside", ~window))}
    scale = atmospheric_scale_height()
    amplitude = float(np.percentile(inputs["depth"], 95) - np.percentile(inputs["depth"], 5))
    multi_delta = np.asarray([float(row["delta_chi2_no_co2_minus_full"]) for row in multiverse])
    cov_delta = np.asarray([float(row["delta_chi2_no_co2_minus_full"]) for row in covariance])
    block_delta = np.asarray([float(row["delta_chi2_no_co2_minus_full"]) for row in blocks])
    summary = [
        ("archive_spectral_bins", int(inputs["archive_bins"][0]), "count", "all Eureka bins"),
        ("common_support_bins", len(wavelength), "count", "entire bin supported by both models"),
        ("excluded_extrapolation_bins", int(inputs["excluded_bins"][0]), "count", "excluded before comparison"),
        ("primary_delta_chi2", base["delta_chi2_no_co2_minus_full"], "descriptive contrast", "center sampling; mean asymmetric error"),
        ("multiverse_designs", len(multiverse), "count", "sampling x error x factor x every origin"),
        ("multiverse_positive_fraction", np.mean(multi_delta > 0), "fraction", "preference for supplied full model"),
        ("multiverse_delta_min", multi_delta.min(), "descriptive contrast", "declared multiverse"),
        ("multiverse_delta_median", np.median(multi_delta), "descriptive contrast", "declared multiverse"),
        ("multiverse_delta_max", multi_delta.max(), "descriptive contrast", "declared multiverse"),
        ("covariance_grid_cells", len(covariance), "count", "amplitude x length grid"),
        ("covariance_delta_min", cov_delta.min(), "descriptive contrast", "declared covariance grid"),
        ("leave_one_bin_delta_min", loo.min(), "descriptive contrast", "all single-bin deletions"),
        ("block_deletion_delta_min", block_delta.min(), "descriptive contrast", "all contiguous block deletions"),
        ("inside_4p1_4p6_delta", regions["inside"]["delta_chi2_no_co2_minus_full"], "descriptive contrast", "models refit inside interval"),
        ("outside_4p1_4p6_delta", regions["outside"]["delta_chi2_no_co2_minus_full"], "descriptive contrast", "models refit outside interval"),
        ("surface_gravity", scale["gravity_m_s2"], "m s-2", "saved system parameters"),
        ("atmospheric_scale_height", scale["scale_height_km"], "km", "illustrative mu=2.3"),
        ("one_scale_height_signal", scale["one_scale_height_ppm"], "ppm", "illustrative"),
        ("robust_spectral_amplitude", amplitude, "ppm", "95th minus 5th percentile"),
        ("robust_amplitude_scale_heights", amplitude / scale["one_scale_height_ppm"], "illustrative scale heights", "not composition"),
    ]
    _write(SUMMARY_FILE, [{"quantity": q, "value": f"{v:.12g}" if isinstance(v, (float, np.floating)) else v,
                           "unit": unit, "interpretation": note} for q, v, unit, note in summary])
    _write(MULTIVERSE_FILE, multiverse); _write(COVARIANCE_FILE, covariance); _write(BLOCK_FILE, blocks)
    _write(CONTRIBUTIONS_FILE, [{"wavelength_micron": f"{w:.10g}",
        "bin_width_micron": f"{inputs['bin_width'][i]:.10g}", "depth_ppm": f"{inputs['depth'][i]:.10g}",
        "error_ppm": f"{inputs['error'][i]:.10g}", "full_model_with_offset_ppm": f"{full_fit['model'][i]:.10g}",
        "no_co2_model_with_offset_ppm": f"{no_fit['model'][i]:.10g}", "delta_chi2_contribution": f"{contribution[i]:.10g}",
        "leave_one_out_delta_chi2": f"{loo[i]:.10g}", "in_4p1_4p6_window": bool(window[i])}
        for i, w in enumerate(wavelength)])

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 9), constrained_layout=True)
    ax = axes[0, 0]; ax.errorbar(wavelength, inputs["depth"], yerr=inputs["error"], fmt="o", ms=3.1,
        color="#17212b", ecolor="#94a3b8", label="Eureka spectrum")
    ax.plot(wavelength, full_fit["model"], color="#0f766e", lw=2.1, label="full + offset")
    ax.plot(wavelength, no_fit["model"], color="#be123c", lw=1.8, label="remove-CO2 + offset")
    ax.axvspan(4.1, 4.6, color="#d97706", alpha=.12, label="4.1–4.6 μm")
    ax.set(xlabel="Wavelength [μm]", ylabel="Transit depth [ppm]", title="Strict common-support comparison"); ax.grid(alpha=.2); ax.legend(frameon=False, fontsize=8)
    ax = axes[0, 1]; ax.bar(wavelength, contribution, width=.035, color=np.where(window, "#d97706", "#64748b")); ax.axhline(0, color="#334155", lw=1)
    ax.set(xlabel="Wavelength [μm]", ylabel="Per-bin Δχ²", title="Primary-design contribution map"); ax.grid(axis="y", alpha=.2)
    ax = axes[1, 0]
    for factor, color in zip(REBIN_FACTORS, ("#0f766e", "#2563eb", "#7c3aed", "#d97706")):
        values = [float(r["delta_chi2_no_co2_minus_full"]) for r in multiverse if r["rebin_factor"] == factor]
        ax.scatter(np.full(len(values), factor) + np.linspace(-.09, .09, len(values)), values, s=24, alpha=.78, color=color)
    ax.axhline(0, color="#334155", lw=1); ax.set(xticks=REBIN_FACTORS, xlabel="Adjacent-bin factor (all origins)", ylabel="Δχ² (remove-CO2 − full)", title=f"Declared multiverse · {len(multiverse)} designs"); ax.grid(axis="y", alpha=.2)
    ax = axes[1, 1]
    for fraction, color in zip(COVARIANCE_FRACTIONS, ("#0f766e", "#0891b2", "#2563eb", "#7c3aed", "#be123c")):
        rows = [r for r in covariance if r["correlated_amplitude_over_median_error"] == fraction]
        ax.plot([r["correlation_length_micron"] for r in rows], [float(r["delta_chi2_no_co2_minus_full"]) for r in rows], "o-", color=color, label=f"amplitude {fraction:g}×")
    ax.axhline(0, color="#334155", lw=1); ax.set(xlabel="Correlation length [μm]", ylabel="Δχ² (remove-CO2 − full)", title="Illustrative covariance sensitivity"); ax.grid(alpha=.2); ax.legend(frameon=False, fontsize=8)
    fig.suptitle("WASP-39 b · supplied-model sensitivity audit", fontsize=15, weight="bold")
    fig.savefig(FIGURE_FILE, dpi=190); plt.close(fig)
    return {"inputs": inputs, "primary": base, "base": base, "leave_one_out": loo,
            "multiverse_rows": multiverse, "covariance_rows": covariance,
            "block_rows": blocks, "region_results": regions, "scale": scale}


if __name__ == "__main__":
    result = main()
    delta = np.asarray([float(r["delta_chi2_no_co2_minus_full"]) for r in result["multiverse_rows"]])
    print(f"WASP-39 b: {len(result['inputs']['wavelength'])} common-support bins; "
          f"primary delta-chi2={result['primary']['delta_chi2_no_co2_minus_full']:.1f}; "
          f"{len(delta)}-design range={delta.min():.1f}-{delta.max():.1f}")
