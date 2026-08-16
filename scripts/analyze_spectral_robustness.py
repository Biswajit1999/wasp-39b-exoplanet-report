"""Robustness tests for the public WASP-39 b NIRSpec/PRISM spectrum.

The calculation compares the archived ScCHIMERA full and no-CO2 model files
under identical vertical-offset treatment.  It stress-tests binning, single-bin
deletion, and illustrative wavelength-correlated covariance.  These are model
comparison diagnostics, not an atmospheric retrieval or a conversion to a
molecular detection significance.
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

import analyze_spectrum as spectrum


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "spectra"
FIGURES = ROOT / "figures"
SUMMARY_FILE = FIGURES / "spectral_robustness_statistics.csv"
CONTRIBUTIONS_FILE = FIGURES / "spectral_bin_contributions.csv"
FIGURE_FILE = FIGURES / "wasp39b_spectral_robustness.png"


def load_inputs() -> dict[str, np.ndarray]:
    table = Table.read(DATA / "eureka_transmission_spectrum.txt", format="ascii.ecsv")
    wavelength = np.asarray(table["wavelength"], dtype=float)
    depth = np.asarray(table["tr_depth"], dtype=float) * 1e6
    error = 0.5 * (
        np.asarray(table["tr_depth_errneg"], dtype=float)
        + np.asarray(table["tr_depth_errpos"], dtype=float)
    ) * 1e6
    full_raw = np.loadtxt(DATA / "scchimeramodel.txt")
    no_co2_raw = np.loadtxt(DATA / "scchimeramodel_no_co2.txt")
    full = np.interp(wavelength, full_raw[:, 0], full_raw[:, 1] * 1e6)
    no_co2 = np.interp(wavelength, no_co2_raw[:, 0], no_co2_raw[:, 1] * 1e6)
    return {
        "wavelength": wavelength,
        "depth": depth,
        "error": error,
        "full": full,
        "no_co2": no_co2,
    }


def gls_offset_fit(
    values: np.ndarray,
    model: np.ndarray,
    covariance: np.ndarray,
) -> dict[str, object]:
    precision = np.linalg.pinv(covariance)
    ones = np.ones(len(values))
    offset = float((ones @ precision @ (values - model)) / (ones @ precision @ ones))
    fitted = model + offset
    residual = values - fitted
    chi_square = float(residual @ precision @ residual)
    return {"offset": offset, "model": fitted, "residual": residual, "chi2": chi_square}


def covariance_matrix(
    wavelength: np.ndarray,
    error: np.ndarray,
    correlated_fraction: float,
    correlation_length_micron: float,
) -> np.ndarray:
    covariance = np.diag(error**2)
    if correlated_fraction <= 0:
        return covariance
    amplitude = correlated_fraction * float(np.median(error))
    separation = np.abs(wavelength[:, None] - wavelength[None, :])
    covariance += amplitude**2 * np.exp(-separation / correlation_length_micron)
    return covariance


def compare_models(
    wavelength: np.ndarray,
    depth: np.ndarray,
    error: np.ndarray,
    full: np.ndarray,
    no_co2: np.ndarray,
    correlated_fraction: float = 0.0,
    correlation_length_micron: float = 0.1,
) -> dict[str, object]:
    covariance = covariance_matrix(
        wavelength, error, correlated_fraction, correlation_length_micron
    )
    full_fit = gls_offset_fit(depth, full, covariance)
    no_fit = gls_offset_fit(depth, no_co2, covariance)
    return {
        "full": full_fit,
        "no_co2": no_fit,
        "delta_chi2_no_co2_minus_full": no_fit["chi2"] - full_fit["chi2"],
    }


def rebin_arrays(inputs: dict[str, np.ndarray], factor: int) -> dict[str, np.ndarray]:
    order = np.argsort(inputs["wavelength"])
    sorted_inputs = {key: np.asarray(value)[order] for key, value in inputs.items()}
    result = {key: [] for key in sorted_inputs}
    count = len(order)
    for start in range(0, count, factor):
        stop = min(start + factor, count)
        if stop - start < max(1, factor // 2):
            continue
        selected = slice(start, stop)
        weights = 1.0 / sorted_inputs["error"][selected] ** 2
        result["wavelength"].append(float(np.average(sorted_inputs["wavelength"][selected], weights=weights)))
        result["depth"].append(float(np.average(sorted_inputs["depth"][selected], weights=weights)))
        result["error"].append(float(np.sqrt(1.0 / np.sum(weights))))
        result["full"].append(float(np.average(sorted_inputs["full"][selected], weights=weights)))
        result["no_co2"].append(float(np.average(sorted_inputs["no_co2"][selected], weights=weights)))
    return {key: np.asarray(value, dtype=float) for key, value in result.items()}


def atmospheric_scale_height() -> dict[str, float]:
    with (ROOT / "data" / "system_parameters.csv").open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle))
    mass = float(row["pl_bmasse"]) * M_earth
    radius = float(row["pl_rade"]) * R_earth
    star_radius = float(row["st_rad"]) * R_sun
    temperature = float(row["pl_eqt"])
    mean_molecular_weight = 2.3
    gravity_quantity = G * mass / radius**2
    gravity = gravity_quantity.to_value("m/s2")
    height_quantity = k_B * (temperature * u.K) / (
        mean_molecular_weight * m_p * gravity_quantity
    )
    height = height_quantity.to_value("km")
    one_height_ppm = (2 * radius * height_quantity / star_radius**2 * 1e6).decompose().value
    return {
        "gravity_m_s2": gravity,
        "temperature_k": temperature,
        "mean_molecular_weight": mean_molecular_weight,
        "scale_height_km": height,
        "one_scale_height_ppm": float(one_height_ppm),
    }


def main() -> dict[str, object]:
    FIGURES.mkdir(exist_ok=True)
    inputs = load_inputs()
    base = compare_models(**inputs)
    full_fit = base["full"]
    no_fit = base["no_co2"]
    per_bin_delta = (
        (np.asarray(no_fit["residual"]) / inputs["error"]) ** 2
        - (np.asarray(full_fit["residual"]) / inputs["error"]) ** 2
    )

    leave_one_out = []
    for index in range(len(inputs["wavelength"])):
        keep = np.arange(len(inputs["wavelength"])) != index
        result = compare_models(**{key: value[keep] for key, value in inputs.items()})
        leave_one_out.append(float(result["delta_chi2_no_co2_minus_full"]))
    leave_one_out = np.asarray(leave_one_out)

    binning_rows = []
    for factor in (1, 2, 4, 8):
        rebinned = rebin_arrays(inputs, factor)
        result = compare_models(**rebinned)
        binning_rows.append({
            "test": "adjacent-bin rebinning",
            "setting": f"factor {factor}",
            "n_bins": len(rebinned["wavelength"]),
            "correlated_fraction": 0.0,
            "correlation_length_micron": 0.0,
            "full_chi2": result["full"]["chi2"],
            "no_co2_chi2": result["no_co2"]["chi2"],
            "delta_chi2_no_co2_minus_full": result["delta_chi2_no_co2_minus_full"],
        })

    covariance_rows = []
    for fraction, length in ((0.0, 0.1), (0.5, 0.1), (1.0, 0.1), (1.0, 0.25)):
        result = compare_models(
            **inputs,
            correlated_fraction=fraction,
            correlation_length_micron=length,
        )
        covariance_rows.append({
            "test": "illustrative correlated covariance",
            "setting": f"amplitude {fraction:.1f} x median error; length {length:.2f} micron",
            "n_bins": len(inputs["wavelength"]),
            "correlated_fraction": fraction,
            "correlation_length_micron": length,
            "full_chi2": result["full"]["chi2"],
            "no_co2_chi2": result["no_co2"]["chi2"],
            "delta_chi2_no_co2_minus_full": result["delta_chi2_no_co2_minus_full"],
        })

    scale = atmospheric_scale_height()
    robust_amplitude = float(np.percentile(inputs["depth"], 95) - np.percentile(inputs["depth"], 5))
    amplitude_scale_heights = robust_amplitude / scale["one_scale_height_ppm"]
    co2_window = (inputs["wavelength"] >= 4.1) & (inputs["wavelength"] <= 4.6)
    co2_window_delta = float(np.sum(per_bin_delta[co2_window]))

    summary_rows = [
        ("spectral_bins", len(inputs["wavelength"]), "count"),
        ("wavelength_min", float(inputs["wavelength"].min()), "micron"),
        ("wavelength_max", float(inputs["wavelength"].max()), "micron"),
        ("full_model_chi_square", full_fit["chi2"], "one fitted offset"),
        ("no_co2_model_chi_square", no_fit["chi2"], "one fitted offset"),
        ("delta_chi2_no_co2_minus_full", base["delta_chi2_no_co2_minus_full"], "diagnostic; not retrieval sigma"),
        ("co2_window_delta_chi2", co2_window_delta, "4.1-4.6 micron contribution"),
        ("leave_one_bin_delta_chi2_min", float(leave_one_out.min()), "diagnostic"),
        ("leave_one_bin_delta_chi2_max", float(leave_one_out.max()), "diagnostic"),
        ("robust_spectral_amplitude_ppm", robust_amplitude, "95th minus 5th percentile"),
        ("surface_gravity", scale["gravity_m_s2"], "m s-2"),
        ("assumed_mean_molecular_weight", scale["mean_molecular_weight"], "proton masses"),
        ("atmospheric_scale_height", scale["scale_height_km"], "km; using equilibrium temperature"),
        ("one_scale_height_transmission_amplitude", scale["one_scale_height_ppm"], "ppm"),
        ("robust_amplitude_in_scale_heights", amplitude_scale_heights, "illustrative H2-He atmosphere"),
    ]
    with SUMMARY_FILE.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value", "unit"])
        for name, value, unit in summary_rows:
            writer.writerow([name, f"{value:.12g}" if isinstance(value, float) else value, unit])
        writer.writerow([])
        fields = [
            "test", "setting", "n_bins", "correlated_fraction",
            "correlation_length_micron", "full_chi2", "no_co2_chi2",
            "delta_chi2_no_co2_minus_full",
        ]
        dict_writer = csv.DictWriter(handle, fieldnames=fields)
        dict_writer.writeheader()
        dict_writer.writerows(binning_rows + covariance_rows)

    with CONTRIBUTIONS_FILE.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "wavelength_micron", "depth_ppm", "error_ppm", "full_model_ppm",
            "no_co2_model_ppm", "delta_chi2_contribution",
            "leave_one_out_delta_chi2", "in_co2_4p1_4p6_window",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for index, wavelength in enumerate(inputs["wavelength"]):
            writer.writerow({
                "wavelength_micron": f"{wavelength:.10g}",
                "depth_ppm": f"{inputs['depth'][index]:.10g}",
                "error_ppm": f"{inputs['error'][index]:.10g}",
                "full_model_ppm": f"{np.asarray(full_fit['model'])[index]:.10g}",
                "no_co2_model_ppm": f"{np.asarray(no_fit['model'])[index]:.10g}",
                "delta_chi2_contribution": f"{per_bin_delta[index]:.10g}",
                "leave_one_out_delta_chi2": f"{leave_one_out[index]:.10g}",
                "in_co2_4p1_4p6_window": bool(co2_window[index]),
            })

    fig, axes = plt.subplots(2, 2, figsize=(12, 8.7), constrained_layout=True)
    ax = axes[0, 0]
    ax.errorbar(
        inputs["wavelength"], inputs["depth"], yerr=inputs["error"],
        fmt="o", ms=3.2, color="#17212b", ecolor="#94a3b8", label="public Eureka spectrum",
    )
    ax.plot(inputs["wavelength"], full_fit["model"], color="#0f766e", lw=2.1, label="full model + offset")
    ax.plot(inputs["wavelength"], no_fit["model"], color="#be123c", lw=1.8, label="no-CO2 model + offset")
    ax.axvspan(4.1, 4.6, color="#f59e0b", alpha=0.12, label="4.1-4.6 micron window")
    ax.set(xlabel="Wavelength [micron]", ylabel="Transit depth [ppm]", title="Identical-offset model comparison")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    colors = np.where(co2_window, "#d97706", "#64748b")
    ax.bar(inputs["wavelength"], per_bin_delta, width=0.035, color=colors, alpha=0.85)
    ax.axhline(0, color="#334155", lw=1)
    ax.set(
        xlabel="Wavelength [micron]", ylabel="Per-bin Delta chi-square",
        title="Where the full model gains or loses fit quality",
    )
    ax.grid(axis="y", alpha=0.2)

    ax = axes[1, 0]
    ax.plot(inputs["wavelength"], leave_one_out, "o-", color="#6d28d9", ms=3.5)
    ax.axhline(base["delta_chi2_no_co2_minus_full"], color="#334155", ls="--", label="all bins")
    ax.set(
        xlabel="Deleted bin wavelength [micron]", ylabel="Delta chi-square after deletion",
        title="Leave-one-bin-out stability",
    )
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    labels = [row["setting"] for row in binning_rows + covariance_rows]
    values = [row["delta_chi2_no_co2_minus_full"] for row in binning_rows + covariance_rows]
    positions = np.arange(len(labels))
    ax.barh(positions, values, color=["#0f766e"] * len(binning_rows) + ["#2563eb"] * len(covariance_rows))
    ax.set(
        yticks=positions, yticklabels=labels,
        xlabel="Delta chi-square (no CO2 - full)",
        title="Binning and covariance stress tests",
    )
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.2)
    fig.suptitle(
        "WASP-39 b: robustness of the supplied full versus no-CO2 comparison",
        fontsize=15, weight="bold",
    )
    fig.savefig(FIGURE_FILE, dpi=190)
    plt.close(fig)

    return {
        "inputs": inputs,
        "base": base,
        "leave_one_out": leave_one_out,
        "binning_rows": binning_rows,
        "covariance_rows": covariance_rows,
        "scale": scale,
        "robust_amplitude_ppm": robust_amplitude,
        "amplitude_scale_heights": amplitude_scale_heights,
        "co2_window_delta_chi2": co2_window_delta,
    }


if __name__ == "__main__":
    result = main()
    print(
        "WASP-39 b: full-versus-no-CO2 diagnostic Delta chi-square="
        f"{result['base']['delta_chi2_no_co2_minus_full']:.1f}; "
        f"leave-one-bin range={result['leave_one_out'].min():.1f}-"
        f"{result['leave_one_out'].max():.1f}"
    )
