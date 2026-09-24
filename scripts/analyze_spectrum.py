from __future__ import annotations
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "spectra"
FIGURES = ROOT / "figures"
STATS_FILE = FIGURES / "spectrum_statistics.csv"

def flat_test(values, errors):
    values, errors = np.asarray(values, float), np.asarray(errors, float)
    good = np.isfinite(values) & np.isfinite(errors) & (errors > 0)
    values, errors = values[good], errors[good]
    weights = 1 / errors**2
    mean = np.sum(weights * values) / np.sum(weights)
    statistic = np.sum(((values - mean) / errors)**2)
    dof = len(values) - 1
    return {"n": len(values), "mean": mean, "chi2": statistic, "dof": dof,
            "p": chi2.sf(statistic, dof)}

def offset_model_test(wavelength, values, errors, model_wavelength, model_values):
    if wavelength.min() < model_wavelength.min() or wavelength.max() > model_wavelength.max():
        raise ValueError("observed wavelengths extend beyond model support")
    model = np.interp(wavelength, model_wavelength, model_values)
    weights = 1 / errors**2
    offset = np.sum(weights * (values - model)) / np.sum(weights)
    statistic = np.sum(((values - model - offset) / errors)**2)
    dof = len(values) - 1
    return {"n": len(values), "offset": offset, "chi2": statistic,
            "dof": dof, "p": chi2.sf(statistic, dof), "model": model + offset}

def write_rows(rows):
    fields = sorted({key for row in rows for key in row})
    with STATS_FILE.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)

from astropy.table import Table

FIGURE_FILE = FIGURES / "wasp39b_published_spectrum.png"

def main():
    FIGURES.mkdir(exist_ok=True)
    table = Table.read(DATA / "eureka_transmission_spectrum.txt", format="ascii.ecsv")
    wavelength = np.asarray(table["wavelength"], float)
    depth = np.asarray(table["tr_depth"], float) * 1e6
    error = .5 * (np.asarray(table["tr_depth_errneg"], float) + np.asarray(table["tr_depth_errpos"], float)) * 1e6
    flat = flat_test(depth, error)
    rows = [{"comparison": "weighted flat", **flat}]
    raw_models = [
        ("ScCHIMERA full", np.loadtxt(DATA / "scchimeramodel.txt")),
        ("ScCHIMERA remove CO2", np.loadtxt(DATA / "scchimeramodel_no_co2.txt")),
    ]
    support_min = max(model[:, 0].min() for _, model in raw_models)
    support_max = min(model[:, 0].max() for _, model in raw_models)
    width = np.asarray(table["bin_width"], float)
    supported = ((wavelength - width / 2) >= support_min) & ((wavelength + width / 2) <= support_max)
    model_wavelength = wavelength[supported]
    model_depth = depth[supported]
    model_error = error[supported]
    models = []
    for label, model in raw_models:
        fit = offset_model_test(model_wavelength, model_depth, model_error,
                                model[:, 0], model[:, 1] * 1e6)
        rows.append({"comparison": label + " + fitted vertical offset",
                     **{key: value for key, value in fit.items() if key != "model"}})
        models.append((label, fit["model"]))
    write_rows(rows)
    fig, ax = plt.subplots(figsize=(9.2, 5.3))
    ax.errorbar(wavelength, depth, yerr=error, fmt="o", ms=3.2, color="#17212b",
                ecolor="#78909c", alpha=.85, label="Eureka reduction")
    for label, values in models:
        ax.plot(model_wavelength, values, lw=2, label=label + " (offset fitted)")
    ax.set(xlabel="Wavelength [micron]", ylabel="Transit depth [ppm]",
           title="WASP-39 b: published JWST NIRSpec PRISM transmission spectrum")
    ax.grid(alpha=.2); ax.legend(frameon=False, fontsize=8); fig.tight_layout()
    fig.savefig(FIGURE_FILE, dpi=190); plt.close(fig)
    return {"flat": flat, "rows": rows, "n": len(depth),
            "common_support_n": int(supported.sum()),
            "excluded_model_bins": int((~supported).sum())}

if __name__ == "__main__":
    result = main(); print(f"WASP-39 b: {result['n']} spectral bins, "
                           f"{result['common_support_n']} in model comparisons; "
                           f"flat-spectrum p={result['flat']['p']:.3g}")
