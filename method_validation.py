"""
method_validation.py

Simulates analytical method validation for the two stability-indicating
assays used in the formulation screen (SEC %aggregation, relative-potency
bioassay), per standard ICH Q2(R1)-style validation parameters:

  Accuracy   % recovery of spiked/reference material at 3 levels (triplicate)
  Precision  intra-day %CV at a mid-range level (6 replicates)
  Linearity  R^2 of response vs. nominal value across 5 levels (triplicate)

Acceptance criteria (typical for a stability-indicating physicochemical /
bioassay method):
  SEC %aggregation assay:  accuracy 90-110%, precision %CV <= 5%, R^2 >= 0.995
  Relative potency assay:  accuracy 80-120%, precision %CV <= 15%, R^2 >= 0.98
  (potency bioassays are inherently higher-variability than a chromatographic
  assay, hence the wider acceptance window - this matches typical cell-based
  bioassay validation practice)

Writes output/method_validation.csv with real pass/fail outcomes computed
from the simulated data (nothing hardcoded).
"""

import os

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(PROJECT_DIR, "output")

RANDOM_SEED = 7

ASSAYS = {
    "SEC %aggregation": dict(
        levels=[0.5, 2.0, 5.0], true_recovery_bias=1.02, recovery_sd=0.03,
        precision_level=2.0, precision_sd_pct=1.8, n_precision=6,
        linearity_levels=[0.2, 1.0, 2.0, 5.0, 8.0], linearity_sd_pct=2.0,
        accuracy_lo=90.0, accuracy_hi=110.0, precision_cv_max=5.0, r2_min=0.995,
    ),
    "Relative potency bioassay": dict(
        levels=[70.0, 100.0, 130.0], true_recovery_bias=1.00, recovery_sd=0.06,
        precision_level=100.0, precision_sd_pct=7.0, n_precision=6,
        linearity_levels=[50.0, 75.0, 100.0, 125.0, 150.0], linearity_sd_pct=3.5,
        accuracy_lo=80.0, accuracy_hi=120.0, precision_cv_max=15.0, r2_min=0.98,
    ),
}


def validate_assay(rng, name, p):
    rows = []

    # --- Accuracy: % recovery at 3 levels, triplicate ---
    accuracy_recoveries = []
    for level in p["levels"]:
        measured = rng.normal(level * p["true_recovery_bias"], level * p["recovery_sd"], 3)
        recovery_pct = measured / level * 100
        accuracy_recoveries.extend(recovery_pct.tolist())
    mean_accuracy = float(np.mean(accuracy_recoveries))
    accuracy_pass = p["accuracy_lo"] <= mean_accuracy <= p["accuracy_hi"]
    rows.append({
        "assay": name, "parameter": "Accuracy", "metric": "mean %recovery (3 levels x triplicate)",
        "result": round(mean_accuracy, 2),
        "acceptance_criterion": f"{p['accuracy_lo']}-{p['accuracy_hi']}% recovery",
        "pass": accuracy_pass,
    })

    # --- Precision: intra-day %CV at mid-level, n replicates ---
    precision_data = rng.normal(p["precision_level"], p["precision_level"] * p["precision_sd_pct"] / 100,
                                 p["n_precision"])
    precision_cv = float(np.std(precision_data, ddof=1) / np.mean(precision_data) * 100)
    precision_pass = precision_cv <= p["precision_cv_max"]
    rows.append({
        "assay": name, "parameter": "Precision", "metric": f"intra-day %CV (n={p['n_precision']})",
        "result": round(precision_cv, 2),
        "acceptance_criterion": f"%CV <= {p['precision_cv_max']}%",
        "pass": precision_pass,
    })

    # --- Linearity: R^2 of measured vs. nominal across 5 levels, triplicate ---
    nominal, measured = [], []
    for level in p["linearity_levels"]:
        reps = rng.normal(level, level * p["linearity_sd_pct"] / 100, 3)
        nominal.extend([level] * 3)
        measured.extend(reps.tolist())
    slope, intercept, r_value, p_value, std_err = stats.linregress(nominal, measured)
    r2 = r_value ** 2
    linearity_pass = r2 >= p["r2_min"]
    rows.append({
        "assay": name, "parameter": "Linearity", "metric": f"R^2 ({len(p['linearity_levels'])} levels x triplicate)",
        "result": round(r2, 5),
        "acceptance_criterion": f"R^2 >= {p['r2_min']}",
        "pass": linearity_pass,
    })

    return rows


def main():
    rng = np.random.default_rng(RANDOM_SEED)
    all_rows = []
    for name, p in ASSAYS.items():
        all_rows.extend(validate_assay(rng, name, p))

    df = pd.DataFrame(all_rows)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    out_path = os.path.join(OUTPUT_DIR, "method_validation.csv")
    df.to_csv(out_path, index=False)

    print("Method validation results:")
    print(df.to_string(index=False))
    overall_pass = bool(df["pass"].all())
    print(f"\nOverall validation status: {'PASS' if overall_pass else 'FAIL'} "
          f"({df['pass'].sum()}/{len(df)} criteria met)")
    return df


if __name__ == "__main__":
    main()
