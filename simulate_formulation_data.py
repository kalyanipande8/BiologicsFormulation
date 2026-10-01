"""
simulate_formulation_data.py

Simulates an accelerated-stability formulation screen for a biologic drug
product (IgG1 mAb) across a full factorial of:

  buffer type      x  pH  x  stabilizer/excipient  x  concentration level

4 buffers x 5 pH levels x 4 stabilizers x 2 concentration levels = 160
formulation variants, each stressed at an accelerated condition (40 C /
75% RH, 4 weeks) with 3 replicate preps measured for:

  pct_aggregation     HMW species by SEC (%)
  pct_acidic_variant   acidic charge variants by iCIEF (%)
  pct_basic_variant    basic charge variants by iCIEF (%)
  pct_main_peak        main (native) charge species by iCIEF (%)
  potency_retention_pct  relative potency vs T0 (%), by cell-based bioassay

The degradation chemistry simulated is literature-plausible for an IgG1:
each buffer has a pH "sweet spot" (aggregation is a quadratic function of
distance from it), stabilizers reduce aggregation in a concentration-
dependent way (trehalose > arginine > sucrose > polysorbate-80 for
thermal-stress aggregation, consistent with published excipient
screening studies), charge-variant formation (deamidation/isomerization)
increases with pH, and potency retention falls with both aggregation and
charge heterogeneity. No real molecule or company data is used.

Writes data/raw/formulation_screen.csv. Uses a fixed random seed.
"""

import os

import numpy as np
import pandas as pd

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(PROJECT_DIR, "data", "raw")

RANDOM_SEED = 42
N_REPLICATES = 3
STRESS_CONDITION = "40C/75%RH, 4 weeks (accelerated)"

# buffer -> (pH optimum, baseline aggregation at optimum pH, curvature)
BUFFERS = {
    "Histidine":  dict(ph_opt=6.0, agg_base=0.90, curvature=0.40),
    "Citrate":    dict(ph_opt=5.5, agg_base=1.30, curvature=0.55),
    "Phosphate":  dict(ph_opt=7.0, agg_base=1.60, curvature=0.60),
    "Acetate":    dict(ph_opt=5.0, agg_base=1.10, curvature=0.50),
}

PH_LEVELS = [5.0, 5.5, 6.0, 6.5, 7.0]

# stabilizer -> (reduction in %agg at high conc, at low conc, unit label, [low, high] conc values)
STABILIZERS = {
    "Sucrose":        dict(reduction_high=0.50, reduction_low=0.25, unit="% w/v", levels=[1.0, 4.0]),
    "Trehalose":      dict(reduction_high=0.75, reduction_low=0.40, unit="% w/v", levels=[1.0, 4.0]),
    "Polysorbate-80": dict(reduction_high=0.30, reduction_low=0.15, unit="% w/v", levels=[0.02, 0.08]),
    "Arginine":       dict(reduction_high=0.65, reduction_low=0.35, unit="mM", levels=[50.0, 150.0]),
}

BASELINE_FORMULATION = dict(buffer="Histidine", pH=6.0, stabilizer="Sucrose", conc_level="low")


def _simulate_variant(rng, buffer, ph, stabilizer, conc_level):
    b = BUFFERS[buffer]
    s = STABILIZERS[stabilizer]
    conc_value = s["levels"][0] if conc_level == "low" else s["levels"][1]
    reduction = s["reduction_low"] if conc_level == "low" else s["reduction_high"]

    agg_mean = b["agg_base"] + b["curvature"] * (ph - b["ph_opt"]) ** 2 - reduction
    agg_mean = max(0.15, agg_mean)
    agg = np.clip(rng.normal(agg_mean, 0.10, N_REPLICATES), 0.05, None)

    # charge variants: acidic species grow with pH (deamidation) and mildly with aggregation
    acidic_mean = 7.5 + 1.4 * (ph - 6.0) + 0.3 * agg_mean
    acidic = np.clip(rng.normal(acidic_mean, 0.35, N_REPLICATES), 1.0, None)
    basic_mean = 3.0 + 0.25 * agg_mean + max(0.0, 6.0 - ph) * 0.2
    basic = np.clip(rng.normal(basic_mean, 0.25, N_REPLICATES), 0.5, None)
    main = 100 - acidic - basic

    # potency retention falls with aggregation and charge heterogeneity
    potency_mean = 100 - 2.2 * agg_mean - 0.35 * max(0.0, acidic_mean - 7.5)
    potency = np.clip(rng.normal(potency_mean, 0.8, N_REPLICATES), 50, 100)

    return agg, acidic, basic, main, potency, conc_value


def simulate():
    rng = np.random.default_rng(RANDOM_SEED)
    rows = []
    variant_id = 1
    for buffer in BUFFERS:
        for ph in PH_LEVELS:
            for stabilizer in STABILIZERS:
                for conc_level in ["low", "high"]:
                    agg, acidic, basic, main, potency, conc_value = _simulate_variant(
                        rng, buffer, ph, stabilizer, conc_level)
                    is_baseline = (buffer == BASELINE_FORMULATION["buffer"] and ph == BASELINE_FORMULATION["pH"]
                                   and stabilizer == BASELINE_FORMULATION["stabilizer"]
                                   and conc_level == BASELINE_FORMULATION["conc_level"])
                    for rep in range(1, N_REPLICATES + 1):
                        rows.append({
                            "variant_id": f"FORM-{variant_id:04d}",
                            "buffer": buffer,
                            "pH": ph,
                            "stabilizer": stabilizer,
                            "conc_level": conc_level,
                            "conc_value": conc_value,
                            "conc_unit": STABILIZERS[stabilizer]["unit"],
                            "stress_condition": STRESS_CONDITION,
                            "replicate": rep,
                            "pct_aggregation": round(agg[rep - 1], 3),
                            "pct_acidic_variant": round(acidic[rep - 1], 3),
                            "pct_basic_variant": round(basic[rep - 1], 3),
                            "pct_main_peak": round(main[rep - 1], 3),
                            "potency_retention_pct": round(potency[rep - 1], 3),
                            "is_baseline": is_baseline,
                        })
                    variant_id += 1

    df = pd.DataFrame(rows)
    os.makedirs(RAW_DIR, exist_ok=True)
    out_path = os.path.join(RAW_DIR, "formulation_screen.csv")
    df.to_csv(out_path, index=False)

    n_variants = df["variant_id"].nunique()
    print(f"Simulated {n_variants} formulation variants "
          f"({len(BUFFERS)} buffers x {len(PH_LEVELS)} pH levels x {len(STABILIZERS)} stabilizers x 2 conc levels), "
          f"{N_REPLICATES} replicates each -> {len(df)} rows")
    print(f"  {out_path}")
    return df


if __name__ == "__main__":
    simulate()
