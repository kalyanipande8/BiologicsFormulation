"""
analyze_formulations.py

Statistical analysis of the formulation screen:

1. Aggregates replicate-level data to per-variant means.
2. Runs a multi-factor ANOVA (Type II sum of squares) on %aggregation as a
   function of buffer, pH, stabilizer, and concentration level, to identify
   which factors significantly drive aggregation.
3. Identifies the optimal formulation (lowest mean %aggregation with
   potency retention as tie-breaker) and computes a stability-improvement
   metric vs. the pre-defined baseline formulation (Histidine pH 6.0 +
   low-concentration sucrose, a common industry-default mAb formulation).

Writes:
  data/processed/variant_summary.csv
  output/anova_table.csv
  output/optimal_formulation.csv
"""

import os

import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

import simulate_formulation_data as sim

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(PROJECT_DIR, "data", "raw")
PROCESSED_DIR = os.path.join(PROJECT_DIR, "data", "processed")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "output")


def summarize_variants(df):
    summary = df.groupby(
        ["variant_id", "buffer", "pH", "stabilizer", "conc_level", "conc_value", "conc_unit", "is_baseline"],
        as_index=False
    ).agg(
        mean_pct_aggregation=("pct_aggregation", "mean"),
        sd_pct_aggregation=("pct_aggregation", "std"),
        mean_pct_acidic_variant=("pct_acidic_variant", "mean"),
        mean_pct_main_peak=("pct_main_peak", "mean"),
        mean_potency_retention_pct=("potency_retention_pct", "mean"),
    )
    return summary


def run_anova(df):
    """Multi-factor ANOVA (Type II) of %aggregation ~ buffer + pH + stabilizer + conc_level,
    on the replicate-level data (not the per-variant means), with pH treated as categorical
    since the screen uses 5 discrete design levels."""
    model = smf.ols(
        "pct_aggregation ~ C(buffer) + C(pH) + C(stabilizer) + C(conc_level)",
        data=df
    ).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    return model, anova_table


def find_optimal(summary):
    baseline = summary.loc[summary["is_baseline"]].iloc[0]
    ranked = summary.sort_values(
        ["mean_pct_aggregation", "mean_potency_retention_pct"], ascending=[True, False]
    )
    optimal = ranked.iloc[0]

    agg_improvement_pct = (baseline["mean_pct_aggregation"] - optimal["mean_pct_aggregation"]) \
        / baseline["mean_pct_aggregation"] * 100
    potency_improvement_pp = optimal["mean_potency_retention_pct"] - baseline["mean_potency_retention_pct"]

    return baseline, optimal, agg_improvement_pct, potency_improvement_pp


def main():
    raw_path = os.path.join(RAW_DIR, "formulation_screen.csv")
    df = pd.read_csv(raw_path)

    summary = summarize_variants(df)
    os.makedirs(PROCESSED_DIR, exist_ok=True)
    summary.to_csv(os.path.join(PROCESSED_DIR, "variant_summary.csv"), index=False)

    model, anova_table = run_anova(df)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    anova_table.to_csv(os.path.join(OUTPUT_DIR, "anova_table.csv"))

    baseline, optimal, agg_improvement_pct, potency_improvement_pp = find_optimal(summary)
    optimal_df = pd.DataFrame([
        {"role": "baseline", **baseline.to_dict()},
        {"role": "optimal", **optimal.to_dict()},
        {"role": "improvement", "mean_pct_aggregation": agg_improvement_pct,
         "mean_potency_retention_pct": potency_improvement_pp},
    ])
    optimal_df.to_csv(os.path.join(OUTPUT_DIR, "optimal_formulation.csv"), index=False)

    print(f"{summary['variant_id'].nunique()} formulation variants analyzed")
    print("\nANOVA (Type II) on %aggregation:")
    print(anova_table.round(5))
    print(f"\nBaseline: {baseline['buffer']} pH {baseline['pH']}, {baseline['stabilizer']} "
          f"({baseline['conc_level']}) -> {baseline['mean_pct_aggregation']:.3f}% aggregation, "
          f"{baseline['mean_potency_retention_pct']:.2f}% potency retention")
    print(f"Optimal:  {optimal['buffer']} pH {optimal['pH']}, {optimal['stabilizer']} "
          f"({optimal['conc_level']}) -> {optimal['mean_pct_aggregation']:.3f}% aggregation, "
          f"{optimal['mean_potency_retention_pct']:.2f}% potency retention")
    print(f"Stability improvement: {agg_improvement_pct:.1f}% reduction in aggregation, "
          f"+{potency_improvement_pp:.2f} percentage points potency retention vs. baseline")

    return summary, model, anova_table, baseline, optimal, agg_improvement_pct, potency_improvement_pp


if __name__ == "__main__":
    main()
