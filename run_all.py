"""
run_all.py

Runs the full pipeline end-to-end:
  1. Simulate the formulation screen (160 variants) -> data/raw/
  2. Summarize + multi-factor ANOVA + identify optimal formulation -> data/processed/, output/
  3. Analytical method validation (accuracy/precision/linearity) -> output/
  4. Generate technical reports -> output/reports/

Run with:  python run_all.py
"""

import simulate_formulation_data
import analyze_formulations
import method_validation
import generate_reports


def main():
    print("--- Simulating formulation screen ---")
    simulate_formulation_data.simulate()

    print("\n--- Statistical analysis (ANOVA, optimal formulation) ---")
    summary, model, anova_table, baseline, optimal, agg_improvement_pct, potency_improvement_pp = \
        analyze_formulations.main()

    print("\n--- Analytical method validation ---")
    validation_df = method_validation.main()

    print("\n--- Generating technical reports ---")
    generate_reports.main(summary, anova_table, baseline, optimal, agg_improvement_pct,
                           potency_improvement_pp, validation_df)


if __name__ == "__main__":
    main()
