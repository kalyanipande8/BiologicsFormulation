import os
import sys

import pandas as pd
import pytest

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_DIR)

import simulate_formulation_data as sim  # noqa: E402
import analyze_formulations as af  # noqa: E402
import method_validation as mv  # noqa: E402


@pytest.fixture(scope="module")
def raw_df():
    return sim.simulate()


@pytest.fixture(scope="module")
def summary(raw_df):
    return af.summarize_variants(raw_df)


def test_at_least_150_variants_simulated(raw_df):
    assert raw_df["variant_id"].nunique() >= 150


def test_each_variant_has_replicates_and_required_readouts(raw_df):
    assert set(raw_df.groupby("variant_id").size().unique()) == {sim.N_REPLICATES}
    for col in ["pct_aggregation", "pct_acidic_variant", "pct_basic_variant",
                "pct_main_peak", "potency_retention_pct"]:
        assert col in raw_df.columns
        assert raw_df[col].notna().all()


def test_exactly_one_baseline_variant(raw_df):
    baseline_variants = raw_df.loc[raw_df["is_baseline"], "variant_id"].unique()
    assert len(baseline_variants) == 1


def test_anova_table_has_expected_factors(raw_df):
    _, anova_table = af.run_anova(raw_df)
    for factor in ["C(buffer)", "C(pH)", "C(stabilizer)", "C(conc_level)"]:
        assert factor in anova_table.index
    assert "PR(>F)" in anova_table.columns


def test_optimal_formulation_genuinely_outperforms_baseline(summary):
    """Sanity check: the statistically-identified optimal formulation must have
    lower mean aggregation AND higher (or equal) potency retention than the
    pre-defined baseline formulation - not just a different formulation."""
    baseline, optimal, agg_improvement_pct, potency_improvement_pp = af.find_optimal(summary)
    assert optimal["mean_pct_aggregation"] < baseline["mean_pct_aggregation"]
    assert agg_improvement_pct > 0
    assert optimal["mean_potency_retention_pct"] >= baseline["mean_potency_retention_pct"]


def test_optimal_is_the_global_minimum_aggregation_variant(summary):
    baseline, optimal, _, _ = af.find_optimal(summary)
    assert optimal["mean_pct_aggregation"] == pytest.approx(summary["mean_pct_aggregation"].min())


def test_method_validation_pass_fail_logic(raw_df):
    df = mv.main()
    # Recompute pass/fail independently from the stored result columns to confirm the logic,
    # rather than trusting the module's own boolean.
    for _, row in df.iterrows():
        crit = row["acceptance_criterion"]
        result = row["result"]
        if row["parameter"] == "Accuracy":
            lo, hi = [float(x.replace("%", "").strip()) for x in crit.replace("% recovery", "").split("-")]
            expected = lo <= result <= hi
        elif row["parameter"] == "Precision":
            max_cv = float(crit.split("<=")[1].replace("%", "").strip())
            expected = result <= max_cv
        else:  # Linearity
            min_r2 = float(crit.split(">=")[1].strip())
            expected = result >= min_r2
        assert bool(row["pass"]) == expected


def test_method_validation_overall_status_is_pass():
    df = mv.main()
    assert df["pass"].all(), f"Validation failures:\n{df[~df['pass']]}"
