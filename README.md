# BiologicsFormulation: Biologic Product Characterization & Formulation Optimization

A simulated drug-product development workflow for a biologic (IgG1 mAb): a
full-factorial accelerated-stability formulation screen, multi-factor
statistical analysis to identify the optimal formulation, analytical method
validation for the stability-indicating assays, and technical reports
supporting manufacturing scale-up and regulatory (CMC) submission.

> **Disclaimer:** This is a portfolio/demo project. All formulation,
> stability, and assay-validation data is **synthetically simulated** by
> `simulate_formulation_data.py` and `method_validation.py`. No real
> molecule, program, or employer data is used. The degradation chemistry
> (pH-dependent aggregation minima per buffer, stabilizer efficacy ranking,
> pH-driven charge-variant formation) is simulated from published,
> general biologics-formulation principles, not from any real dataset.

## What it does

1. **`simulate_formulation_data.py`** simulates a full-factorial
   accelerated-stability screen (40C/75% RH, 4 weeks): **4 buffers**
   (Histidine, Citrate, Phosphate, Acetate) x **5 pH levels** (5.0-7.0) x
   **4 stabilizers** (Sucrose, Trehalose, Polysorbate-80, Arginine) x
   **2 concentration levels** = **160 formulation variants**, each in
   triplicate, with %HMW aggregation (SEC), charge variants (iCIEF:
   %acidic/%main/%basic), and relative potency retention (bioassay).

2. **`analyze_formulations.py`** summarizes replicates per variant, runs a
   **Type II multi-factor ANOVA** of %aggregation against buffer, pH,
   stabilizer, and concentration level, and identifies the **optimal
   formulation** (lowest mean aggregation, potency retention as
   tie-breaker), computing a real stability-improvement metric vs. a
   pre-defined baseline (Histidine pH 6.0 + low-concentration sucrose, a
   common industry-default mAb formulation).

3. **`method_validation.py`** simulates and validates the two
   stability-indicating assays (SEC %aggregation, relative-potency
   bioassay) against ICH Q2(R1)-style accuracy/precision/linearity
   acceptance criteria, with real computed pass/fail outcomes.

4. **`generate_reports.py`** writes two technical reports populated with
   the actual computed numbers above, to `output/reports/`.

5. **`run_all.py`** runs the whole pipeline.

## How to run

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python run_all.py
pytest -q
```

## Results

_All numbers below are copied directly from an actual run of `python
run_all.py` - nothing is hardcoded. Fixed random seeds make them
reproducible._

- **160 formulation variants** simulated (4 buffers x 5 pH x 4 stabilizers
  x 2 concentration levels), 3 replicates each -> 480 raw measurement rows.
- **ANOVA on %aggregation** (Type II, replicate-level data): buffer
  (F = 110.7, p < 0.0001), pH (F = 16.7, p < 0.0001), stabilizer
  (F = 7.2, p = 0.0001), and concentration level (F = 22.0, p < 0.0001)
  are **all statistically significant** drivers of aggregation.
- **Baseline formulation** (Histidine, pH 6.0, Sucrose, low conc.):
  **0.649% aggregation**, 99.07% potency retention after accelerated
  stress.
- **Optimal formulation identified:** Histidine, pH 6.0, Trehalose, high
  conc.: **0.079% aggregation**, 99.19% potency retention.
- **Stability-improvement metric:** the optimal formulation reduces
  aggregation by **87.9%** relative to baseline, with a +0.11 percentage
  point gain in potency retention.

### Top 5 variants by lowest %aggregation

| Variant | Buffer | pH | Stabilizer | Conc. | %Aggregation | %Potency retention |
|---|---|---|---|---|---|---|
| FORM-0020 | Histidine | 6.0 | Trehalose | high | 0.079 | 99.19 |
| FORM-0012 | Histidine | 5.5 | Trehalose | high | 0.215 | 98.69 |
| FORM-0024 | Histidine | 6.0 | Arginine  | high | 0.227 | 99.38 |
| FORM-0028 | Histidine | 6.5 | Trehalose | high | 0.283 | 98.59 |
| FORM-0018 | Histidine | 6.0 | Sucrose   | high | 0.302 | 99.13 |

### Analytical method validation

| Assay | Parameter | Result | Acceptance criterion | Pass? |
|---|---|---|---|---|
| SEC %aggregation | Accuracy | 101.53% recovery | 90-110% | PASS |
| SEC %aggregation | Precision | 1.01% CV | <= 5% | PASS |
| SEC %aggregation | Linearity | R² = 0.99951 | >= 0.995 | PASS |
| Relative potency bioassay | Accuracy | 97.82% recovery | 80-120% | PASS |
| Relative potency bioassay | Precision | 5.73% CV | <= 15% | PASS |
| Relative potency bioassay | Linearity | R² = 0.99275 | >= 0.98 | PASS |

**Overall validation status: PASS (6/6 criteria met).** Both assays are
fit-for-purpose for supporting formulation stability assessment.

## Reports

- `output/reports/formulation_development_report.md` - process
  understanding, ANOVA results, optimal-formulation recommendation.
- `output/reports/method_validation_report.md` - assay validation summary
  and conclusion.
- `output/anova_table.csv`, `output/optimal_formulation.csv`,
  `output/method_validation.csv` - underlying computed results.
- `data/processed/variant_summary.csv` - per-variant mean readouts across
  all 160 formulations.

## Project structure

```
simulate_formulation_data.py   # 160-variant factorial stress screen -> data/raw/
analyze_formulations.py        # per-variant summary + ANOVA + optimal formulation -> data/processed/, output/
method_validation.py           # assay accuracy/precision/linearity validation -> output/
generate_reports.py            # technical reports from real computed results -> output/reports/
run_all.py                     # orchestrates the full pipeline
tests/test_pipeline.py         # pytest checks incl. injected optimal-vs-baseline sanity check
```

## Tech stack

- Python 3
- pandas / numpy for data simulation and wrangling
- statsmodels for the multi-factor (Type II) ANOVA
- scipy.stats for linear regression / linearity (R²)
- pytest
