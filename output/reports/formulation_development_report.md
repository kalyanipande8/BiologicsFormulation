# Formulation Development Report: Biologic Drug Product Stability Screen

## 1. Purpose

This report documents the characterization of a biologic drug product (IgG1 mAb) across **160 formulation variants** varying buffer type, pH, stabilizer/excipient identity, and concentration, under an accelerated stress condition (40C/75% RH, 4 weeks). The objective was to identify the buffer/pH/stabilizer combination that minimizes aggregation and charge-variant formation while maximizing potency retention, to support formulation lock, manufacturing scale-up, and regulatory (CMC) submission.

## 2. Design of Experiments

A full-factorial screen was executed: 4 buffers (Histidine, Citrate, Phosphate, Acetate) x 5 pH levels (5.0-7.0) x 4 stabilizers (Sucrose, Trehalose, Polysorbate-80, Arginine) x 2 concentration levels (low/high) = **160 unique variants**, each measured in triplicate for %HMW aggregation (SEC), charge variants (iCIEF), and relative potency.

## 3. Statistical Analysis (Multi-Factor ANOVA)

A Type II ANOVA was fit to %aggregation as a function of buffer, pH, stabilizer, and concentration level, using all replicate-level observations (160 variants x 3 replicates each):

```
                  sum_sq     df          F  PR(>F)
C(buffer)      118.52763    3.0  110.73061  0.0000
C(pH)           23.77211    4.0   16.65625  0.0000
C(stabilizer)    7.67607    3.0    7.17112  0.0001
C(conc_level)    7.83823    1.0   21.96784  0.0000
Residual       166.98464  468.0        NaN     NaN
```

Factors significant at p < 0.05: **C(buffer), C(pH), C(stabilizer), C(conc_level)**. This confirms that buffer identity, solution pH, and stabilizer choice are statistically meaningful drivers of aggregation under accelerated stress, supporting a design-space-based formulation strategy rather than single-factor selection.

## 4. Optimal Formulation

**Baseline formulation** (common industry-default): Histidine, pH 6.0, Sucrose (low conc.) -> 0.649% aggregation, 99.07% potency retention after stress.

**Optimal formulation identified:** Histidine, pH 6.0, Trehalose (high conc.) -> 0.079% aggregation, 99.19% potency retention after stress.

**Stability-improvement metric:** optimal vs. baseline reduces aggregation by **87.9%** and improves potency retention by **0.11 percentage points**.

### Top 5 variants by lowest %aggregation

```
variant_id    buffer  pH stabilizer conc_level  mean_pct_aggregation  mean_potency_retention_pct
 FORM-0020 Histidine 6.0  Trehalose       high              0.078667                   99.186333
 FORM-0012 Histidine 5.5  Trehalose       high              0.215000                   98.688667
 FORM-0024 Histidine 6.0   Arginine       high              0.227333                   99.382667
 FORM-0028 Histidine 6.5  Trehalose       high              0.282667                   98.593000
 FORM-0018 Histidine 6.0    Sucrose       high              0.302000                   99.129000
```

## 5. Process Understanding and Recommendations

- Aggregation follows a quadratic relationship with pH within each buffer system, consistent with colloidal/conformational stability being maximized near each buffer's characteristic pH optimum.
- Trehalose at the high concentration level consistently outperformed sucrose and polysorbate-80 as a stabilizer against thermally-induced aggregation, consistent with its higher glass-transition temperature and preferential-exclusion behavior.
- Charge-variant (acidic species) formation increases with pH, consistent with base-catalyzed deamidation kinetics; this creates a trade-off that the optimal formulation balances against aggregation risk.
- **Recommendation:** advance the optimal formulation identified above to confirmatory forced-degradation and real-time/accelerated stability studies, and use it as the basis for the drug product formulation design space in the CMC section of the regulatory submission.
