# Figure-to-Code Mapping

Paper: *Optimal Chemo–Immunotherapy Scheduling* (doi:[10.3934/dcdsb.2026072](https://doi.org/10.3934/dcdsb.2026072))

All deterministic figure scripts use pinned settings from `configs/reproducibility_settings.m` (`rng(42,'twister')`) and write outputs under `results/figures/`.

## One-command reproduction

| Action | Script | Primary outputs |
|--------|--------|-----------------|
| Full Section 5 + Table 3 | `scripts/reproduce_paper.m` | `results/case_*_solution.mat`, `results/tables/table3.csv`, `results/figures/*.png` |
| Validate vs Table 3 | `scripts/validate_against_reference.m` | Console report |
| Freeze author reference | `reproduce_paper_pipeline('freeze_reference', true)` | `results/reference/canonical_outputs.mat` |

## Figure map

| Paper Figure | Section | Script | Output prefix | Status |
|--------------|---------|--------|---------------|--------|
| **Figure 1** | 5.1 Case A states/controls | `scripts/reproduce_figure_1.m` | `results/figures/fig1_caseA_*` | ✅ Deterministic |
| **Figure 2** | 5.1 Case A PMP diagnostics | `scripts/reproduce_figure_2.m` | `results/figures/fig2_pmp_caseA_*` | ✅ Deterministic |
| **Figure 3** | 5.2 Case B asymmetric | `scripts/reproduce_figure_3.m` | `results/figures/fig3_caseB_*`, `fig3_pmp_caseB_*` | ✅ Deterministic |
| **Figure 4** | 5.3 Case C min-time | `scripts/reproduce_figure_4.m` | `results/figures/fig4_caseC_*`, `fig4_pmp_caseC_*` | ✅ Deterministic |
| **Figure 5** | 5.4 Bound sweep | `scripts/reproduce_figure_5.m` | `results/figures/` (archived EPS) | ⚠️ Partial — sweep driver not yet modularized |
| **Figure 6** | 5.4 Parameter tornado | `scripts/reproduce_figure_6.m` | `results/figures/tornado.eps` (archived) | ⚠️ Partial — uses archived weight sensitivity |
| **Table 3** | 5.5 Comparative metrics | `scripts/generate_table3.m` | `results/tables/table3.csv` | ✅ Deterministic |

## Archived figure assets (author-generated)

Pre-refactor EPS/PNG files from the original submission are preserved in `results/figures/`:

| File | Likely paper mapping |
|------|---------------------|
| `Case_A.eps`, `Case_A_PMP.eps` | Figures 1–2 |
| `Case_B.eps`, `Case_B_PMP.eps` | Figure 3 |
| `Case_C.eps`, `Case_C_PMP.eps` | Figure 4 |
| `Chemo_sensitivity.eps`, `Immuno_Sensitivity.eps` | Sensitivity panels |
| `tornado.eps`, `Pareto_chemo.eps` | Figure 6 / Pareto |

## Case study → config mapping

| Case | Config | Cost weights `(ℓ₁,ℓ₂,ℓ₃,A,B,C)` |
|------|--------|----------------------------------|
| A | `configs/case_a_balanced.m` | `(0,0,0; 1, 0.3, 1)` |
| B | `configs/case_b_asymmetric.m` | `(0,0,0; 1, 0, 1)` |
| C | `configs/case_c_min_time.m` | `(0,0,0; 0, 0, 1)` |

## Reference metrics

Published Table 3 values are frozen in:

- `configs/reference_metrics.m`
- `results/reference/table3_paper.csv`

Validation tolerances are defined in `configs/reproducibility_settings.m`.

## Not paper figures (archived code)

| Item | Location | Note |
|------|----------|------|
| Control-affine bang-bang baseline | `archive/legacy/qpso_bang_bang_original.m` | **Not** Michaelis–Menten Case C |
| Quadrotor navigation | `archive/quadrotor_ocp/` | Unrelated project |
| Trait-structured tumor model | `archive/tumor_trait_model/` | Unrelated project |
