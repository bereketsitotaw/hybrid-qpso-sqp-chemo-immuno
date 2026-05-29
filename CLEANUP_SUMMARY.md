# Cleanup Summary

**Project:** QPSO_SQP reproducibility audit  
**Date:** 2026-05-28

---

## Canonical release selection

The publication-faithful implementation is the **modular chemo-immunotherapy solver** in `src/chemo_immuno/`, refactored from:

- `archive/legacy/qpso_optimal_original.m` (formerly `QPSO/qpso_optimal.m`)

---

## Archived files (preserved, not deleted)

### `archive/legacy/`

| Original path | Archived as | Reason |
|---------------|-------------|--------|
| `QPSO/qpso_optimal.m` | `qpso_optimal_original.m` | Source of truth — refactored into `src/` |
| `QPSO/qpso_bang_banng.m` | `qpso_bang_bang_original.m` | Control-affine comparison baseline (not paper Case C) |
| `QPSO/param_sensitivity.m` | `param_sensitivity_original.m` | Section 5.4 sensitivity (reduced grid) |
| `QPSO/qpso_chemo_immu.m` | `qpso_chemo_immu.m` | Outdated: linear kinetics, fixed T |
| `QPSO/Sensitivity_optimal.m` | `Sensitivity_optimal.m` | Experimental A/B ratio sweep |
| `solveOCP.asv` | `solveOCP.asv` | MATLAB editor auto-save |

### `archive/quadrotor_ocp/`

| File | Reason |
|------|--------|
| `solveOCP.m`, `ocp_transcription.m`, `qpso_knots_search.m` | Unified quadrotor OCP framework — different research thread |
| `define_quadrotor_nav_problem.m`, `run_quadcopter.m` | Quadrotor experiment driver |
| `plot_quadrotor_results.m`, `runPMPDiagnostics.m` | Quadrotor visualization / diagnostics |

### `archive/tumor_trait_model/`

| File | Reason |
|------|--------|
| `hybrid_qpso.m` | Trait-structured tumor QPSO; references missing `config.m`, `tumor_dynamics.m`, `compute_metrics.m` |

---

## Merged functionality

| Before (duplicated) | After (canonical) |
|---------------------|-------------------|
| Inline functions in `qpso_optimal.m` | `dynamics.m`, `objective_func.m`, `nonlinear_constraints.m`, `qpso_cost_func.m` |
| Monolithic QPSO loop | `run_qpso_phase.m` |
| Monolithic fmincon block | `run_sqp_refinement.m` |
| Inline PMP block | `compute_pmp_diagnostics.m` |
| Hard-coded Case A weights | `configs/case_a_balanced.m` (+ B, C variants) |
| Ad-hoc plotting | `plot_chemo_immuno_results.m` |

---

## Refactoring decisions

1. **Scientific behavior preserved:** Governing equations, collocation defects, QPSO update rule, and PMP formulas unchanged from `qpso_optimal_original.m`.
2. **Table 2 defaults applied in configs:** `Np=50`, `kmax=100`, `ρ=1e6`, SQP tolerances `1e-6` — original script sometimes used larger swarm/iter counts; documented in `REPRODUCIBILITY_REPORT.md`.
3. **Case C uses MM kinetics:** Not merged with `qpso_bang_bang_original.m` (different dynamical model).
4. **No permanent deletion:** All `.m` sources moved to `archive/` except new modular files.
5. **Figures relocated:** EPS/PNG outputs → `results/figures/`.

---

## Files removed from project root

Root directory now contains only release artifacts (`README.md`, reports, configs folder, etc.). Former root `.m` files relocated to `archive/` or `src/`.

Duplicate PDF at root (`10.3934_dcdsb.2026072.pdf`) retained alongside `paper/published_paper.pdf` for convenience; canonical copy is under `paper/`.

---

## Not merged (intentionally separate)

- Quadrotor OCP solver (`solveOCP` framework) — future standalone repo candidate
- Trait-based tumor model — prior CDC/QPSO work, different paper
- Control-affine bang-bang script — comparison to [Ledzewicz et al., 2024], not this paper's Case C
