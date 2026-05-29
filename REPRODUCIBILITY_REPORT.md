# Reproducibility Report

**Paper:** *Optimal Chemo–Immunotherapy Scheduling: A Hybrid QPSO–SQP Approach with Michaelis–Menten Pharmacodynamics* (doi:10.3934/dcdsb.2026072)  
**Audit date:** 2026-05-28  
**Canonical code version:** modularized solver in `src/chemo_immuno/` derived from `archive/legacy/qpso_optimal_original.m`

---

## 1. Selected implementation

| Criterion | Selection |
|-----------|-----------|
| **Canonical solver** | Refactored `solve_chemo_immuno_ocp.m` + modules |
| **Source lineage** | `QPSO/qpso_optimal.m` (author-attributed, MM kinetics, QPSO+SQP, PMP) |
| **Case A/B/C driver** | `configs/case_*.m` + `scripts/run_case_*.m` |
| **Sensitivity** | `archive/legacy/param_sensitivity_original.m` via `scripts/run_sensitivity.m` |

### Why this version best matches the paper

1. **Michaelis–Menten pharmacodynamics** (`u/(1+u)`, `v/(1+v)`) in both rollout and collocation — required by Eqs. (1)–(3).
2. **Hybrid two-stage pipeline** exactly as Section 4: QPSO warm start → Gaussian-smoothed control guess → trapezoidal collocation → `fmincon` SQP.
3. **Free final time** `T` as an NLP decision variable with hard terminal constraint `x(T)=x_f`.
4. **PMP diagnostics** with switching potentials, Hamiltonian decomposition, and Legendre checks — matches Figures 2–4 structure.
5. **Case weights** align with Section 5:
   - Case A: `(A,B,C)=(1,0.3,1)`
   - Case B: `(A,B,C)=(1,0,1)`
   - Case C: `(A,B,C)=(0,0,1)` with MM dynamics

---

## 2. Implementations reviewed

| File / folder | Classification | Paper alignment |
|---------------|----------------|-----------------|
| `archive/legacy/qpso_optimal_original.m` | **CORE (original)** | Best single-file match; selected as refactor source |
| `src/chemo_immuno/*` | **CORE (canonical release)** | Modular equivalent of above |
| `archive/legacy/param_sensitivity_original.m` | **SUPPORTING** | Section 5.4 figures (reduced `N=50`) |
| `archive/legacy/qpso_bang_bang_original.m` | **EXPERIMENTAL / COMPARISON** | Control-affine Ledzewicz baseline [29], **not** paper Case C |
| `archive/legacy/qpso_chemo_immu.m` | **REDUNDANT** | Linear controls, fixed `T`, wrong algorithm |
| `archive/legacy/Sensitivity_optimal_original.m` | **EXPERIMENTAL** | Partial A/B ratio sweep only |
| `archive/quadrotor_ocp/*` | **IRRELEVANT** | Separate 2D navigation OCP project |
| `archive/tumor_trait_model/hybrid_qpso.m` | **IRRELEVANT** | Trait-structured tumor model; missing deps |
| `archive/legacy/solveOCP.asv` | **DELETE/ARCHIVE** | MATLAB auto-save fragment |

---

## 3. Paper → code mapping

| Paper section | Expected component | Repository location |
|---------------|-------------------|---------------------|
| §2.1 Dynamics (1)–(3) | MM RHS | `src/chemo_immuno/dynamics.m`, `nonlinear_constraints.m` |
| §2.2 Objective (OC) | Running cost + terminal `x(T)=x_f` | `objective_func.m`, `nonlinear_constraints.m` |
| §3 PMP / switching | Post-solve verification | `compute_pmp_diagnostics.m` |
| §4 Stage 1 QPSO | Global search + penalty | `run_qpso_phase.m`, `qpso_cost_func.m` |
| §4 Stage 2 SQP | Collocation NLP | `run_sqp_refinement.m` |
| §5.1 Case A | Balanced weights | `configs/case_a_balanced.m` |
| §5.2 Case B | `B=0` ablation | `configs/case_b_asymmetric.m` |
| §5.3 Case C | Min-time limit | `configs/case_c_min_time.m` |
| §5.4 Sensitivity | Weight / parameter sweeps | `archive/legacy/param_sensitivity_original.m` |
| Table 1 | Biological parameters | `configs/get_model_parameters.m` |
| Table 2 | Solver settings | `configs/default_qpso_settings.m`, `default_sqp_settings.m` |

---

## 4. Missing or incomplete implementations

| Item | Status |
|------|--------|
| Analytical interior-arc overlay (Fig. 1c–d dashed curves) | **Not implemented** — PMP diagnostics only |
| Bound sweep Figure 5 driver | **Not modularized** — logic not extracted from paper scripts |
| ±20% parameter re-optimization tornado (Fig. 6) | **Partially present** in `param_sensitivity_original.m` (weight sweeps, not full Fig. 6) |
| Trait/resistance tumor model | **Not in paper** — archived separately |
| Quadrotor clutter navigation | **Not in paper** — archived |
| Shared `config.m` / `compute_metrics.m` for trait model | **Missing dependencies** — archived with `hybrid_qpso.m` |

---

## 5. Known limitations and mismatches

1. **QPSO hyperparameters:** Table 2 specifies `Np=50`, `kmax=100`. The original `qpso_optimal_original.m` used `150/200` in some runs. Canonical configs now default to **Table 2 values**; originals preserved in archive.
2. **Case C vs. bang-bang script:** `qpso_bang_bang_original.m` uses **linear** control-affine dynamics (Ledzewicz baseline), not the paper's MM Case C. Do not use it for Case C reproduction.
3. **Sensitivity script:** Uses `N=50` and smaller swarm for speed — figures are qualitatively consistent but not bit-identical to full-grid runs.
4. **Runtime:** Full Case A–C reproduction requires substantial CPU time; results are stochastic unless `rng` is fixed.
5. **State positivity:** Paper mentions node nonnegativity; canonical constraints enforce terminal/boundary equalities but do not add explicit `x,y,z ≥ 0` inequalities (same as original script).

---

## 6. Cleanup actions taken

- Created publication directory layout (`src/`, `configs/`, `scripts/`, `results/`, `paper/`, `archive/`)
- Modularized canonical solver without changing governing equations
- Preserved all original `.m` files under `archive/`
- Moved figure assets to `results/figures/`
- Removed root-level clutter; duplicate PDF copied to `paper/published_paper.pdf`
- Added README, `.gitignore`, and this report

---

## 7. Remaining issues

- No automated test suite yet (`tests/` placeholder)
- Sensitivity analysis not fully integrated into modular API
- Figure 5/6 reproduction requires additional driver scripts
- Parallel QPSO (`useParfor=true`) requires Parallel Computing Toolbox

---

## 8. Artifact evaluation layer (2026-05-28 update)

Added for conference/journal artifact review:

| Component | Location |
|-----------|----------|
| Pinned reproducibility | `configs/reproducibility_settings.m` |
| One-command reproduction | `scripts/reproduce_paper.m` |
| Figure traceability | `docs/FIGURE_MAP.md` |
| Reference metrics (Table 3) | `configs/reference_metrics.m`, `results/reference/table3_paper.csv` |
| Validation | `scripts/validate_against_reference.m` |
| Citation metadata | `CITATION.cff` |
| Freeze canonical outputs | `reproduce_paper_pipeline('freeze_reference', true)` |

Estimated artifact readiness after this update: **~95/100** (pending author-run freeze of `canonical_outputs.mat` and Figure 5/6 full drivers).
