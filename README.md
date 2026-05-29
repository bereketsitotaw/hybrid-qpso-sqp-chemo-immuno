# Hybrid QPSO–SQP for Optimal Chemo–Immunotherapy Scheduling

[![DOI](https://img.shields.io/badge/DOI-10.3934%2Fdcdsb.2026072-blue)](https://doi.org/10.3934/dcdsb.2026072)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2024a-orange)](requirements.txt)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Reproducibility](https://img.shields.io/badge/Reproducibility-pinned%20settings-informational)](configs/reproducibility_settings.m)
[![Version](https://img.shields.io/badge/release-v1.0.0--paper-lightgrey)](CITATION.cff)

Reproducible MATLAB code for the paper:

> **Optimal Chemo–Immunotherapy Scheduling: A Hybrid QPSO–SQP Approach with Michaelis–Menten Pharmacodynamics**  
> Bereket Sitotaw Kidane, Md Samiul Haque Motayed, and Shuo Wang  
> *Discrete and Continuous Dynamical Systems – Series B* ([doi:10.3934/dcdsb.2026072](https://doi.org/10.3934/dcdsb.2026072))

## Quick start (one command)

```matlab
cd /path/to/QPSO_SQP
run('scripts/reproduce_paper.m')
```

This runs Cases A/B/C, exports figures, writes `results/tables/table3.csv`, and validates against published Table 3 metrics.

To freeze author reference trajectories after a successful run:

```matlab
reproduce_paper_pipeline('freeze_reference', true);
```

## Reproducibility notes

The published results were generated using:

| Setting | Value |
|---------|-------|
| MATLAB | R2024a |
| Optimization Toolbox | v9.5 (update in `configs/reproducibility_settings.m` if needed) |
| QPSO random seed | `rng(42,'twister')` |
| QPSO rollout (`ode45`) relative tolerance | `1e-8` |
| QPSO rollout (`ode45`) absolute tolerance | `1e-10` |
| QPSO swarm size / iterations | `50` / `100` (Table 2) |
| Terminal warm-start penalty | `ρ = 1e6` |
| SQP tolerances | `1e-6` (constraints and optimality) |
| Collocation grid | `N = 100` |

All pinned values live in [`configs/reproducibility_settings.m`](configs/reproducibility_settings.m) and are applied via [`apply_reproducibility_settings`](src/utils/apply_reproducibility_settings.m).

**Note:** Paper Table 2 documents `ode45` tolerances of `1e-3` for the QPSO warm-start phase. The pinned author settings above use tighter tolerances for artifact reproduction; both are documented in the config file.

To reproduce paper figures deterministically:

```matlab
rng(42,'twister');
run('scripts/reproduce_paper.m')
```

Due to stochastic QPSO initialization, small deviations in terminal time and control effort may occur if random seeds, ODE tolerances, or MATLAB versions differ. Validation tolerances are defined in `reproducibility_settings.m` and checked by [`validate_against_reference.m`](scripts/validate_against_reference.m).

## Overview

This repository implements a **two-stage hybrid optimizer** for a three-state tumor–immune optimal control problem with **Michaelis–Menten saturating pharmacodynamics**:

1. **Stage 1 (QPSO):** global warm start via Quantum Particle Swarm Optimization with `ode45` rollout and a terminal tumor penalty.
2. **Stage 2 (SQP):** direct trapezoidal collocation refined with MATLAB `fmincon` (SQP algorithm).

## Figure-to-code traceability

See [`docs/FIGURE_MAP.md`](docs/FIGURE_MAP.md) for the full mapping of paper figures/tables to scripts and output paths.

| Paper artifact | Script | Output |
|----------------|--------|--------|
| Figure 1 | `scripts/reproduce_figure_1.m` | `results/figures/fig1_caseA_*` |
| Figure 2 | `scripts/reproduce_figure_2.m` | `results/figures/fig2_pmp_caseA_*` |
| Figure 3 | `scripts/reproduce_figure_3.m` | `results/figures/fig3_caseB_*` |
| Figure 4 | `scripts/reproduce_figure_4.m` | `results/figures/fig4_caseC_*` |
| Table 3 | `scripts/generate_table3.m` | `results/tables/table3.csv` |

## ⚠️ Archived code warning

`archive/legacy/qpso_bang_bang_original.m` is preserved for **historical reference only** and does **not** correspond to the published Michaelis–Menten formulation used in the paper. It implements a **control-affine baseline** inspired by earlier Ledzewicz-type formulations and should **not** be used to reproduce published results or Case C in this paper.

Similarly, `archive/quadrotor_ocp/` and `archive/tumor_trait_model/` are unrelated research threads.

## Repository structure

```
├── paper/published_paper.pdf
├── src/chemo_immuno/          Canonical solver
├── src/utils/                 Plotting, validation, reproducibility helpers
├── configs/                   Case studies + pinned settings
├── scripts/                   reproduce_paper.m, figure scripts, validators
├── results/
│   ├── figures/               Generated + archived EPS/PNG
│   ├── tables/                table3.csv
│   └── reference/             Frozen metrics + optional canonical_outputs.mat
├── archive/                   Legacy / unrelated code (preserved)
├── docs/FIGURE_MAP.md
├── CITATION.cff
├── REPRODUCIBILITY_REPORT.md
└── CLEANUP_SUMMARY.md
```

## Requirements

- **MATLAB** R2020b or later (**R2024a** for exact author environment)
- **Optimization Toolbox** (`fmincon`, SQP)
- **Parallel Computing Toolbox** (optional)

See [`requirements.txt`](requirements.txt).

## Usage

### Individual cases

```matlab
run('scripts/run_case_a.m')   % Section 5.1
run('scripts/run_case_b.m')   % Section 5.2
run('scripts/run_case_c.m')   % Section 5.3
```

### Validate against paper Table 3

```matlab
bootstrap_project();
% after running cases:
validate_against_reference(summaries, 'verbose', true);
```

Or run the full pipeline (includes validation):

```matlab
run('scripts/reproduce_paper.m')
```

### Programmatic API

```matlab
bootstrap_project();
cfg = apply_reproducibility_settings(case_a_balanced());
[solution, summary] = solve_chemo_immuno_ocp(cfg);
plot_chemo_immuno_results(solution, summary, 'results/figures');
```

## Evaluation metrics

Key outputs: final time `T*`, objective `J*`, terminal tumor `x(T)`, integrated exposures `∫u dt` and `∫v dt`, PMP stationarity residuals, and Table 3 comparison vs `results/reference/table3_paper.csv`.

## Citation

GitHub renders citation metadata from [`CITATION.cff`](CITATION.cff). BibTeX:

```bibtex
@article{kidane2026chemoimmuno,
  title   = {Optimal Chemo--Immunotherapy Scheduling: A Hybrid {QPSO--SQP} Approach with Michaelis--Menten Pharmacodynamics},
  author  = {Kidane, Bereket Sitotaw and Motayed, Md Samiul Haque and Wang, Shuo},
  journal = {Discrete and Continuous Dynamical Systems -- Series B},
  year    = {2026},
  doi     = {10.3934/dcdsb.2026072}
}
```

## Contact

- Bereket Sitotaw Kidane — UT Arlington
- Corresponding author: Shuo Wang

See [`REPRODUCIBILITY_REPORT.md`](REPRODUCIBILITY_REPORT.md) for the full audit trail.
