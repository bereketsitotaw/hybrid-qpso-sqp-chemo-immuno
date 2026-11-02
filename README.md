# Hybrid QPSO–SQP Chemo-Immunotherapy Code

[![DOI](https://img.shields.io/badge/DOI-10.3934%2Fdcdsb.2026072-blue)](https://doi.org/10.3934/dcdsb.2026072)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2024a-orange)](requirements.txt)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

MATLAB code for the hybrid **Quantum PSO + SQP** method in:

> *Optimal Chemo–Immunotherapy Scheduling: A Hybrid QPSO–SQP Approach with Michaelis–Menten Pharmacodynamics*  
> B. S. Kidane, M. S. H. Motayed, S. Wang — [doi:10.3934/dcdsb.2026072](https://doi.org/10.3934/dcdsb.2026072)

**Authors:** Bereket Sitotaw Kidane, Md Samiul Haque Motayed, Shuo Wang (UT Arlington)

## Method overview

Michaelis–Menten saturation in both therapy channels yields a non-control-affine OCP. A QPSO warm start plus trapezoidal collocation / `fmincon` SQP recovers mixed boundary–interior dosing; PMP diagnostics corroborate the structure (DCDS-B Sections 2–5).

## Quick start

```matlab
cd hybrid-qpso-sqp-chemo-immuno
main          % Cases A, B, C
main('a')     % one case only
main('test')  % smoke tests (no full solve)
main('smoke') % tiny non-paper plumbing solve
```

Outputs: `results/generated/case_*_solution.mat` (and legacy `results/*.mat`), `results/table3.csv`, `results/figures/*.png`

## Requirements

- MATLAB R2020b+ (R2024a tested)
- Optimization Toolbox (`fmincon`, SQP)

## Cases (Section 5)

| Command | Description | `(A, B, C)` |
|---------|-------------|-------------|
| `main('a')` | Balanced objectives | `(1, 0.3, 1)` |
| `main('b')` | Asymmetric regularization | `(1, 0, 1)` |
| `main('c')` | Minimum-time limit | `(0, 0, 1)` |

Table 1 parameters: `configs/get_model_parameters.m`.  
Table 3 pin: `configs/reference_metrics.m`.

## Repository layout

```
main.m / init_project.m
configs/               Table 1/2 settings, cases, asserts, repro pins
src/chemo_immuno/      QPSO + SQP + PMP
src/utils/             Plots, Table 3, validation
tests/                 MATLAB smoke tests
docs/REPRODUCIBILITY.md
results/generated/     Run outputs (preferred)
results/reference/     Optional frozen mats + pin notes
.github/workflows/     Layout / publication-constant CI
```

## Reproducibility

See [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md). Seed `rng(42,'twister')`. Default path asserts publication Table 1 / case weights before each solve.

## Citation

```bibtex
@article{kidane2026chemoimmuno,
  title   = {Optimal Chemo--Immunotherapy Scheduling: A Hybrid {QPSO--SQP} Approach with Michaelis--Menten Pharmacodynamics},
  author  = {Kidane, Bereket Sitotaw and Motayed, Md Samiul Haque and Wang, Shuo},
  journal = {Discrete and Continuous Dynamical Systems -- Series B},
  year    = {2026},
  doi     = {10.3934/dcdsb.2026072}
}
```

Also see [CITATION.cff](CITATION.cff).

## License

MIT — see [LICENSE](LICENSE).

## Issues

Bug reports and reproduction questions: GitHub Issues on this repository.
