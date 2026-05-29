# Hybrid QPSO–SQP Chemo-Immunotherapy Code

[![DOI](https://img.shields.io/badge/DOI-10.3934%2Fdcdsb.2026072-blue)](https://doi.org/10.3934/dcdsb.2026072)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2024a-orange)](requirements.txt)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

MATLAB code for the hybrid **Quantum PSO + SQP** method in:

> *Optimal Chemo–Immunotherapy Scheduling: A Hybrid QPSO–SQP Approach with Michaelis–Menten Pharmacodynamics*  
> B. S. Kidane, M. S. H. Motayed, S. Wang — [doi:10.3934/dcdsb.2026072](https://doi.org/10.3934/dcdsb.2026072)

**Authors:** Bereket Sitotaw Kidane, Md Samiul Haque Motayed, Shuo Wang (UT Arlington)

## Quick start

```matlab
cd hybrid-qpso-sqp-chemo-immuno
main          % Cases A, B, C
main('a')     % one case only
```

Outputs: `results/case_*_solution.mat`, `results/table3.csv`, `results/figures/*.png`

## Requirements

- MATLAB R2020b+ (R2024a tested)
- Optimization Toolbox (`fmincon`, SQP)

## Cases (Section 5)

| Command | Description | `(A, B, C)` |
|---------|-------------|-------------|
| `main('a')` | Balanced objectives | `(1, 0.3, 1)` |
| `main('b')` | Asymmetric regularization | `(1, 0, 1)` |
| `main('c')` | Minimum-time limit | `(0, 0, 1)` |

## Repository layout

```
main.m                 Run here
init_project.m         Path setup
configs/               Parameters and case studies
src/chemo_immuno/      QPSO + SQP solver
src/utils/             Plots and Table 3 check
results/               Generated locally (not in git)
```

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

## License

MIT — see [LICENSE](LICENSE).
