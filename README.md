# Hybrid QPSO–SQP Chemo-Immunotherapy Code

MATLAB implementation of the hybrid **Quantum PSO + SQP** method for optimal chemo-immunotherapy scheduling with **Michaelis–Menten** pharmacodynamics.

**Paper:** [doi:10.3934/dcdsb.2026072](https://doi.org/10.3934/dcdsb.2026072)

## Requirements

- MATLAB R2020b+ (R2024a tested)
- Optimization Toolbox (`fmincon`, SQP)

## Quick start

```matlab
cd hybrid-qpso-sqp-chemo-immuno   % repo folder
main                              % Cases A, B, C
main('a')                         % Case A only
```

Outputs go to `results/` (`.mat` solutions, `table3.csv`, PNG figures).

## Cases

| Call | Section | Weights `(A,B,C)` |
|------|---------|-------------------|
| `main('a')` | Balanced objectives | `(1, 0.3, 1)` |
| `main('b')` | Asymmetric regularization | `(1, 0, 1)` |
| `main('c')` | Minimum-time limit | `(0, 0, 1)` |

## Layout

```
main.m              Entry point
init_project.m      Path setup
configs/            Model + case parameters
src/chemo_immuno/   QPSO, SQP, dynamics, PMP
src/utils/          Plotting, I/O, validation
results/            Generated outputs (gitignored)
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
