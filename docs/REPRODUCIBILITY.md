# Reproducibility

Companion notes for *Optimal Chemo–Immunotherapy Scheduling: A Hybrid QPSO–SQP Approach with Michaelis–Menten Pharmacodynamics* (DCDS-B).

**DOI:** [10.3934/dcdsb.2026072](https://doi.org/10.3934/dcdsb.2026072)

## Environment

- MATLAB R2020b+ (tested R2024a / R2024b)
- Optimization Toolbox (`fmincon`, SQP)
- Parallel Computing Toolbox optional (`parfor` in QPSO; off by default)

Pinned settings: `configs/reproducibility_settings.m`.

## Publication defaults (critical)

Section 5 cases use $\ell=(0,0,0)$ and:

| Case | $(A,B,C)$ | Command |
|------|-----------|---------|
| A Balanced | $(1, 0.3, 1)$ | `main('a')` |
| B Asymmetric | $(1, 0, 1)$ | `main('b')` |
| C Min-time | $(0, 0, 1)$ | `main('c')` |

Table 1 biological parameters live in `configs/get_model_parameters.m`.  
Table 3 reference metrics live in `configs/reference_metrics.m`.

`assert_publication_defaults` runs on every paper/smoke solve and fails if Table 1 or case weights drift.

## ODE tolerances

- Paper Table 2 lists QPSO ODE `RelTol`/`AbsTol` $=10^{-3}$.
- Default reproduction path pins tighter author values (`1e-8` / `1e-10`) via `apply_reproducibility_settings` for rollout stability. The Table 2 pair is stored as `paper_table2_ode_*` in `reproducibility_settings.m`.

## Determinism

- `run_cases` calls `rng(42,'twister')` before solves.
- QPSO remains stochastic; identical MATLAB/toolbox builds help but bit-stability is not guaranteed.
- `useParfor` defaults to `false`.

## Reproduce paper cases

```matlab
cd hybrid-qpso-sqp-chemo-immuno
main            % Cases A, B, C + Table 3 check
main('a')       % Case A only
```

Outputs:

- `results/generated/case_*_solution.mat` (preferred)
- `results/case_*_solution.mat` (legacy mirror)
- `results/table3.csv`
- `results/figures/*.png`

## Reference pin

`configs/reference_metrics.m` encodes published Table 3:

| Metric | A | B | C |
|--------|---|---|---|
| $T$ | 4.917 | 4.723 | 4.417 |
| $D_u$ | 2.771 | 2.528 | 3.401 |
| $D_v$ | 1.785 | 4.723 | 4.417 |

`main('validate')` checks the pin without solving. Full paper-scale `.mat` archives under `results/reference/` are optional follow-up pins after a trusted `main('all')` run.

## Smoke / tests (not paper-quality)

```matlab
main('test')       % structural + publication pin checks
main('smoke')      % tiny QPSO/SQP plumbing on Case A
main('validate')   % Table 3 config pin only
```

## CI

`.github/workflows/repo-smoke.yml` checks layout and that Table 1 / Table 3 constants remain present. Full MATLAB solves are local (license required).
