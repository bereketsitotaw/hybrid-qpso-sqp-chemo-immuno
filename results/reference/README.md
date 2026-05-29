# Frozen Reference Outputs

This directory stores **published reference metrics** and optionally **author-frozen simulation outputs** for validation.

## Static reference (always present)

| File | Description |
|------|-------------|
| `table3_paper.csv` | Table 3 values transcribed from the published paper |
| `table3_computed.csv` | Generated after running `reproduce_paper.m` |

## Optional frozen trajectories (author-generated)

After a successful reproduction run on the pinned author environment:

```matlab
reproduce_paper_pipeline('freeze_reference', true);
```

This creates:

| File | Contents |
|------|----------|
| `canonical_outputs.mat` | Full state/control trajectories + summaries for Cases A/B/C |
| `freeze_metadata.mat` | MATLAB version, reproducibility settings, timestamp |

## Validation

```matlab
bootstrap_project();
load('results/case_a_solution.mat'); % etc.
summaries = struct('case_a', summary, ...);
validate_against_reference(summaries, 'verbose', true);
```

Or simply:

```matlab
run('scripts/reproduce_paper.m')
```

Validation compares `T`, `Du`, `Dv`, and `x(T)` against `configs/reference_metrics.m` within tolerances in `configs/reproducibility_settings.m`.

## Note on stochastic QPSO

Small deviations in terminal time and drug exposure are expected if MATLAB version, RNG seed, or ODE tolerances differ from the pinned author settings documented in `README.md`.
