# Reference results

Published Table 3 metrics are pinned in `configs/reference_metrics.m` (doi:10.3934/dcdsb.2026072).

Optional binary pins (after a trusted paper-scale `main('all')` on the author machine):

- `case_a_solution.mat`
- `case_b_solution.mat`
- `case_c_solution.mat`
- `figures/` — corresponding PNGs

Until those mats are archived here, use `main('validate')` for the numeric Table 3 self-check and compare new solves with `validate_against_reference`.

Do **not** treat `results/generated/` or legacy `results/case_*_solution.mat` as frozen references.
