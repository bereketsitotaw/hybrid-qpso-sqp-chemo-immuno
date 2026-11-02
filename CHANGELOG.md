# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added
- `CHANGELOG.md`
- `tests/` smoke suite (`main('test')`)
- `docs/REPRODUCIBILITY.md`
- `.github/workflows/repo-smoke.yml` layout/publication-pin CI
- `configs/assert_publication_defaults.m`, `configs/apply_smoke_overrides.m`
- `main('smoke')`, `main('validate')` entry points
- `results/generated/` and `results/reference/` layout notes
- SQP exitflag / constraint-violation diagnostics in Stage 2

### Changed
- Expanded `.gitignore` for MATLAB / OS artifacts and generated results under `results/generated/`
- `run_cases` writes preferred generated mats + legacy mirror; asserts publication defaults
- README expanded with method overview, tests/smoke, reproducibility link

### Fixed
- (none — publication Table 1 / case weights already matched DCDS-B; no scientific restoration required)
