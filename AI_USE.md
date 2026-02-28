# AI Use Acknowledgment

The code and documentation in the `lpt` package were created with the
assistance of **Claude**, an AI assistant developed by
[Anthropic](https://www.anthropic.com).

## Scope

AI assistance was used throughout the development of this package,
including:

- **Package architecture** — design of the S3 class structure, function
  interfaces, and internal data flow
- **R code** — implementation of estimation functions (`lpt()`,
  `calibrate_B()`, `estimate_dose_slope()`), plot helpers, and S3
  methods (`print`, `summary`, `plot`)
- **Documentation** — roxygen2 docstrings, the package vignette, and
  this README
- **Tests** — unit test scaffolding in `tests/testthat/`

## Verification

All AI-generated code was reviewed, tested, and validated by the package
author against the econometric methodology described in:

> Dealberti, F. (2026). *Lipschitz Parallel Trends for Continuous
> Difference-in-Differences.*

The theoretical content, identification results, and empirical application
to the SRU dataset are the author's original work. AI assistance was
limited to translating methodology into R code and documentation.

## Model

- **AI system:** Claude (claude-sonnet-4-6) by Anthropic
- **Interface:** Claude Code (CLI)
