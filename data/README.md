# Benchmark data

Reference tables used by the regression tests (`tests/benchmark.cpp`) and the study programs (`studies/main`).

- `hoppet/`, `pegasus/`: the Les Houches toy PDFs evolved to Q = 100 GeV with hoppet and pegasus, at the standard x values (columns: x, g, u, d, s, c, b, t, ū, d̄, s̄, c̄, b̄). `-mu05`/`-mu20` are μ_R²/μ_F² = 0.5 and 2. The hoppet N3LO tables use the approximate non-singlet P3.
- `hoppet/*-ffn4-*`: the same in the fixed flavour number scheme with n_f = 4, made with `hoppet/lh_hoppet.cpp` (hoppet 1.2.0); the regression tests use them for μ_R ≠ μ_F.

The hoppet variable-flavour tables with μ_R ≠ μ_F (`out-nlo-mu05/20`, `out-n3lo-mu05/20`) follow hoppet's own treatment of μ_R ≠ μ_F across the heavy-quark thresholds, which differs from the Les Houches benchmark prescription used by pegasus (and Candia) by about 1%, so the tests compare those cases with the pegasus tables instead.
