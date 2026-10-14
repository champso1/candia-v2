# Benchmark data

Reference tables used by the regression tests (`tests/benchmark.cpp`) and the study programs (`studies/main`).

- `hoppet/`, `pegasus/`: the Les Houches toy PDFs evolved to Q = 100 GeV with hoppet and pegasus, at the standard x values (columns: x, g, u, d, s, c, b, t, ū, d̄, s̄, c̄, b̄). `-mu05`/`-mu20` are μ_R²/μ_F² = 0.5 and 2. The hoppet N3LO table uses the approximate non-singlet P3.
- `fhmruvv/`: `FHMRUVV_VFNS.dat`, N3LO benchmark values in the variable-flavour scheme (columns: Q, x, xu_v, xd_v, xL-, 2xL+, xs+, xc+, xb+, xg). Not used yet.
