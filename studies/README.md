# Studies

Programs that produce benchmark tables and plots. Unlike the regression tests in `tests/`, they are not run automatically and have no pass/fail result: their output is meant to be looked at. They are built with `-DCANDIA_BUILD_STUDIES=ON`; `main` and `pdflumiunc` also need `-DCANDIA_WITH_LHAPDF=ON`. Each program is run from its directory in the build tree, where CMake copies the scripts and templates it needs.

## main

`evolve` evolves the Les Houches toy PDFs to Q = 100 GeV; `read` puts the resultant datafile into a formatted table; `compare` compares the results with e.g. hoppet or pegasus and formats the differences into a table; `precision` compares two candia outputs and records differences in a table; `ratio` outputs a new datafile with the ratios of two datafiles for each PDF and each x-val; `scaleratio` is the same as `read` but accepts also mu /= 1.0 and generates a composite table with all outputs.

| Program | Usage | Output |
|---|---|---|
| `evolve` | `<order> <iterations> <trunc_idx> <mur2_muf2> <use_trunc> <debug> [title]` | `data/<name>.dat`, `log/<name>.log` |
| `compare` | `<candia-file> <other-file> [<format> <type> <diff-type>]` | `diffs/diffs-other-<name>.pdf` |
| `read` | `<data-file> [<origin> <format> <type>]` | `tables/<name>.pdf` |
| `precision` | `<candia-file-1> <candia-file-2> <title>` | `diffs/diffs-candia-<title>.pdf` |
| `ratio` | `<base-file> <file1> <title1> [<fileN> <titleN> ...]` | `ratios/<title>.dat` |
| `scaleratio` | `<file-mu1> <file-mu0.5> <file-mu2> [<format> <type>] <title>` | `tables/<title>.pdf` |

Note that to generate any tables or plots, one needs `pdflatex` available on the command line.

There is also a backup script (`backup.sh [--remove-files]`) which archives the output directories.

## p3ns-exact

Compares N3LO evolutions with the exact non-singlet P3 against the approximations of arXiv:1707.08315 (central, IMOD 1 and IMOD 2), for P3nsm and P3nsp separately.

## plot-generation

`plot-generation <n3lo-file> <nnlo-file> <output>` writes the N3LO/NNLO ratio of b+, c+, s+, u_v, d_v and g from two `evolve` data files; `gnuplot -c plot-scripts/ratio.gplt <output> <pdf>` plots it. `gnuplot -c plot-scripts/main.gplt <evolve data file>` plots the valence, sea and strange distributions of one evolution.

## lhapdf-datafile

`lhapdf-datafile` writes the LHAPDF set `testpdf`. It also writes the charm subtraction PDFs at Q = 4 GeV to `out-c.dat` and the bottom ones at Q = 10 GeV to `out-b.dat`, which `gnuplot plot.gplt` plots.

## numerical-precision

Measures the numerical error of a single convolution E⊗A, for splitting functions from P0qq to P3gg and smooth test functions shaped like the toy gluon or u valence, at x = 1e-4 and 1e-2 (logarithmic region of the grid), 0.3 (linear) and 0.9 (quadratic). It compares three results with a high-precision reference that does not use Candia's mappings, nodes or grid:

- ε_GL: Candia's mappings and Gauss–Legendre nodes with the exact function at the nodes (quadrature error only);
- ε_total: the precomputed kernel acting on the function at the grid points, as in the evolution;
- ε_interp: the difference of the two, i.e. the interpolation error alone.

Sweep A varies `num_gauss_points` at the default grid (`sweepA-<kernel>.dat`); Sweep B scales all grid sizes (log, linear and quadratic parts together) at 200 Gauss–Legendre points (`sweepB-<kernel>.dat`) and prints the fitted slope of ε_interp against the grid spacing h. Run `numerical-precision` (no arguments, about a second), then `gnuplot plot.gplt` for `sweepA.pdf` and `sweepB.pdf` (one page per kernel).

## pdflumiunc

`evolvepdfset` evolves CT25NNLO at N3LO into the LHAPDF set `CT25aN3LO`, written to the working directory; `pdflumiunc`, run from the same directory, then compares its gg luminosity at √s = 13.6 TeV with MSHT20 aN3LO and NNPDF4.0 aN3LO (all these sets must be installed, of course), and `plot.gplt` plots the ratios.
