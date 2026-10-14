# Studies

Programs that produce benchmark tables and plots. Unlike the regression tests in `tests/`, they are not run automatically and have no pass/fail result: their output is meant to be looked at. They are built with `-DCANDIA_BUILD_STUDIES=ON`; `main`, `lhapdf-datafile` and `pdflumiunc` also need `-DCANDIA_WITH_LHAPDF=ON`. Each program is run from its directory in the build tree, where CMake copies the scripts and templates it needs.

## main

Evolves the Les Houches toy PDFs to Q = 100 GeV and compares the result with the hoppet and pegasus benchmark tables (`data/`, copied to `hoppet/` and `pegasus/` in the build directory). The tables are LaTeX documents compiled with `pdflatex`.

| Program | Usage | Output |
|---|---|---|
| `evolve` | `<order> <iterations> <trunc_idx> <mur2_muf2> <use_trunc> <debug> [title]` | `data/<name>.dat`, `log/<name>.log` |
| `compare` | `<candia-file> <other-file> [<format> <type> <diff-type>]` | `diffs/diffs-other-<name>.pdf` |
| `read` | `<data-file> [<origin> <format> <type>]` | `tables/<name>.pdf` |
| `precision` | `<candia-file-1> <candia-file-2> <title>` | `diffs/diffs-candia-<title>.pdf` |
| `ratio` | `<base-file> <file1> <title1> [<fileN> <titleN> ...]` | `ratios/<title>.dat` |
| `scaleratio` | `<file-mu1> <file-mu0.5> <file-mu2> [<format> <type>] <title>` | `tables/<title>.pdf` |

Running a program without arguments prints what the options mean. `evolve` uses the approximate (ImodAvg) non-singlet P3, as in the hoppet N3LO table; other P3 variants (e.g. for `plot-scripts/approx-vs-exact.gplt`) are made by changing its `setP3ApproximationTypes` call. `backup.sh [--remove-files]` archives the output directories.

## p3ns-exact

Compares N3LO evolutions with the exact non-singlet P3 against the approximations of arXiv:1707.08315 (central, IMOD 1 and IMOD 2), for P3nsm and P3nsp separately. Run `p3ns-exact` (no arguments), then `plot.sh` for the PDFs.

## plot-generation

`plot-generation <n3lo-file> <nnlo-file> <output>` writes the N3LO/NNLO ratio of b+, c+, s+, u_v, d_v and g from two `evolve` data files; `gnuplot -c plot-scripts/ratio.gplt <output> <pdf>` plots it.

## lhapdf-datafile

`lhapdf-datafile` writes an N3LO LHAPDF grid of the toy PDFs; `subpdfs` computes the charm subtraction PDFs from it; `plot.gplt` plots them.

## pdflumiunc

`evolvepdfset` evolves CT25NNLO at N3LO into the LHAPDF set `CT25aN3LO`; `pdflumiunc` then compares its gg luminosity at √s = 13.6 TeV with MSHT20 aN3LO and NNPDF4.0 aN3LO (these sets must be installed), and `plot.gplt` plots the ratio.
