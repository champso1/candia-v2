# Studies

Programs that produce benchmark tables and plots. Unlike the regression tests in `tests/`, they are not automatically run or anything, and are meant to be either one-off codes for table generation, or in the case of `main` a full set of evolution and table/plot generation utilities. They are built whenever `-DCANDIA_BUILD_STUDIES=ON`; `main`, `lhapdf-datafile` and `pdflumiunc` need LHAPDF.

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

## pdflumiunc

`evolvepdfset` evolves CT25NNLO at N3LO into the LHAPDF set `CT25aN3LO`, written to the working directory; `pdflumiunc`, run from the same directory, then compares its gg luminosity at √s = 13.6 TeV with MSHT20 aN3LO and NNPDF4.0 aN3LO (all these sets must be installed, of course), and `plot.gplt` plots the ratios.
