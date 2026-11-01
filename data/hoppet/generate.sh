#!/bin/bash
# Regenerates every hoppet table in this directory with lh_hoppet.cpp.
# usage: generate.sh <path to hoppet-config of a hoppet-v2 built with HOPPET_BUILD_EXACT_N3LO_SPLIT=ON>
# e.g.   data/hoppet/generate.sh ~/software/hoppet/install/bin/hoppet-config

set -euo pipefail

hoppet_config=${1:?usage: generate.sh <path to hoppet-config>}
here=$(cd "$(dirname "$0")" && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

g++ -O2 $("$hoppet_config" --cxxflags) "$here/lh_hoppet.cpp" -o "$work/lh_hoppet" $("$hoppet_config" --libs) -lgfortran
# hoppet prints a banner (and expected warnings) on every run
run() { "$work/lh_hoppet" "$@" > /dev/null 2>&1; }

# variable flavour number scheme, mu_R = mu_F
run 1 1.0 0 approx "$here/out-lo.dat"
run 2 1.0 0 approx "$here/out-nlo.dat"
run 3 1.0 0 approx "$here/out-nnlo.dat"
run 4 1.0 0 approx "$here/out-n3lo-approx.dat"
run 4 1.0 0 exact  "$here/out-n3lo-exact.dat"

# variable flavour number scheme, mu_R != mu_F (hoppet's own treatment across the thresholds, see README.md)
for mu in 05 20; do
	ratio=$([[ $mu == 05 ]] && echo 0.5 || echo 2.0)
	run 2 $ratio 0 approx "$here/out-nlo-mu$mu.dat"
	run 4 $ratio 0 approx "$here/out-n3lo-approx-mu$mu.dat"
done

# fixed flavour number scheme (n_f = 4), mu_R != mu_F
for mu in 05 20; do
	ratio=$([[ $mu == 05 ]] && echo 0.5 || echo 2.0)
	run 1 $ratio 4 approx "$here/out-lo-ffn4-mu$mu.dat"
	run 2 $ratio 4 approx "$here/out-nlo-ffn4-mu$mu.dat"
	run 3 $ratio 4 approx "$here/out-nnlo-ffn4-mu$mu.dat"
	run 4 $ratio 4 approx "$here/out-n3lo-approx-ffn4-mu$mu.dat"
	run 4 $ratio 4 exact  "$here/out-n3lo-exact-ffn4-mu$mu.dat"
done
