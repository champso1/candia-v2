#!/usr/bin/sh

set -e

gnuplot ./p3nsm-approx-vs-exact-qns.gplt
gnuplot ./p3nsm-approx-vs-exact-lm.gplt
gnuplot ./p3nsp-approx-vs-exact-qns.gplt
gnuplot ./p3nsp-approx-vs-exact-lm.gplt
