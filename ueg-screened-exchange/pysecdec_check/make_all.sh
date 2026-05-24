#!/bin/sh
# Compile the generated pySecDec packages (a few minutes with -j4).
set -e
cd "$(dirname "$0")"
for d in fp_cross me_euler; do make -C "$d" -j4; done
