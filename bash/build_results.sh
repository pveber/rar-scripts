#!/bin/bash

# produce annotations for chipseq regions
mkdir -p results/chipseq/annotations/
./rar-scripts/_build/annotations.native
R --vanilla < rar-scripts/R/conservation_in_jclusterings.R
