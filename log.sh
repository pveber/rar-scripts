#!/bin/bash

# get justin's clusters
mkdir -p results/chipseq/clustering/justin
cd results/chipseq/clustering/justin
echo "get Justin's clusterings" && exit 0
tar xvfz clustall.tar.gz
cd -

# get the original panRAR regions selection
mkdir -p results/chipseq/regions/
echo "get original panRAR regions in panRAR_regions.bed" && exit 0

# produce annotations for chipseq regions
mkdir -p results/chipseq/annotations/

