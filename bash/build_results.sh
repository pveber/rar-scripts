#!/bin/bash

# produce annotations for chipseq regions
mkdir -p results/chipseq/annotations/
mkdir -p results/rnaseq
./rar-scripts/_build/annotations.native

R --vanilla < rar-scripts/R/heatmaps_of_jclusterings.R
R --vanilla < rar-scripts/R/profiles_of_jclusterings.R
R --vanilla < rar-scripts/R/conservation_in_jclusterings.R

mkdir -p results/chipseq/clustering_maison
R --vanilla < rar-scripts/R/chipseq_clustering_kmeans.R

mkdir -p results/rnaseq/clustering_maison
R --vanilla < rar-scripts/R/rnaseq_clustering_kmeans.R
