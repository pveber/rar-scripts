#!/bin/bash

function die() {
    ECODE=$?
    echo !!! Error \#${ECODE} !!!
    echo -e "$1"
    exit ${ECODE}
}

# Asserts File Exists
function afe() {
    if [ ! -f $1 ]; then
	echo "Missing file $1";
	exit 1;
    fi
}



# Assumed input files
afe resources/chipseq/PanRAR_allregions_pmt.tsv
afe resources/chipseq/PanRAR_regions.bed
afe resources/chipseq/PanRAR_regions.fa

mkdir -p results/chipseq/annotations

# Detect spurious regions in the PanRAR selection
./scripts/ocaml/_build/repeat_proportion.native resources/chipseq/PanRAR_regions.fa > results/chipseq/annotations/PanRAR_repeats.tsv

# produce annotations for chipseq regions
mkdir -p results/chipseq/annotations/
mkdir -p results/rnaseq
./scripts/ocaml/_build/annotations.native 

./scripts/ocaml/_build/motif_prevalence.native resources/chipseq/regions/PanRAR_regions.fa results/chipseq/regions/PanRAR_regions_motif_prevalence
./scripts/ocaml/_build/motif_prevalence.native resources/chipseq/regions/Random_regions.fa results/chipseq/regions/Random_regions_motif_prevalence
./scripts/ocaml/_build/nhr_scan_fasta.native resources/chipseq/regions/PanRAR_regions.fa results/chipseq/regions/PanRAR_regions_nhr_scan
./scripts/ocaml/_build/nhr_scan_fasta.native resources/chipseq/regions/Random_regions.fa results/chipseq/regions/Random_regions_nhr_scan
./scripts/ocaml/_build/chipseq_track_annotation_main.native resources/chipseq/regions/PanRAR_regions.bed results/chipseq/annotations/PanRAR_regions_chipseq_es.tsv
./scripts/ocaml/_build/chipseq_track_annotation_main.native resources/chipseq/regions/Random_regions.bed results/chipseq/annotations/Random_regions_chipseq_es.tsv






# produce visualizations of justin clusterings for chipseq data
R --vanilla < scripts/R/heatmaps_of_jclusterings.R
R --vanilla < scripts/R/profiles_of_jclusterings.R
R --vanilla < scripts/R/conservation_in_jclusterings.R

# produce a naive clustering of chipseq data
mkdir -p results/chipseq/clustering_maison
R --vanilla < scripts/R/chipseq_clustering_kmeans.R

# mkdir -p results/rnaseq/clustering_maison
# R --vanilla < scripts/R/rnaseq_clustering_kmeans.R
