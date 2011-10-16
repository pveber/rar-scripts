open Batteries
open Printf
open Oregon
open Archive.Infix

let motif_rank library bed = 
  let peakseq = Ucsc.fasta_of_bed `mm9 bed in 
  let ctrlseq = Fasta.enumerate (Fasta.shuffle peakseq) in
  Motif_library.rank library peakseq ctrlseq

let motif_roc_curve item bed = 
  let peakseq = Ucsc.fasta_of_bed `mm9 bed in 
  let ctrlseq = Fasta.enumerate (Fasta.shuffle peakseq) in
  Motif_library.roc_curve item peakseq ctrlseq
  
let rank_on_selected_regions = motif_rank Selected_motifs.drerir Chipseq_regions.selected_regions_bed

let rank_for_time t = Chipseq_regions.(
  sprintf "rank_on_selected_regions_at_%s.pdf" (string_of_time t), 
  motif_rank Selected_motifs.drerir (Chipseq_regions.selected_regions_bed_for_time t)
)

let archive = Archive.make (
  [] +| ("rank_on_selected_regions.pdf", rank_on_selected_regions)
     +| ("roc_curve_DR0.pdf", motif_roc_curve (Selected_motifs.balmer_drn 0) Chipseq_regions.selected_regions_bed)
     +| ("roc_curve_DR5.pdf", motif_roc_curve (Selected_motifs.balmer_drn 5) Chipseq_regions.selected_regions_bed)
     +| ("roc_curve_DR2.pdf", motif_roc_curve (Selected_motifs.balmer_drn 2) Chipseq_regions.selected_regions_bed)
     +| ("roc_curve_DR1.pdf", motif_roc_curve (Selected_motifs.balmer_drn 1) Chipseq_regions.selected_regions_bed)
     +| ("roc_curve_IR3.pdf", motif_roc_curve (Selected_motifs.balmer_irn 3) Chipseq_regions.selected_regions_bed)
     +| ("roc_curve_DR4.pdf", motif_roc_curve (Selected_motifs.balmer_drn 4) Chipseq_regions.selected_regions_bed)
     +| rank_for_time `t0
     +| rank_for_time `t2
     +| rank_for_time `t24
     +| rank_for_time `t48
)
