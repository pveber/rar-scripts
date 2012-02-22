open Batteries
open Printf
open Biocaml
open Nhre
open Misc

let selected_motifs = [ `direct, 0 ; `direct, 1 ; `direct, 2 ; `direct, 5 ]

let fa =
  Fasta.enum_of_file Sys.argv.(1)
  |> snd
  |> List.of_enum

let sequences = List.map snd fa
let control_sequences = List.map string_shuffle sequences
let control_sequences = List.map (fun s -> random_dna_seq (String.length s)) (*string_shuffle*) sequences

let theta_estimation_step tol hexamer_counts = 
  let matches = 
    List.map (matches_of_sequence selected_motifs hexamer_counts 0.) sequences
  and control_matches = 
    List.map (matches_of_sequence selected_motifs hexamer_counts 0.) control_sequences
  in
  let theta = score_threshold tol control_matches in
  let _ = printf "theta = %f\ttreatment = %f\tcontrol = %f\n%!" theta (score_above theta matches) (score_above theta control_matches) in
  theta, matches, control_matches

let _ = 
  print_endline "Balmer threshold:" ;
  theta_estimation_step 0.1 Nhre.balmer_counts

let _ = 
  print_endline "Estimated counts threshold:" ;
  theta_estimation_step 0.1 Nhre.estimated_counts
