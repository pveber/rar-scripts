open Batteries
open Biocaml
open Printf

type motif = [`direct | `everted | `inverted] * int
type matches = (int * int * motif) list array

let string_of_motif = function
    `direct, n -> sprintf "DR%d" n
  | `everted, n -> sprintf "ER%d" n
  | `inverted, n -> sprintf "IR%d" n

let selected_motifs = List.concat [
  List.init 10 (fun i -> `direct, i) ;
  List.init 10 (fun i -> `everted, i) ;
  List.init 10 (fun i -> `inverted, i) ;
]

let balmer_freqs = [|
  [| 0.654 ; 0.045 ; 0.262 ; 0.039 |] ;
  [| 0.019 ; 0.01  ; 0.935 ; 0.036 |] ;
  [| 0.042 ; 0.013 ; 0.673 ; 0.272 |] ;
  [| 0.013 ; 0.074 ; 0.133 ; 0.78  |] ;
  [| 0.01  ; 0.819 ; 0.113 ; 0.058 |] ;
  [| 0.893 ; 0.01  ; 0.068 ; 0.029 |]
|]

let balmer_counts = 
  Array.map 
    (Array.map (fun f -> int_of_float (float 309 *. f)))
    balmer_freqs

let nhre_matrix (orientation, spacer) seq = 
  let bg = Pwm.background_of_sequence seq 0.1 in
  Pwm.tandem ~orientation ~spacer balmer_counts balmer_counts bg

let fa =
  Fasta.enum_of_file Sys.argv.(1)
  |> snd
  |> List.of_enum

let sequences = List.map snd fa |> Array.of_list

let matches_of_sequence s = 
  List.map 
    (fun motif ->
       let mat = nhre_matrix motif s in
       let rcmat = Pwm.reverse_complement mat in 
       List.append
	 (Pwm.fast_scan mat s 5.)
	 (Pwm.fast_scan rcmat s 5.)
       |> List.map (fun (i,s) -> motif, i, s))
    selected_motifs
  |> List.concat
      
let matches = 
  Array.map matches_of_sequence sequences

