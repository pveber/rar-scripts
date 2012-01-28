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

let selected_motifs = [ `direct, 0 ; `direct, 1 ; `direct, 2 ; `direct, 5 ]

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

let nhre_matrix hexamer_counts (orientation, spacer) seq = 
  let bg = Pwm.background_of_sequence seq 0.1 in
  Pwm.tandem ~orientation ~spacer hexamer_counts hexamer_counts bg

let fa =
  Fasta.enum_of_file Sys.argv.(1)
  |> snd
  |> List.of_enum

let sequences = List.map snd fa
let control_sequences = List.map (* FIXME *) identity sequences

let matches_of_sequence hexamer_counts s = 
  let annot motif sense (i,s) -> (motif, sense, i, s) in
  selected_motifs
  |> List.map 
      (fun motif ->
	 let mat = nhre_matrix hexamer_counts motif s in
	 let rcmat = Pwm.reverse_complement mat in 
	 List.(
	   append
	     (Pwm.fast_scan mat s 5.   |> map (annot motif `sense))
	     (Pwm.fast_scan rcmat s 5. |> map (annot motif `antisense))))
  |> List.concat
      
let quantile tol xs = assert false

let score_threshold matches tol =
  matches 
  |> List.map (List.map (fun (_,_,s) -> s))
  |> List.concat
  |> quantile tol

let rec estimation hexamer_counts =
  let matches = 
    List.map (matches_of_sequence hexamer_counts) sequences
  and control_matches = 
    List.map (matches_of_sequence hexamer_counts) control_sequences
  in
  let theta = score_threshold control_matches 0.1 in
  let selected_matches = 
    List.map
      (List.filter (fun (_,_,x) -> x > theta))
      matches
  in
  let fragments = 
    List.map2
      (fun seq matches ->
	List.map (fun (motif,i,_) -> extract motif seq i) matches)
      sequences selected_matches
    |> List.concat
    
