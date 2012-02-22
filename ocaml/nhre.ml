open Batteries
open Biocaml
open Printf

type motif = [`direct | `everted | `inverted] * int
type matches = (int * int * motif) list array

let string_of_motif = function
    `direct, n -> sprintf "DR%d" n
  | `everted, n -> sprintf "ER%d" n
  | `inverted, n -> sprintf "IR%d" n

let all_motifs = List.concat [
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


let counts_of_freqs n freqs = 
  Array.map 
    (Array.map (fun f -> int_of_float (float n *. f)))
    freqs

let balmer_counts = counts_of_freqs 309 balmer_freqs

let estimated_freqs = [|
  [| 0.776 ; 0.000 ; 0.224 ; 0.000 |] ;
  [| 0.000 ; 0.000 ; 1.000 ; 0.000 |] ;
  [| 0.000 ; 0.000 ; 0.707 ; 0.293 |] ;
  [| 0.000 ; 0.095 ; 0.066 ; 0.839 |] ;
  [| 0.000 ; 0.966 ; 0.034 ; 0.000 |] ;
  [| 1.000 ; 0.000 ; 0.000 ; 0.000 |] ;
|]


let estimated_freqs = [|
  [| 0.681 ; 0.020 ; 0.251 ; 0.049 |] ;
  [| 0.055 ; 0.000 ; 0.919 ; 0.026 |] ;
  [| 0.023 ; 0.039 ; 0.707 ; 0.231 |] ;
  [| 0.032 ; 0.075 ; 0.107 ; 0.786 |] ;
  [| 0.007 ; 0.837 ; 0.091 ; 0.065 |] ;
  [| 0.879 ; 0.046 ; 0.036 ; 0.039 |] ;
|]

let estimated_counts = counts_of_freqs 1000 estimated_freqs

let matrix hexamer_counts (orientation, spacer) = 
  let bg = Pwm.flat_background () in
  Pwm.tandem ~orientation ~spacer hexamer_counts hexamer_counts bg


let matches_of_sequence selected_motifs hexamer_counts theta s = 
  let annot motif sense (i,s) = (motif, sense, i, s) in
  selected_motifs
  |> List.map 
      List.(fun motif ->
	let mat = matrix hexamer_counts motif in
	let rcmat = Pwm.reverse_complement mat in 
	append
	  (Pwm.fast_scan mat s theta   |> map (annot motif `sense))
	  (Pwm.fast_scan rcmat s theta |> map (annot motif `antisense)))
  |> List.concat

let print_matrix m = 
  print_endline "  A       C      G     T" ;
  Array.iter
    (fun p -> 
      Array.iter (printf "%.3f  ") p ;
      print_newline ())
    m
