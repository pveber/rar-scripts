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

let estimated_freqs = [|
  [| 0.763 ; 0.006 ; 0.220 ; 0.011 |] ; 
  [| 0.012 ; 0.000 ; 0.976 ; 0.012 |] ;
  [| 0.033 ; 0.006 ; 0.660 ; 0.301 |] ;
  [| 0.019 ; 0.094 ; 0.073 ; 0.814 |] ;
  [| 0.005 ; 0.904 ; 0.060 ; 0.031 |] ;
  [| 0.950 ; 0.008 ; 0.023 ; 0.019 |] 
|]

