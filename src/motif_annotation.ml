open Batteries

let quantile tol xs = 
  let open Array in
  let t = of_enum xs in
  sort Pervasives.compare t ;
  let n = length t in
  t.(int_of_float (float n *. tol))

let score_threshold tol matches =
  matches 
  |> List.map (List.fold_left (fun accu (_,_,_,s) -> max accu s) neg_infinity)
  |> List.enum
  |> quantile (1. -. tol)

let random_dna_seq n = 
  String.init 
    n 
    (fun i -> match Random.int 4 with 0 -> 'A' | 1 -> 'C' | 2 -> 'G' | 3 -> 'T' | _ -> assert false)

let threshold_of_sequences tol sequences =
  let control_sequences = List.map (fun s -> random_dna_seq (String.length s)) sequences in 
  let control_matches : 'a list list = 
    List.map (Nhre.scan_string [`direct, 0] 0.) control_sequences
  in
  score_threshold tol control_matches


let of_sequences motifs tol sequences = 
  let sequences = List.of_enum sequences in
  let theta = threshold_of_sequences tol sequences in 
  List.enum sequences
  /@ Nhre.scan_string motifs theta




















