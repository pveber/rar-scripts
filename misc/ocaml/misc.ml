open Batteries

(** randomly permute a string (taken from Core Jane St. *)
let string_shuffle s =
  let s = String.copy s in
  let swap i j =
    let tmp = s.[i] in
    s.[i] <- s.[j];
    s.[j] <- tmp 
  in
  for i = String.length s downto 2 do
    swap (i - 1) (Random.int i)
  done ;
  s

let random_dna_seq n = 
  String.init 
    n 
    (fun i -> match Random.int 4 with 0 -> 'A' | 1 -> 'C' | 2 -> 'G' | 3 -> 'T' | _ -> assert false)

let quantile tol xs = 
  let open Array in
  let t = of_enum xs in
  sort compare t ;
  let n = length t in
  t.(int_of_float (float n *. tol))

let proportion_above theta xs = 
  let open Array in
  let t = of_enum xs in
  (filter (( < ) theta) t |> length |> float) /. (length t |> float)

let score_threshold tol matches =
  matches 
  |> List.map (List.fold_left (fun accu (_,_,_,s) -> max accu s) neg_infinity)
  |> List.enum
  |> quantile (1. -. tol)

let score_above theta matches =
  matches 
  |> List.map (List.fold_left (fun accu (_,_,_,s) -> max accu s) neg_infinity)
  |> List.enum
  |> proportion_above theta
