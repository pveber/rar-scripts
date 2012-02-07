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

let print_matrix m = 
  print_endline "  A       C      G     T" ;
  Array.iter
    (fun p -> 
      Array.iter (printf "%.3f  ") p ;
      print_newline ())
    m

let nhre_matrix hexamer_counts (orientation, spacer) = 
  let bg = Pwm.flat_background () in
  Pwm.tandem ~orientation ~spacer hexamer_counts hexamer_counts bg

let fa =
  Fasta.enum_of_file Sys.argv.(1)
  |> snd
  |> List.of_enum

let sequences = List.map snd fa

(** randomly permute a string (taken from Core JAne St. *)
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

let control_sequences = List.map (fun s -> random_dna_seq (String.length s)) (*string_shuffle*) sequences

let matches_of_sequence hexamer_counts theta s = 
  let annot motif sense (i,s) = (motif, sense, i, s) in
  selected_motifs
  |> List.map 
      List.(fun motif ->
	let mat = nhre_matrix hexamer_counts motif in
	let rcmat = Pwm.reverse_complement mat in 
	append
	  (Pwm.fast_scan mat s theta   |> map (annot motif `sense))
	  (Pwm.fast_scan rcmat s theta |> map (annot motif `antisense)))
  |> List.concat
      
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
  |> quantile tol

let score_above theta matches =
  matches 
  |> List.map (List.fold_left (fun accu (_,_,_,s) -> max accu s) neg_infinity)
  |> List.enum
  |> proportion_above theta


let base_complement = function 
    'a' -> 't'
  | 't' -> 'a'
  | 'c' -> 'g'
  | 'g' -> 'c'
  | 'n' -> 'n'
  | 'A' -> 'T'
  | 'T' -> 'A'
  | 'C' -> 'G'
  | 'G' -> 'C'
  | 'N' -> 'N'
  | c -> printf "%c\n" c ; assert false

let reverse_complement s = 
  let n = String.length s in 
  let rc = String.make n ' ' in
  for i = 0 to n - 1 do
    rc.[i] <- base_complement s.[n - 1 - i]
  done ;
  rc

let int_of_nucleotide = function
  | 'a' | 'A' -> 0
  | 'c' | 'C' -> 1
  | 'g' | 'G' -> 2
  | 't' | 'T' -> 3
  | _ -> 4


let extract_hexamers (kind,spacer) sense seq i = 
  let sub sense i = 
    let s = String.sub seq i 6 in
    if sense = `sense then s else reverse_complement s
  in
  match kind with
      `direct -> [ sub sense i ; sub sense (i + 6 + spacer) ]
    | _ -> assert false (* not implemented *)

let compile_fragments fragments = 
  let r = Array.(init 6 (fun _ -> make 4 0)) in
  let foreach s =
    String.iteri 
      (fun i c -> 
	let k = int_of_nucleotide c in 
	r.(i).(k) <- r.(i).(k) + 1)
      s
  in
  List.iter foreach fragments ;
  r

let theta_estimation_step tol hexamer_counts = 
  let matches = 
    List.map (matches_of_sequence hexamer_counts 5.) sequences
  and control_matches = 
    List.map (matches_of_sequence hexamer_counts 5.) control_sequences
  in
  let theta = score_threshold (1. -. tol) control_matches in
  let _ = printf "theta = %f\ttreatment = %f\tcontrol = %f\n%!" theta (score_above theta matches) (score_above theta control_matches) in
  theta, matches, control_matches

let estimation_step tol hexamer_counts =
  let theta, matches, control_matches = theta_estimation_step tol hexamer_counts in
  let selected_matches = 
    List.map
      (List.filter (fun (_,_,_,x) -> x > theta))
      matches
  in
  let fragments = 
    List.map2
      (fun seq matches ->
	List.map (fun (motif,sense,i,_) -> extract_hexamers motif
	  sense seq i) matches |> List.concat)
      sequences selected_matches
    |> List.concat
  in
  compile_fragments fragments

let rec estimation n tol hexamer_counts = 
  if n <= 0 then Tuple3.first (theta_estimation_step tol hexamer_counts), hexamer_counts
  else estimation (n - 1) tol (estimation_step tol hexamer_counts)

let freq_matrix counts = 
  let profile x = 
    let n = Array.fold_left ( + ) 0 x in 
    Array.map (fun k -> float k /. float n) x
  in
  Array.map profile counts

let theta, estimated_counts = estimation 1 0.1 balmer_counts
let _ = print_matrix (freq_matrix estimated_counts)

let rec seq_enum n =
  if n <= 0 then [ "" ]
  else (
    let e = seq_enum (n - 1) in
    List.concat [
      List.map (fun s -> s ^ "A") e ;
      List.map (fun s -> s ^ "C") e ;
      List.map (fun s -> s ^ "G") e ;
      List.map (fun s -> s ^ "T") e ;
    ]
  )

(* Comptage du nombre de dodecamères qui passent le seuil trouvé *)
let () =
  let seqs = seq_enum 12 in
  let mat = nhre_matrix estimated_counts (`direct, 0) in
  let l =
    List.filter
      (fun s -> Pwm.fast_scan mat s theta <> [])
      seqs in
  printf "nb seqs compatible with theta = %f: %d\n" theta (List.length l)

let () = 
  let mat = nhre_matrix estimated_counts (`direct, 1) in
  List.iter
    (fun seq -> 
      Pwm.fast_scan mat seq 0. 
      |> List.iter (fun (_,x) -> printf "%s\t%f\n" seq x))
    [ 
      "AGGTCAgAGGTCA";
      "CGGTCAgAGGTCA";
      "ACGTCAgAGGTCA";
      "AGCTCAgAGGTCA";
      "AGGACAgAGGTCA";
      "AGGTAAgAGGTCA";
      "AGGTCCgAGGTCA";
    ]	  



