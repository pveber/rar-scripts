open Batteries
open Biocaml
open Printf
open Misc

let selected_motifs = [ `direct, 0 ; `direct, 1 ; `direct, 2 ; `direct, 5 ]

let fa =
  Fasta.enum_of_file Sys.argv.(1)
  |> snd
  |> List.of_enum

let sequences = List.map snd fa
let control_sequences = List.map (fun s -> random_dna_seq (String.length s)) (*string_shuffle*) sequences

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
    List.map (Nhre.matches_of_sequence selected_motifs hexamer_counts 0.) sequences
  and control_matches = 
    List.map (Nhre.matches_of_sequence selected_motifs hexamer_counts 0.) control_sequences
  in
  let theta = score_threshold tol control_matches in
  let tpr = score_above theta matches
  and fpr = score_above theta control_matches in
  theta, tpr, fpr, matches, control_matches


let random_profile _ = 
  let a = Random.float 1. in 
  let c = Random.float (1. -. a) in
  let g = Random.float (1. -. a -. c) in
  let t = 1. -. a -. c -. g in
  [| a ; c ; g ; t |]

let random_hexamer_freqs () = 
  Array.init 6 random_profile

let matrix_copy t = Array.(map copy t)

let rec random_pair () =
  let k1, k2 = Random.int 4, Random.int 4 in
  if k1 = k2 then random_pair ()
  else k1, k2

let count_mutation t = 
  let r = matrix_copy t in
  let j = Random.int (Array.length r)
  and k1, k2 = random_pair () in
  let delta = min (Random.int 5) r.(j).(k1) in
  r.(j).(k1) <- r.(j).(k1) - delta ;
  r.(j).(k2) <- r.(j).(k2) + delta ;
  r

let count_mutation_enum t = 
  (Array.range t)
  /@ (fun j ->
    Enum.filter_map 
      (fun (k1, k2) -> 
	let delta = min (Random.int 5) t.(j).(k1) in
	if delta > 0 then 
	  let r = matrix_copy t in
	  r.(j).(k1) <- r.(j).(k1) - delta ;
	  r.(j).(k2) <- r.(j).(k2) + delta ;
	  Some r
	else None)
      (List.enum [ 0,1 ; 0,2 ; 0,3 ; 1,0 ; 1,2 ; 1,3 ; 2,0 ; 2,1 ; 2,3 ]))
  |> Enum.concat

let freq_matrix counts = 
  let profile x = 
    let n = Array.fold_left ( + ) 0 x in 
    Array.map (fun k -> float k /. float n) x
  in
  Array.map profile counts
  

let estimation tol hexamer_counts = 
  for k = 1 to 100 do 
    let theta, tpr, fpr, _, _ = theta_estimation_step tol hexamer_counts in
    let current = ref hexamer_counts 
    and opt = ref tpr in
    print_endline "!!! RESTART !!!" ;
    printf "new AUC is %f\n" !opt ;
    Nhre.print_matrix (freq_matrix !current) ;
    print_newline () ;
    for i = 1 to 10000 do
      let m = count_mutation !current in
      let theta, tpr, _, _, _ = theta_estimation_step tol m in
      if tpr > !opt then (
	current := m ;
	opt := tpr ;
	printf "new TPR is %f (theta = %f)\n" tpr theta ;
	Nhre.print_matrix (freq_matrix m) ;
	print_newline ()
      )
    done
  done



module PriorityQueue = struct
  type 'a t = (float * 'a) list

  let rec add_aux m x = function
    | [] -> [ x, m ]
    | (y,n) :: t when m = n -> (max x y, n) :: t
    | h :: t -> h :: (add_aux m x t)

  let add m x pq = 
    List.sort (flip compare) (add_aux m x pq)

end

      
let estimation2 tol hexamer_counts = 
  let score_init = Tuple5.second (theta_estimation_step tol hexamer_counts) in
  let queue = ref [ score_init , hexamer_counts ] 
  and seen = ref Set.empty 
  and opt = ref score_init 
  and current = ref hexamer_counts in
  for k = 1 to 100000 do
    printf "%d %d\n" (List.length !queue) (Set.cardinal !seen) ;
    let _, m = List.hd !queue in
    queue := List.tl !queue ;
    let theta, tpr, _, _, _ = theta_estimation_step tol m in
    seen := Set.add m !seen ;
    if tpr > !opt then (
      current := m ;
      opt := tpr ;
      printf "new TPR is %f (theta = %f)\n" tpr theta ;
      Nhre.print_matrix (freq_matrix m) ;
      print_newline () ;
    ) ;
    Enum.iter
      (fun n -> 
	if not (Set.mem n !seen) 
	then queue := PriorityQueue.add n tpr !queue)
      (count_mutation_enum m)
  done   

let _ = estimation 0.1 Nhre.estimated_counts
