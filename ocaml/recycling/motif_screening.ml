open Batteries
open Printf
open Oregon

let rank l =
  List.sort ~cmp:(fun (_,x) (_,y) -> compare y x) l

let output_results oc l =
  let ranking = rank l in
  List.iter
    (fun (item,auc) -> 
       fprintf oc "%.3f\t\t\t%s\n" auc (Motif_library.string_of_item item))
    (List.take 40 ranking)
       

let interactive library bed =
  let fa = Ucsc.fasta_of_bed `mm9 bed in
  let ctrl = Fasta.shuffle fa in
  let res = ref [] in
  Array.iter
    (fun item -> 
       let auc = Motif_library.roc_auc (Motif_library.eval_motif' item fa#path ctrl#path) in
       res := (item, auc) :: !res ;
       output_results stdout !res ;
       printf "%s just gave %g\n%!" (Motif_library.string_of_item item) auc)
    library ;
  rank !res

