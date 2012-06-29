(*
 * for a fasta input under the convention that lower case means repeat, computes the proportion of repeats in each sequence and outputs it line by line
 *)

open Batteries
open Biocaml
open Printf

let islowercase = function
  | 'a'..'z' -> true
  | _ -> false

let repeat_proportion (_,seq) = 
  let n = 
    String.enum seq // islowercase |> Enum.count
  in float n /. float (String.length seq)


let () = 
  snd (Fasta.enum_of_file Sys.argv.(1))
  /@ repeat_proportion
  |> Enum.iter (printf "%f\n")

