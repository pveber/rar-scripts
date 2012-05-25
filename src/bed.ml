open Batteries
open Biocaml

let enum ?(header = false) fn =
  let e = Tsv.enum fn in
  if header then Enum.drop 1 e ;
  e /@ (fun r -> r.(0), Range.make (int_of_string r.(1)) (int_of_string r.(2)))




















