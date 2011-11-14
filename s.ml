open Batteries
open Oregon
open Oregon.Target
open Oregon.Target.Infix 
open Sle

let less ?(p = []) x = sh "less %s" (String.concat "/" (x#path :: p))

let evince pdf = sh "evince %s" pdf#path

let bed_to_value bed = 
  V.make (object
    method id = "S.bed_to_value"
    method deps = [] ++ bed
    method build = 
      Tsv.enum bed /@ (fun x -> x#loc) |> Array.of_enum
  end)

	   
