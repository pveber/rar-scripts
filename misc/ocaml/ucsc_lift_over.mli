open Batteries
open Genome
open Oregon

val conversion' : org_from:string -> org_to:string -> Location.t Enum.t -> Location.t Enum.t

val conversion : org_from:string -> org_to:string -> 'a Bed.file -> Bed.basic_file
