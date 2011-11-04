(** identify loci supporting protein binding across several conditions *)


open Oregon
open Genome
open Oregon.Target

val make : 
  summit:('a -> int) -> (** relative coordinate of the summit inside a peak *)
  location:('a -> Location.t) ->
  'a list ->            (** 'a is the type of peaks *)
  (string * int * int * int) list (** chromosome name, center of the region, width of the region, number of peaks involved *)

val of_macs_targets : 
  radius:int ->
  ?pvalue:float ->
  Macs.peak Bed.file list -> Location.t list value
