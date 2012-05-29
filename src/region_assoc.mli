open Batteries
open Biocaml
open GenomeMap

val closest : 
  ('a -> 'chr location) ->
  ('chr, 'b) LMap.t ->
  'a Enum.t -> ('a * ('chr location * 'b)) Enum.t 
(** [closest f map a] returns an enumeration mapped from [a] where
    each member of [a] is annotated with the closest element in
    [map]. *)

val range_pos : from:Range.t -> Range.t -> int

(* val score : *)
(*   ('chr * int) list -> *)
(*   ('a -> 'chr location) -> *)
(*   'a Enum.t -> *)
(*   ('b -> 'chr location) -> *)
(*   'b Enum.t -> *)
(*   ('a * 'b * float) Enum.t *)

val score :
  (string * int) list ->
  ('a -> string location) ->
  'a Enum.t ->
  ('b -> string location) ->
  'b Enum.t ->
  ('a * 'b * float) Enum.t





















