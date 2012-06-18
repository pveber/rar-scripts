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

val neighbours :
  ('a -> string location) ->
  int ->
  (string, 'b) LMap.t ->
  'a Enum.t -> ('a * (string location * 'b) array) Enum.t 
(** [neighbours f dmax map a] returns an enumeration mapped from [a]
    where each membre of [a] is annotated with all elements in [map]
    not further than [dmax] bp. *)


val range_pos : from:Range.t -> Range.t -> int

val score :
  ('chr * int) array ->
  ('a -> 'chr location) ->
  'a Enum.t ->
  ('b -> 'chr location) ->
  'b Enum.t ->
  ('a * 'b * float) Enum.t

    







