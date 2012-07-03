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
  ('a -> 'chr location) ->
  int ->
  ('chr, 'b) LMap.t ->
  'a -> ('chr location * 'b) array
(** [neighbours f dmax map e] returns an annotation of [e] with all
    elements in [map] not further than [dmax] bp. *)


val range_pos : from:Range.t -> Range.t -> int
val stranded_range_pos : from:(Range.t * [`Sense | `Antisense]) -> Range.t -> int

val score :
  ('chr * int) array ->
  ('a -> 'chr location) -> 'a Enum.t ->
  ('b -> 'chr location) -> 'b Enum.t ->
  ('a * 'b * float) Enum.t
(** [score fa a fb b] returns an enumeration of all pairs of a's and
    b's weighted by a score which reflects how usual the position of b
    with respect to a is, compared to a control *)

(**TODO 

val score_stranded_version :
   ('chr * int) array ->
   ('a -> 'chr location * [`Sense | `Antisense]) -> 'a Enum.t ->
   ('b -> 'chr location) -> 'b Enum.t ->
   ('a * 'b * float) Enum.t

*)






















