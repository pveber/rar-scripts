open Batteries
open Biocaml
open GenomeMap
open Guizmin_bioinfo.MBSchema

val closest : 
  ('a -> Location.t) ->
  (string, 'b) LMap.t ->
  'a Enum.t -> ('a * (Location.t * 'b)) Enum.t 
(** [closest f map a] returns an enumeration mapped from [a] where
    each member of [a] is annotated with the closest element in
    [map]. *)

val neighbours :
  ('a -> Location.t) ->
  int ->
  (string, 'b) LMap.t ->
  'a -> (Location.t * 'b) array
(** [neighbours f dmax map e] returns an annotation of [e] with all
    elements in [map] not further than [dmax] bp. *)

val score :
  (string * int) array ->
  ('a -> string location) -> 'a Enum.t ->
  ('b -> string location) -> 'b Enum.t ->
  ('a * 'b * float) Enum.t
(** [score fa a fb b] returns an enumeration of all pairs of a's and
    b's weighted by a score which reflects how usual the position of b
    with respect to a is, compared to a control *)

(** TODO 

val score_stranded_version :
   (string * int) array ->
   ('a -> string location * [`Sense | `Antisense]) -> 'a Enum.t ->
   ('b -> string location) -> 'b Enum.t ->
   ('a * 'b * float) Enum.t

*)
