open Batteries
open Biocaml
open GenomeMap

val closest : 
  ('a -> 'chr location) ->
  ('chr, 'b) LMap.t ->
  'a Enum.t -> ('a * ('chr location * 'b)) Enum.t 

val score :
  ('chr * int) list ->
  ('a -> 'chr location) ->
  'a Enum.t ->
  ('b -> 'chr location) ->
  'b Enum.t ->
  ('a * 'b * float) Enum.t





















