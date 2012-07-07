open Batteries
open Guizmin_bioinfo.MBSchema

val uniform : 
  ?seed:int ->
  chrom_size:(string * int) array -> 
  ('b -> Location.t) ->
  'b array -> int -> Location.t Enum.t




















