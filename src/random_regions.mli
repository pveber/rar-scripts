open Batteries

val uniform : 
  ?seed:int ->
  chrom_size:('a * int) array -> 
  ('b -> 'a Location.t) ->
  'b array -> int -> 'a Location.t Enum.t




















