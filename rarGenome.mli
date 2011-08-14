open Batteries
open Pervasives
open Genome

(*
module Location : sig
  type t = private {
    chr : string ;
    st : int ;
    ed : int
  }

  val make : string -> int -> int -> t

  val to_string : t -> string 
(** String representation of a location, as <chr>:<start>-<end> *)
    
  val length : t -> int
    
  val of_string : string -> t
(** Parses a string representation of a genomic location, of the
    form <string>:<int>-<int>, like chr1:23-45. Reciprocal to the
    function {!to_string}.  *)

  val upstream : up:int -> down:int -> [`Sense | `Antisense] -> t -> t

  val included_in : t -> t -> bool

  val intersection : t -> t -> bool

  val inter : t -> t -> t
    
  val dist : t -> t -> int
(** Both locations should be on the same chromosome, throws
    [Invalid_argument "Ucsc.Location.dist"] otherwise *)

  val position : from:t -> t -> int
  val compare : t -> t -> int
end
*)

module Selection : sig
  type t = (string, ISet.t) PMap.t
  val of_locations : Location.t Enum.t -> t
  val inter : t -> t -> t
  val length : t -> int
end

