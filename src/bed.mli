open Batteries
open Biocaml

val enum : ?header:bool -> string -> (string * Range.t) Enum.t
(** [header] default value is [false] *)
