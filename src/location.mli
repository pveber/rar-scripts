open Biocaml

type 'a t = 'a * Range.t

val make : 'a -> int -> int -> 'a t
val size : 'a t -> int

val upstream : up:int -> down:int -> [`Sense | `Antisense] -> 'a t -> 'a t

(* (\* This file is part of guizmin. *)

(*     guizmin is free software: you can redistribute it and/or modify *)
(*     it under the terms of the GNU General Public License as published by *)
(*     the Free Software Foundation, either version 3 of the License, or *)
(*     (at your option) any later version. *)

(*     guizmin is distributed in the hope that it will be useful, *)
(*     but WITHOUT ANY WARRANTY; without even the implied warranty of *)
(*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the *)
(*     GNU General Public License for more details. *)

(*     You should have received a copy of the GNU General Public License *)
(*     along with guizmin.  If not, see <http://www.gnu.org/licenses/>. *)
(* *\) *)
(* open Batteries *)

(* type t = private { *)
(*   chr : string ; *)
(*   st : int ; *)
(*   ed : int *)
(* } *)
(* val move : t -> int -> int -> t *)
(* val relmove : t -> int -> int -> t *)



(* (\** [relative loc x a b] is the location obtained considering [a] and  *)
(*     [b] as relative coordinates inside [loc] *\) *)
(* val relative : t -> int -> int -> t *)

(* (\** [zoom l z] is a location of width [z *. float (Location.width l)] *\) *)
(* val zoom : t -> float -> t *)
  
(* (\** String representation of a location, as <chr>:<start>-<end> *\) *)
(* val to_string : t -> string  *)
(* val length : t -> int *)
  
(* (\** Parses a string representation of a genomic location, of the *)
(*     form <string>:int-int, like chr1:23-45. Reciprocal to the *)
(*     function {!to_string}.  *\) *)
(* val of_string : string -> t *)
(* val included_in : t -> t -> bool *)
(* val intersection : t -> t -> bool *)

(* val inter : t -> t -> t *)

(* val convex_hull : t list -> t *)
  
(* (\** Both locations should be on the same chromosome, throws *)
(*     [Invalid_argument "Ucsc.Location.dist"] otherwise *\) *)
(* val dist : t -> t -> int *)

(* val position : from:t -> t -> int *)
(* val compare : t -> t -> int *)
(* val enum_merge : t Enum.t -> t Enum.t -> t Enum.t *)
(* val islands : t Enum.t -> t Enum.t *)

(* type 'a annotation = 'a  *)
(* constraint 'a = < loc : t ; .. > *)
