(* This file is part of guizmin.

    guizmin is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    guizmin is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with guizmin.  If not, see <http://www.gnu.org/licenses/>.
*)
(** 
    Channel-based access to R statistical software 
    
    This interface is limited to sending commands, no
    reading from R is supported.
*)

(** Type of an R process *)
type t

val make : unit -> t

(** Send command to the process -- printf style *)
val c : t -> ('a,t,unit) format -> 'a

(** Defines a vector in the R process *)
val v : t -> string -> float array -> unit
val iv : t -> string -> int array -> unit
val sv : t -> string -> string array -> unit

val factor : t -> string -> string array -> unit

(** Defines a data frame in the R process *)
val frame : t -> string -> (string * [ `float of float array 
			             | `string of string array
				     | `factor of string array ]) list -> unit

val matrix : t -> string -> float array array -> unit

val pdf : t -> string -> unit
val devoff : t -> unit

(** 
    Should be called when dealing with the process is finished.
    Only returning from this call can you be sure that the 
    commands have been executed by the R process
*)
val close : t -> unit

