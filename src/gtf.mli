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
open Batteries
open Sle
open Biocaml


type item = {
  loc : Location.t ;
  strand : [`plus | `minus] option ;
  kind : [ `exon | `gene | `CDS | `start_codon | `stop_codon ] ;
  attr : (string * string) list
}

val unparse_line : 
  seqname:string -> source:string -> feature:string -> 
  start:int -> end_:int -> 
  ?score:int option -> ?strand:[`plus | `minus] option ->
  ?frame:int option -> 
  group:(string * [`S of string | `I of int | `F of float]) list ->
  unit IO.output -> unit

val gene_id : item -> string
val transcript_id : item -> string

val parse_line : (string -> string) -> string array -> item 
(** first argument is for renaming chromosome id  *)

val enum : ?chr:(string -> string) -> string -> item Enum.t


(* val items_by_id : file -> (string * string, item) hrel *)

(* (\** Extract exon items of a gff and sort them in a hash table by gene id and transcript id *\) *)
(* val transcripts : (item list -> 'a) -> file -> (string * string, 'a Transcript.t) hfun *)

(* val target_type : item Tsv.ty *)













