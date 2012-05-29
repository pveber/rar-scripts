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
open Printf
open Sle

let option f = function
    None -> "."
  | Some x -> f x

let string_of_strand = function
    `plus -> "+"
  | `minus -> "-"

let format_attrib = function
| `S s -> Sle.string_quote s
| `I i -> string_of_int i
| `F f -> sprintf "%g" f 
    
let unparse_line ~seqname ~source ~feature ~start ~end_
    ?(score = None) ?(strand = None) ?(frame = None)
    ~group oc =
  Tsv.unparse_row oc [| seqname ; 
			source ;
			feature ;
			string_of_int start ;
			string_of_int end_ ;
			option string_of_int score ;
			option string_of_strand strand ;
			option string_of_int frame ;
			List.map (fun (k,v) -> sprintf "%s %s" k (format_attrib v)) group 
                        |> String.concat "; " |]


type item = {
  loc : string Location.t ;
  strand : [`plus | `minus] option ;
  kind : [ `exon | `gene | `CDS | `start_codon | `stop_codon ] ;
  attr : (string * string) list
}

let pair_of_array = function
    [| x ; y |] -> x, Sle.string_unquote y
  | _ -> raise (Invalid_argument "Gtf.pair_of_array")

let parse_attributes attr_field = 
  let fields = String.nsplit attr_field "; " in
    List.map (fun s -> pair_of_array (Sle.string_split_noeps ' ' s)) fields

let kind_of_string = function
    "exon" -> `exon
  | "gene" -> `gene
  | "CDS"  -> `CDS
  | "start_codon" -> `start_codon
  | "stop_codon" -> `stop_codon
  | s -> (print_endline s ;assert false)
    
let parse_line chr col = {
  loc = Location.make (chr col.(0)) (int_of_string col.(3)) (int_of_string col.(4)) ;
  strand = (match col.(6) with "+" -> Some `plus | "-" -> Some `minus | _ -> None) ;
  kind = kind_of_string col.(2) ;
  attr = parse_attributes col.(8)
}

let gene_id item = List.assoc "gene_id" item.attr
let transcript_id item = List.assoc "transcript_id" item.attr

let enum ?(chr = identity) fn =
  Tsv.enum fn /@ parse_line chr


(* let items gtf =  *)
(*   let accu = Hashtbl.create 1024 in *)
(*   let aux col () =  *)
(*     let item = parse_line col in *)
(*       accu <+- ((gene_id item, transcript_id item), item) *)
(*   in *)
(*     Tsv.fold_lines aux gtf () ; *)
(*     accu *)

(* let transcript_of_items f _ items = *)
(*   let items = List.filter (fun i -> i.kind = `exon) items in *)
(*   let items = List.sort (fun x y -> Location.compare x.loc y.loc) items in *)
(*     { Transcript.annot = f items ; *)
(*       Transcript.exons = List.map (fun i -> i.loc) items } *)

(* let transcripts f gtf =  *)
(*   let items =  items gtf in  *)
(*     Hashtbl.map *)
(*       (transcript_of_items f) *)
(*       items *)

(* open Sle.Infix *)

(* let items gtf =  *)
(*   let accu = Hashtbl.create 1024 in *)
(*   let aux item () = accu <+- ((gene_id item, transcript_id item), item) *)
(*   in *)
(*     Tsv.fold aux gtf () ; *)
(*     accu *)

(* let transcript_of_items f _ items = *)
(*   let items = List.filter (fun i -> i.kind = `exon) items in *)
(*   let items = List.sort (fun x y -> Location.compare x.loc y.loc) items in *)
(*     { Transcript.annot = f items ; *)
(*       Transcript.exons = List.map (fun i -> i.loc) items } *)

(* let transcripts f gtf =  *)
(*   let items = items gtf in  *)
(*   Hashtbl.map *)
(*     (transcript_of_items f) *)
(*     items *)
