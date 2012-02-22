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

type t = unit IO.output

let make () = 
  Unix.open_process_out "tee delme.R | R --no-save > /dev/null"

let c oc fmt = kfprintf (fun oc -> fprintf oc "\n%!") oc fmt 

let format_sequence f oc data = 
  match Enum.get data with
      None -> ()
    | Some h -> 
	f oc h ;
	Enum.iter
	  (fun x -> output_char oc ',' ; f oc x)
	  data

let _c f oc data = 
  output_string oc "c(" ;
  format_sequence f oc (Array.enum data) ;
  output_string oc ")"

let format_float oc x = 
  if classify_float x = FP_nan
  then output_string oc "NA"
  else output_string oc (string_of_float x)

let bind f oc ident value = 
  fprintf oc "%s <- " ident ;
  f oc value ;
  output_string oc "\n"

let format_int oc x = output_string oc (string_of_int x)

let v  = bind (_c format_float)
let iv = bind (_c format_int)

let format_string oc s = 
  output_char oc '"' ;
  output_string oc s ;
  output_char oc '"'

let _factor oc data = 
  output_string oc "factor(" ;
  _c format_string oc data ;
  output_string oc ")"


let sv = bind (_c format_string)

let factor = bind _factor

let format_frame_item oc (label,value) =
  fprintf oc "%s=" label ;
  (match value with 
       `float data -> _c format_float oc data
     | `string data -> _c format_string oc data
     | `factor data -> _factor oc data)
  
	
let format_frame oc l = 
  output_string oc "data.frame(" ;
  format_sequence format_frame_item oc (List.enum l) ;
  output_string oc ")"

let frame = bind format_frame

let format_matrix oc m = 
  output_string oc "matrix(c(" ;
  format_sequence 
    format_float oc 
    (Enum.concat (Enum.map Array.enum (Array.enum m))) ;
  fprintf oc "),nrow=%d,ncol=%d)" (Array.length m) (Array.length m.(0))

let matrix = bind format_matrix

let pdf oc path = 
  c oc "pdf('%s')" path

let devoff oc = c oc "dev.off()"

let close oc = 
  flush oc ;
  ignore (Unix.close_process_out oc)

