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

type row = string array

let is_comment l = l = "" || l.[0] = '#' 

let row_of_string = string_split '\t'
let string_of_row r = 
  String.concat (String.make 1 '\t') (Array.to_list r)

let unparse_row oc = function
    [||] -> output_char oc '\n'
  | labels -> 
      output_string oc labels.(0) ;
      for i = 1 to Array.length labels - 1 do
	output_char oc '\t' ; 
	output_string oc labels.(i) ; 
      done ;
      output_char oc '\n' 

let enum fn =
  File.lines_of fn
  // (fun x -> not (is_comment x))
  /@ row_of_string
