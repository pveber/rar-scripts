open Batteries
open Printf
open Biocaml

type 'a t = 'a * Range.t

let make chr st ed =
  if st > ed then (
    printf "Location.make: incorrect coordinates %d %d\n" st ed ;
    assert false
  ) ;
  (chr, Range.make st ed)

let upstream ~up ~down strand (chr, { Range.lo ; hi }) =
  match strand with
    | `Sense ->
	make chr (max 0 (lo - up)) (max 0 (lo + down))
    | `Antisense ->
	make chr (max 0 (hi - down)) (max 0 (hi + up))

let size (_,r) = Range.size r

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
(* open Printf *)
(* open Sle *)

(* type t = { *)
(*   chr : string ; *)
(*   st : int ; *)
(*   ed : int *)
(* } *)
(* let to_string l = sprintf "%s:%d-%d" l.chr l.st l.ed *)
(* let move l st ed = { l with st = st ; ed = ed } *)
(* let relmove l st ed = { l with st = l.st + st ; ed = l.ed + ed } *)
(* let relative l st ed = { l with st = l.st + st ; ed = l.st + ed } *)
(* let zoom l factor =  *)
(*   let size = int_of_float (float (l.ed - l.st + 1) *. factor) in  *)
(*   let center = (l.st + l.ed) / 2  *)
(*   in move l (center - size / 2) (center + size / 2) *)
(* let length l = l.ed - l.st + 1 *)


(* let of_string s =  *)
(*   Scanf.sscanf  *)
(*     s "%s@:%d@-%d"  *)
(*     make *)



(* let included_in s s' = *)
(*   if s.chr <> s'.chr then false  *)
(*   else (s.st >= s'.st && s.ed <= s'.ed) *)

(* let intersection s s' =  *)
(*   if s.chr <> s'.chr then false else ( *)
(*     let p s s' = (s.ed >= s'.st) && (s.ed <= s'.ed) *)
(*     in p s s' || p s' s *)
(*   ) *)

(* let inter s s' =  *)
(*   if not (intersection s s') then raise (Invalid_argument "Ucsc.Location.inter") *)
(*   else (make s.chr (max s.st s'.st) (min s.ed s'.ed)) *)


(* let convex_hull = function *)
(*     [] -> raise (Invalid_argument "Location.convex_hull: empty list") *)
(*   | h :: t ->  *)
(*       List.fold_left  *)
(* 	(fun accu l ->  *)
(* 	   if accu.chr = l.chr then *)
(* 	     make l.chr (min l.st accu.st) (max l.ed accu.ed) *)
(* 	   else *)
(* 	     raise (Invalid_argument "Location.convex_hull: locations on several chromosomes")) *)
(* 	h t *)
  
(* let dist s s' =  *)
(*   if s.chr <> s'.chr then raise (Invalid_argument "Ucsc.Location.dist") *)
(*   else ( *)
(*     if intersection s s' then 0 *)
(*     else min (abs (s'.st - s.ed)) (abs (s.st - s'.ed)) *)
(*   ) *)

(* let position ~from loc =  *)
(*   if loc.chr <> from.chr then raise (Invalid_argument "Ucsc.Location.position") *)
(*   else ( *)
(*     if intersection from loc then 0 *)
(*     else ( *)
(*       let a, b = loc.ed - from.st, loc.st - from.ed in *)
(*       if abs a < abs b then a else b *)
(*     ) *)
(*   ) *)
    
(* let compare (x : t) (y : t) = compare x y *)

(* open Batteries *)

(* let enum_merge x y =  *)
(*   Enum.merge (fun x y -> compare x y < 0) x y *)

    
(* let islands x = *)
(*   let next () =  *)
(*     let rec island loc = match Enum.peek x with  *)
(* 	Some loc' ->  *)
(* 	  if intersection loc loc' *)
(* 	  then (Enum.junk x ; island (make loc.chr loc.st loc'.ed)) *)
(* 	  else loc  *)
(*       | None -> raise Enum.No_more_elements in *)
(*     match Enum.get x with  *)
(* 	Some loc -> island loc *)
(*       | None -> raise Enum.No_more_elements in *)
(*   Enum.from next *)

(* type 'a annotation = 'a  *)
(* constraint 'a = < loc : t ; .. > *)
