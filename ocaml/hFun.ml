type ('a, 'b) t = 
    Open of ('a -> 'b) * ('a,'b) Hashtbl.t
  | Closed of ('a,'b) Hashtbl.t

let ht = function 
  | Open (_,ht) -> ht
  | Closed ht -> ht

module Infix = struct
  let ( <-- ) hfun (x, v) = Hashtbl.replace (ht hfun) x v
    
  let ( $ ) hfun x = match hfun with
    | Open (f,t) -> 
      begin
        try Hashtbl.find t x 
	with Not_found -> (
	  let y = f x in 
	  Hashtbl.replace t x y ;
	  y
	)
      end
    | Closed t -> Hashtbl.find t x
end

open Infix

let make ?(init = BatEnum.empty ()) f = 
  let r  = Open (f, Hashtbl.create 251) in
  BatEnum.iter 
    (fun x -> r <-- (x, f x))
    init ;
  r

let close = function
  | Open (f,t) -> Closed t
  | x -> x

let domain hfun = 
  BatHashtbl.keys (ht hfun)

let image hfun = 
  BatHashtbl.values (ht hfun)

let cache f = 
  let hfun = make f in
  fun x -> hfun $ x

(* let of_enum f support =  *)
(*   let r = Hashtbl.create 1024 in *)
(*   Enum.iter (fun x -> r <-- (x, f x)) support ; *)
(*   r *)

(* let make2 domain image =  *)
(*   let n = List.length domain  *)
(*   and m = List.length image in *)
(*   let _ = if n <> m then raise (Invalid_argument "hfun_extension") in *)
(*   let r = Hashtbl.create n in *)
(*   List.iter2 (fun x y -> r <-- (x, y)) domain image ; *)
(*   r *)

(* let preimage f support =  *)
(*   let r = Hashtbl.create 251 in *)
(*   Enum.iter (fun x -> r <-- (f x, x)) support ; *)
(*   r *)

(* let compose f hfun = Hashtbl.map (fun _ v -> f v) hfun *)

(* let image hfun = Hashtbl.values hfun *)
