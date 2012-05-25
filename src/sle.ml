open Batteries

let sp = Printf.sprintf

let shell s =
  if Sys.command s != 0 
  then failwith (sp "shell call failed:\n  %s\n" s)
    
let sh fmt = 
  Printf.ksprintf shell fmt

let load fn = 
  let ic = open_in fn in
  let v  = Marshal.from_channel ic in
  close_in ic ; v

let save fn v = 
  let oc = open_out fn in
  Marshal.to_channel oc v [] ;
  close_out oc

let with_tmp_filename f =
  with_dispose ~dispose:Sys.remove f (Filename.temp_file "rar" "")


let count_occurences ch s =
  let accu = ref 0 in
  for i = 0 to String.length s - 1 do
    if s.[i] = ch then incr accu
  done ;
  !accu 
    
let string_split sep x = 
  let n = String.length x in 
  let m = count_occurences sep x + 1 in 
  let res = Array.make m "" in
  let rec search k i j =
    if j >= n then res.(k) <- String.sub x i (j - i)
    else (
      if x.[j] = sep then (
	res.(k) <- String.sub x i (j - i) ;
	search (k + 1) (j + 1) (j + 1) 
      )
      else search k i (j + 1)
    )
  in
  search 0 0 0 ;
  res

let string_split_noeps sep x = 
  let n = String.length x in 
  let rec search i j accu =
    if j >= n then ((String.sub x i (j - i)) :: accu)
    else (
      if x.[j] = sep then
	let accu' = 
	  if i = j then accu
	  else (String.sub x i (j - i)) :: accu
	in
	search (j + 1) (j + 1) accu'
      else
	search i (j + 1) accu
    )
  in Array.of_list (List.rev (search 0 0 []))

let string_unquote s = 
  if 
    String.length s >= 2 && 
    s.[0] = '"' && 
    s.[String.length s - 1] = '"'
  then String.sub s 1 (String.length s - 2)
  else s
    
let string_quote s = "\"" ^ s ^ "\""



















