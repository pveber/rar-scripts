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




















