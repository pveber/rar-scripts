open Batteries
open ParserCo

type entry = {
  ac : string ;
  id : string ;
  name : string ;
  mat : float array array
}

let source path = 
  Source.of_enum (File.lines_of path) 0 (fun _ n -> n + 1)

let ( >>+ ) p f =
  p >>= (fun x -> try f x with _ -> fail)

let parse_line l = 
  let n = String.length l in
  if n > 3 then
    String.sub l 0 2, String.sub l 4 (n - 4)
  else
    raise (Invalid_argument ("Transfac.parse_line: " ^ l))

let line = 
  any >>+ (fun l -> return (parse_line l)) 

let xx = exactly "XX"

let ss = exactly "//"

let doesnt_start_with prefix = 
  satisfy (fun l -> not (String.starts_with l prefix))

let mapnext prefix f =
  ignore_zero_plus (doesnt_start_with prefix) >>= fun () ->
  line >>= fun (_,x) ->
  return (f x)

let skip_past prefix = 
  ignore_zero_plus (doesnt_start_with prefix) >>> any >>= fun _ -> return ()

let space = Str.regexp "[ \t]+"
let split l = Array.of_list (Str.split space l)

let count_matrix_line (_,l) = 
  Array.sub (split l) 0 4 |> Array.map float_of_string

let count_matrix = 
  one_plus line >>= fun l ->
  xx >>= fun _ ->
  return (List.enum l |> Enum.map count_matrix_line |> Array.of_enum)

let entry = 
  mapnext "AC" identity >>= fun ac ->
  mapnext "ID" identity >>= fun id ->
  mapnext "NA" identity >>= fun name ->
  skip_past "P0" >>= fun () ->
  count_matrix >>= fun mat ->
  ignore_zero_plus (none_of [ "//" ]) >>= fun () ->
  ss >>= fun _ ->
  return { ac ; id ; name ; mat }
  
let file = 
  times 3 any >>>
  one_plus entry >>= fun l ->
  eof >>= fun () -> 
  return l

let of_file path = match run file (source path) with
    Ok x -> x
  | Bad _ -> assert false
