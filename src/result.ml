open Batteries
open Sle

type 'a t = {
  path : string ;
  build : string -> unit ;
  read : string -> 'a ;
}
type 'a serializer = (string -> 'a -> unit) * (string -> 'a)
    
let make ~path ~build ~read = 
  { path ; build ; read }
  
let marshal_serializer = save, load

let value ?(serializer = marshal_serializer) path f = 
  make 
    ~path 
    ~build:(fun path -> (fst serializer) path (f ()))
    ~read:(snd serializer) 

let file path f =
  make 
    ~path
    ~build:f
    ~read:identity

let exists r = Sys.file_exists r.path

let build r = 
  if not (exists r) then (
    sh "mkdir -p %s" (Filename.dirname r.path) ;
    let tmp = Filename.temp_file "" "" in 
    let () = r.build tmp in
    sh "mv %s %s" tmp r.path
  )

let get r =
  build r ; 
  r.read r.path





















