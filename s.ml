open Sle

let less ?(p = []) x = sh "less %s" (String.concat "/" (x#path :: p))
