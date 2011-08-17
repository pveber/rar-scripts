open Sle

let less ?(p = []) x = sh "less %s" (String.concat "/" (x#path :: p))

let evince pdf = sh "evince %s" pdf#path
