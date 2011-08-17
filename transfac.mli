type entry = {
  ac : string ;
  id : string ;
  name : string ;
  mat : float array array
}

val of_file : string -> entry list
