type entry = {
  ac : string ;
  id : string ;
  name : string ;
  mat : float array array
}

val of_file : string -> entry list

val motif_library_item_of_entry : entry -> Oregon.Motif_library.item 

val motif_library_of_file : string -> Oregon.Motif_library.t
