open Batteries
open Printf
open Genome
open Oregon
open Target.Infix

let fetch ~org_from ~org_to ~path = 
  Sle.sh "cd %s && wget ftp://hgdownload.cse.ucsc.edu/goldenPath/%s/liftOver/%sTo%s.over.chain.gz" path org_from org_from (String.capitalize org_to)


let import_aux_to_file locz =
  let fn = Filename.temp_file "ore" ".locations" in
  locz /@ Location.to_string
  |> File.write_lines fn ;
  fn
    
let chain_file ~org_from ~org_to =
  let path = sprintf "resources/liftOver/%s" org_from in
  let file = sprintf "%s/%sTo%s.over.chain" path org_from (String.capitalize org_to) in
  if not (Sys.file_exists file) then (
    Sle.sh "mkdir -p %s" path ;
    fetch ~org_from ~org_to ~path ;
    Sle.sh "cd %s && gunzip %s.gz" path (Filename.basename file)
  ) ;
  file
    
let conversion' ~org_from ~org_to locz = 
  let oldf = import_aux_to_file locz and
      newf = Filename.temp_file "ore" "" and
      unmp = Filename.temp_file "ore" "" in
  let resz = 
    Sle.sh "liftOver -positions %s %s %s %s" oldf (chain_file ~org_from ~org_to) newf unmp ;
    File.lines_of newf
    /@ Location.of_string
    |> List.of_enum
  in (* crade, mais il faut que les fichiers soient lus avant d'être effacés *)
  Sle.sh "rm -f %s %s %s liftOver_*" oldf newf unmp ;
  List.enum resz

let conversion ~org_from ~org_to bed = Target.F.make 
  (object
    method id = sprintf "Ucsc_lift_over.conversion[%s,%s]" org_from org_to
    method deps = [] ++ bed
    method build path =
      (Tsv.enum bed /@ (fun x -> x#loc) |> conversion' ~org_from ~org_to)
      /@ Bed.unparser
      /@ Tsv.string_of_line
      |> File.write_lines path
    method ty = Bed.basic_ty false
   end)
