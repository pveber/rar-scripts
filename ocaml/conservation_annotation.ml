open Batteries
open Genome
open Oregon
open Printf
open Scanf

let digit = function
'0'..'9' -> true
  | _ -> false

let hgWiggle_samechr u v = 
  not (String.length u > 0 && digit u.[0] && String.length v > 0 && not (digit v.[0]))

let hgWiggle_output e = 
  (e // (fun s -> String.length s = 0 || s.[0] <> '#')
      |> Enum.group_by hgWiggle_samechr)
  /@ (fun e -> 
    let first = Option.get (Enum.get e) in
    let chr = sscanf first "variableStep chrom=%s" Std.identity in
    chr, Enum.map ((fun x -> String.split x "\t") |- Tuple2.mapn int_of_string float_of_string) e |> List.of_enum)
  |> Enum.fold (fun accu (chr,posz) -> PMap.add chr posz accu) PMap.empty
      


let conservation_scores locz = 
  Sle.with_tmp (fun fn -> 
    locz /@ (Bed.unparser |- Tsv.string_of_line) |> (File.write_lines fn) ;
    Sle.shout "cd data/conservation/mm9 && hgWiggle -bedFile=%s -db=mm9 phastCons30wayPlacental" fn
    |> Array.enum
    |> hgWiggle_output)
    
let make locz fn = 
  let data = conservation_scores (Enum.clone locz) in
  Enum.map 
    (fun loc -> Location.(
      (try PMap.find loc.chr data |> List.enum with Not_found -> Enum.empty ()) 
      //@ (fun (i, x) -> if i >= loc.st && i <= loc.ed then Some x else None)
			     |> Enum.fold ( +. ) 0. ))
    (Enum.clone locz)
  |> (fun x -> Enum.combine (locz, x))
  |> Enum.map (fun (loc,x) -> Array.append (Bed.unparser loc) [| string_of_float x |])
  |> Enum.map Tsv.string_of_line
  |> File.write_lines fn
