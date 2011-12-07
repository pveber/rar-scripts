open Batteries
open Oregon
open Genome
open Oregon.Target
open Oregon.Target.Infix
open Sle.Infix
open Printf

type sample = B.Rnaseq.sample
type condition = string * sample list * sample list

type annot = {
  id : string ;
  baseMeanA : float array ; (* for all conditions *)
  baseMeanB : float array ;
  padj      : float array ;
}

let hfun_of_deseq x = 
  Tsv.to_fun
    ~from:(fun g -> g#id)
    ~into:(fun g -> g)
    x

let make design genes = 
  let deseqs = List.map (fun (_,x,y) -> B.Rnaseq.DESeq.run x y) design in V.make
  (object
     method id = "Rnaseq_annotation.make[r1]"
     method deps = [] ++* deseqs ++ genes
     method build = 
       let deseqs = List.map hfun_of_deseq deseqs |> Array.of_list in
       let genes = genes#value in
       let design = Array.of_list design in
       design,
       Array.map
	 (fun g -> 
	    let data = Array.map (fun deseq -> g.Gene.id --> deseq) deseqs in 
	    { id = g.Gene.id ;
	      baseMeanA = Array.map (fun x -> x#baseMeanA) data ;
	      baseMeanB = Array.map (fun x -> x#baseMeanB) data ;
	      padj = Array.map (fun x -> x#padj) data ; })
	 genes
   end)

let label l (c,_,_) = 
  sprintf "%s %s" c l

let string_of_float f = 
  if classify_float f = FP_nan then "NA"
  else string_of_float f

let tsv fn target = 
  let conditions, annots = target#value in
  let header = Array.concat [
    [| "gene id" |] ;
    Array.map (label "baseMeanA") conditions ;
    Array.map (label "baseMeanB") conditions ;
    Array.map (label "padj") conditions ;
  ]
  in
  Array.enum annots 
  /@ (fun a -> Array.concat [
	[| a.id |] ;
	Array.map string_of_float a.baseMeanA ;
	Array.map string_of_float a.baseMeanB ;
	Array.map string_of_float a.padj
      ])
  |> Enum.append (Enum.singleton header)
  |> Enum.map Tsv.string_of_line
  |> File.write_lines fn

      
