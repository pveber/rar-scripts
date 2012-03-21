OCB=ocamlbuild -classic-display -cflags -I,+oregon -tag annot

all: repeat_proportion.native estimate_matrix_threshold.native motif_prevalence.native \
     hexamer_estimation_local_search.native annotations.native nhr_scan_fasta.native \
     nhr_pattern.native chipseq_track_annotation_main.native

top:
	$(OCB) -lflags -I,+oregon,sle.cmo,genome.cmo,oregon.cmo,labs.cmo,gzmConfig.cmo -tag debug rar.top


%.byte: %.ml
	$(OCB) -lflags -I,+oregon,sle.cmo,genome.cmo,oregon.cmo,labs.cmo,gzmConfig.cmo -tag debug $@

%.native: %.ml
	$(OCB) -lflags -I,+oregon,sle.cmx,genome.cmx,oregon.cmx,labs.cmx,gzmConfig.cmx -tag debug $@

annotations.native: annotations.ml
	$(OCB) -lflags -I,+oregon,sle.cmx,genome.cmx,oregon.cmx,labs.cmx,gzmConfig.cmx -tag debug $@

chipseq_track_annotation_main.native: chipseq_track_annotation_main.ml
	$(OCB) -lflags -I,+oregon,sle.cmx,genome.cmx,oregon.cmx,labs.cmx,gzmConfig.cmx -tag debug $@

nhr_genome_scan.native: nhr_genome_scan.ml
	$(OCB) -lflags -I,+oregon,sle.cmx,genome.cmx,oregon.cmx,labs.cmx,gzmConfig.cmx -tag debug $@

test.byte: test.ml
	$(OCB) -lflags -I,+oregon,sle.cmo,genome.cmo,oregon.cmo,labs.cmo,gzmConfig.cmo -tag debug $@

test.native: test.ml
	$(OCB) -lflags -I,+oregon,sle.cmx,genome.cmx,oregon.cmx,labs.cmx,gzmConfig.cmx -tag debug $@

clean:
	$(OCB) -clean
