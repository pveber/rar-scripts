OCB=ocamlbuild -classic-display -cflags -I,+oregon 

top:
	$(OCB) -lflags -I,+oregon,sle.cmo,genome.cmo,oregon.cmo,labs.cmo,gzmConfig.cmo -tag debug rar.top

%.byte: %.ml
	$(OCB) -lflags -I,+oregon,sle.cmo,genome.cmo,oregon.cmo,labs.cmo,gzmConfig.cmo -tag debug $@

%.native: %.ml
	$(OCB) -lflags -I,+oregon,sle.cmx,genome.cmx,oregon.cmx,labs.cmx,gzmConfig.cmx -tag debug $@

clean:
	$(OCB) -clean
