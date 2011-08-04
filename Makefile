OCB=ocamlbuild -classic-display -cflags -I,+oregon -lflags -I,+oregon,sle.cmo,genome.cmo,oregon.cmo,labs.cmo,gzmConfig.cmo

top:
	$(OCB) rar.top

%.byte:
	$(OCB) $@

clean:
	$(OCB) -clean
