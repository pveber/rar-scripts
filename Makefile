OCB=ocamlbuild -classic-display

%.byte:
	$(OCB) $@
