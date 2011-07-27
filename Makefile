OCB=ocamlbuild -classic-display

top:
	$(OCB) rar.top

%.byte:
	$(OCB) $@

clean:
	$(OCB) -clean
