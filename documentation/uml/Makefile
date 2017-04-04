# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki

DIAGRAMS = $(patsubst %.puml,%.pdf,$(wildcard *.puml))

.PHONY: all
all: $(DIAGRAMS)

.SECONDEXPANSION:

%.pdf: %.svg
	inkscape $< --export-pdf=$@

%.svg: %.puml $$(shell sed -n -e 's/!include[ ]*\([^ \n]*\)/\1/p' $$*.puml)
	plantuml -tsvg $<

.PHONY: clean
clean:
	-rm -f *.svg *.pdf
