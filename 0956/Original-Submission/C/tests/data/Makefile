################################################################################
# Targets
################################################################################
#.PHONY: ../${LIBNAME}
#	cd ..; make clean; make

GV_FILES := $(wildcard *.gv)
SVG_FILES := $(notdir $(GV_FILES:.gv=.svg))

#.SUFFIXES:
#.SUFFIXES: .gv .svg

%.svg: %.gv
	dot -Tsvg $< -o $*.svg

svg: $(SVG_FILES)

clean:
	rm -f $(GV_FILES) $(SVG_FILES)
