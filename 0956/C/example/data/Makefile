GV_FILES := $(wildcard *.gv)
SVG_FILES := $(notdir $(GV_FILES:.gv=.svg))
.PHONY: trees clean

%.svg : %.gv
	dot -Tsvg $< -o $*.svg

trees: $(SVG_FILES)

clean:
	${RM} -f ${GV_FILES} ${SVG_FILES}
