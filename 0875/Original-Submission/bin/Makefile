
default:
	make all

DIFF = diff
RM = rm -rf
D = ./

truss1test:
	-@${D}dsdp5 truss1.dat-s -dlogsummary 0 -print 1 > truss1.tmp 2>&1; \
	if ( ${DIFF} truss1.tmp output.truss1 ) then echo "  No error: truss1 "; \
	else echo "  Possible Error (truss1): Differences between this output and reference output:"; fi; \
	${RM} truss1.tmp

mcp100test:
	-@${D}dsdp5 mcp100.dat-s -dlogsummary 0 -print 1 > mcp100.tmp 2>&1; \
	if ( ${DIFF} mcp100.tmp output.mcp100 ) then echo "  No error: mcp100 "; \
	else echo "  Possible Error (mcp100): Differences between this output and reference output:"; fi; \
	${RM} mcp100.tmp

maxcuttest:
	-@${D}maxcut graph1 > maxcut.tmp 2>&1; \
	if ( ${DIFF} maxcut.tmp output.maxcut ) then echo "  No error: Maxcut "; \
	else echo "  Possible Error (maxcut): Differences between this output and reference output:"; fi; \
	${RM} maxcut.tmp

thetatest:
	-@${D}theta graph1 > theta.tmp 2>&1; \
	if ( ${DIFF} theta.tmp output.theta ) then echo "  No Error: Theta "; \
	else echo "  Possible Error (theta): Differences between this output and reference output:"; fi; \
	${RM} theta.tmp

output:
	${RM} output.truss1 output.mcp100 output.maxcut output.theta;
	${D}dsdp5 truss1.dat-s -dlogsummary 0 -print 1 > output.truss1;
	${D}dsdp5 mcp100.dat-s -dlogsummary 0 -print 1 > output.mcp100;
	${D}maxcut graph1 > output.maxcut;
	${D}theta graph1 > output.theta;

test:
	-make truss1test
	-make mcp100test
	-make maxcuttest
	-make thetatest

all:
	make test

clean:
	${RM} theta maxcut color stable dsdp5

