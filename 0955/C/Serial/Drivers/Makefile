CC	:= gcc
# CC	:= icc

all:	poissinv_timing poissinv_check

poissinv_timing:	poissinv_timing.c poissinv.h Makefile
	${CC} -O3 poissinv_timing.c -o poissinv_timing -I. -lm

poissinv_check:		poissinv_check.c poissinv_quad.o poissinv.h Makefile
	${CC} -O3 poissinv_check.c poissinv_quad.o -o poissinv_check -I. -lm -lquadmath

poissinv_quad.o:	poissinv_quad.c
	${CC} -O3 -c poissinv_quad.c

clean:
	rm -f poissinv_timing poissinv_check

