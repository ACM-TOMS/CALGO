
INC	:= -I$(CUDA_HOME)/include -I.
LIB	:= -L$(CUDA_HOME)/lib64 -lcudart -lm

NVCCFLAGS	:= -O3 -arch=sm_30 --ptxas-options=-v --use_fast_math

all:	poissinv_timing poissinv_check

poissinv_timing:	poissinv_timing.cu poissinv_cuda.h Makefile_sm30
	nvcc poissinv_timing.cu -o poissinv_timing $(INC) $(NVCCFLAGS) $(LIB)

poissinv_check:		poissinv_check.cu poissinv_quad.o poissinv_cuda.h Makefile_sm30
	nvcc poissinv_check.cu poissinv_quad.o -o poissinv_check $(INC) $(NVCCFLAGS) $(LIB) -lquadmath

poissinv_quad.o:	poissinv_quad.c
	gcc -c poissinv_quad.c

clean:
	rm -f poissinv_timing poissinv_check
