# Makefile for driver programs
# Makefile for GNU gcc compiler


TALBOTdir = ../../TalbotSuiteDE
EXdir = ..

# ----------------compile and link---------------------------------
# 	with -fopenmp and -lgomp options
OMP:
	gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wno-unused-variable -Wno-comment ./OMP_main.c \
	                                   ./OMP_LTsamples_fun.c ./SEQ_LTsamples_fun.c ./LTfun2.c \
	                                   $(EXdir)/LTfun.c $(EXdir)/ILTfun.c $(EXdir)/LTsings.c  \
	              $(TALBOTdir)/COM/COM_Talbot_pack.c $(TALBOTdir)/COM_DE/COM_Talbot_pack_DE.c \
	        $(TALBOTdir)/FUN_DE/SEQ_Talbot_pack_DE.c $(TALBOTdir)/FUN_DE/OMP_Talbot_pack_DE.c

	gcc -o ex0.exe *.o -lm -lrt -lgomp

#	echo "\n\tRUN AS: ./ex0.exe test jFUN  where  test=1,2,3  and  jFUN=0,1,...,4\n"



# ----------------"clean" rule-------------------------------------
clean: 
	/bin/rm -f  ./*.exe
	/bin/rm -f  ./*.o
