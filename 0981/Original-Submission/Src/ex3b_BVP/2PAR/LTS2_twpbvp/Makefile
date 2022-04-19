# Makefile for driver programs
# Makefile for GNU gcc compiler


# Talbot Suite folder
TALBOTdir = ../../../TalbotSuiteDE

# Example folder
EXdir = ../..


# ----------------compile and link---------------------------------
# 	with -fopenmp and -lgomp options (OpenMP parallel)
#	The extension of FORTRAN files has been changed into .F to enable the preprocessing directives


acc:
	gfortran -c -fopenmp twpbvp.F prtfun.F testfun.F
	gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wno-unused-but-set-variable -Wno-comment ./OMP_main_ACCURACY.c ./OMP_LTsamples_twpbvp.c \
	                                                                                             $(EXdir)/ILTfun2.c $(EXdir)/LTsings2.c \
				   $(TALBOTdir)/COM/COM_Talbot_pack.c $(TALBOTdir)/COM_DE/COM_Talbot_pack_DE.c $(TALBOTdir)/FUN_DE/OMP_Talbot_pack_DE.c
#   gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wextra ./OMP_main_ACCURACY.c ./OMP_LTsamples_twpbvp.c $(EXdir)/ILTfun2.c $(EXdir)/LTsings2.c $(TALBOTdir)/COM/COM_Talbot_pack.c $(TALBOTdir)/COM_DE/COM_Talbot_pack_DE.c $(TALBOTdir)/FUN_DE/OMP_Talbot_pack_DE.c
	gcc -o ex3_acc.exe ./OMP_main_ACCURACY.o ./OMP_LTsamples_twpbvp.o ./ILTfun2.o ./LTsings2.o ./COM_Talbot_pack.o ./COM_Talbot_pack_DE.o ./OMP_Talbot_pack_DE.o twpbvp.o prtfun.o testfun.o -lm -lgfortran -lgomp


time:
	gfortran -c -fopenmp twpbvp.F prtfun.F testfun.F
	gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wno-unused-but-set-variable -Wno-comment ./OMP_main_TIMES.c ./OMP_LTsamples_twpbvp.c $(EXdir)/LTsings2.c \
									$(TALBOTdir)/COM/COM_Talbot_pack.c $(TALBOTdir)/COM_DE/COM_Talbot_pack_DE.c $(TALBOTdir)/FUN_DE/OMP_Talbot_pack_DE.c
#	gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wextra ./OMP_main_TIMES.c ./OMP_LTsamples_twpbvp.c $(EXdir)/LTsings2.c $(TALBOTdir)/COM/COM_Talbot_pack.c $(TALBOTdir)/COM_DE/COM_Talbot_pack_DE.c $(TALBOTdir)/FUN_DE/OMP_Talbot_pack_DE.c
	gcc -o ex3_time.exe ./OMP_main_TIMES.o ./OMP_LTsamples_twpbvp.o ./LTsings2.o ./COM_Talbot_pack.o ./COM_Talbot_pack_DE.o ./OMP_Talbot_pack_DE.o twpbvp.o prtfun.o testfun.o -lm -lgfortran -lrt -lgomp



# ----------------"clean" rule-------------------------------------
clean: 
	/bin/rm -f  ./*.exe
	/bin/rm -f  ./*.o
