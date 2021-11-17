# Makefile for driver programs
# Makefile for GNU gcc compiler


# Talbot Suite folder
TALBOTdir = ../../../TalbotSuiteDE

# Example folder
EXdir = ../..


# ----------------compile and link---------------------------------
# 	with -fopenmp and -lgomp options (OpenMP parallel)


acc:
	gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wno-unused-but-set-variable -Wno-comment ./OMP_main_ACCURACY.c ./OMP_LTsamples_fun.c \
	                                                    ./internalMeshPts.c $(EXdir)/LTfun2.c $(EXdir)/ILTfun2.c $(EXdir)/LTsings2.c \
	            $(TALBOTdir)/COM/COM_Talbot_pack.c $(TALBOTdir)/COM_DE/COM_Talbot_pack_DE.c $(TALBOTdir)/FUN_DE/OMP_Talbot_pack_DE.c
#   gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wextra ./OMP_main_ACCURACY.c ./OMP_LTsamples_fun.c ./internalMeshPts.c $(EXdir)/LTfun2.c $(EXdir)/ILTfun2.c $(EXdir)/LTsings2.c $(TALBOTdir)/COM/COM_Talbot_pack.c $(TALBOTdir)/COM_DE/COM_Talbot_pack_DE.c $(TALBOTdir)/FUN_DE/OMP_Talbot_pack_DE.c
	gcc -o ex5_acc.exe ./OMP_main_ACCURACY.o ./OMP_LTsamples_fun.o ./internalMeshPts.o ./LTfun2.o ./ILTfun2.o ./LTsings2.o ./COM_Talbot_pack.o ./COM_Talbot_pack_DE.o ./OMP_Talbot_pack_DE.o -lm -lgomp


time:
	gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wno-unused-but-set-variable -Wno-comment ./OMP_main_TIMES.c ./OMP_LTsamples_fun.c \
	                                                                    ./internalMeshPts.c $(EXdir)/LTfun2.c $(EXdir)/LTsings2.c \
	         $(TALBOTdir)/COM/COM_Talbot_pack.c $(TALBOTdir)/COM_DE/COM_Talbot_pack_DE.c $(TALBOTdir)/FUN_DE/OMP_Talbot_pack_DE.c
#	gcc -c -std=gnu99 -fopenmp -Wall -pedantic -Wextra ./OMP_main_TIMES.c ./OMP_LTsamples_fun.c ./internalMeshPts.c $(EXdir)/LTfun2.c $(EXdir)/LTsings2.c $(TALBOTdir)/COM/COM_Talbot_pack.c $(TALBOTdir)/COM_DE/COM_Talbot_pack_DE.c $(TALBOTdir)/FUN_DE/OMP_Talbot_pack_DE.c
	gcc -o ex5_time.exe ./OMP_main_TIMES.o ./OMP_LTsamples_fun.o ./internalMeshPts.o ./LTfun2.o ./LTsings2.o ./COM_Talbot_pack.o ./COM_Talbot_pack_DE.o ./OMP_Talbot_pack_DE.o -lm -lrt -lgomp



# ----------------"clean" rule-------------------------------------
clean: 
	/bin/rm -f  ./*.exe
	/bin/rm -f  ./*.o
