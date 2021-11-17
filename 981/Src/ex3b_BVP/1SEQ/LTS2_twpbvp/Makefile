# Makefile for driver programs
# Makefile for GNU gcc compiler


# Talbot Suite folder
TALBOTdir = ../../../TalbotSuiteDE
# Example folder
EXdir = ../..


# ----------------compile and link---------------------------------
# 	without -fopenmp and -lgomp options (sequential)


acc:
	gfortran -c twpbvp.f prtfun.f testfun.f
	gcc -c -std=gnu99 -Wall -pedantic -Wextra ./SEQ_main_ACCURACY.c ./SEQ_LTsamples_twpbvp.c $(EXdir)/ILTfun2.c $(EXdir)/LTsings2.c \
			   $(TALBOTdir)/COM/COM_Talbot_pack.c $(TALBOTdir)/COM_DE/COM_Talbot_pack_DE.c $(TALBOTdir)/FUN_DE/SEQ_Talbot_pack_DE.c
	gcc -o ex3_acc.exe ./SEQ_main_ACCURACY.o ./SEQ_LTsamples_twpbvp.o ./twpbvp.o ./prtfun.o ./testfun.o ./ILTfun2.o ./LTsings2.o \
												./COM_Talbot_pack.o ./COM_Talbot_pack_DE.o ./SEQ_Talbot_pack_DE.o -lm -lgfortran


time:
	gfortran -c twpbvp.f prtfun.f testfun.f
	gcc -c -std=gnu99 -Wall -pedantic -Wextra ./SEQ_main_TIMES.c ./SEQ_LTsamples_twpbvp.c $(EXdir)/LTsings2.c \
		$(TALBOTdir)/COM/COM_Talbot_pack.c $(TALBOTdir)/COM_DE/COM_Talbot_pack_DE.c $(TALBOTdir)/FUN_DE/SEQ_Talbot_pack_DE.c
	gcc -o ex3_time.exe ./SEQ_main_TIMES.o ./SEQ_LTsamples_twpbvp.o ./twpbvp.o ./prtfun.o ./testfun.o ./LTsings2.o ./COM_Talbot_pack.o ./COM_Talbot_pack_DE.o ./SEQ_Talbot_pack_DE.o -lm -lrt -lgfortran



# ----------------"clean" rule-------------------------------------
clean: 
	/bin/rm -f  ./*.exe
	/bin/rm -f  ./*.o
