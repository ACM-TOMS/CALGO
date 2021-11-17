# Makefile for driver programs
# Makefile for GNU gcc compiler


# Talbot Suite folder
TALBOTdir = ../../../TalbotSuiteDE
# Example folder
EXdir = ../..


# ----------------compile and link---------------------------------
# 	without -fopenmp and -lgomp options (sequential)

acc:
	gcc -c -std=gnu99 -Wall -pedantic -Wno-unused-parameter ./SEQ_main_ACCURACY.c ./SEQ_LTsamples_ode.c ./ode.c $(EXdir)/ILTfun2.c $(EXdir)/LTsings2.c \
	                $(TALBOTdir)/COM/COM_Talbot_pack.c $(TALBOTdir)/COM_DE/COM_Talbot_pack_DE.c $(TALBOTdir)/FUN_DE/SEQ_Talbot_pack_DE.c
	gcc -o ex4_acc.exe ./SEQ_main_ACCURACY.o ./SEQ_LTsamples_ode.o ./ode.o ./ILTfun2.o ./LTsings2.o ./COM_Talbot_pack.o ./COM_Talbot_pack_DE.o ./SEQ_Talbot_pack_DE.o -lm -lrt

time:
	gcc -c -std=gnu99 -Wall -pedantic -Wno-unused-parameter ./SEQ_main_TIMES.c ./SEQ_LTsamples_ode.c ./ode.c $(EXdir)/LTsings2.c \
	                $(TALBOTdir)/COM/COM_Talbot_pack.c $(TALBOTdir)/COM_DE/COM_Talbot_pack_DE.c $(TALBOTdir)/FUN_DE/SEQ_Talbot_pack_DE.c
	gcc -o ex4_time.exe ./SEQ_main_TIMES.o ./SEQ_LTsamples_ode.o ./ode.o ./LTsings2.o ./COM_Talbot_pack.o ./COM_Talbot_pack_DE.o ./SEQ_Talbot_pack_DE.o -lm -lrt


# ----------------"clean" rule-------------------------------------
clean: 
	/bin/rm -f  ./*.exe
	/bin/rm -f  ./*.o
