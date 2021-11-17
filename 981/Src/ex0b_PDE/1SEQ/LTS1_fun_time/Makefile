# Makefile for driver programs
# Makefile for GNU gcc compiler


# Talbot Suite folder
TALBOTdir = ../../../TalbotSuiteDE
# Example folder
EXdir = ../..


# ----------------compile and link---------------------------------
# 	without -fopenmp and -lgomp options (sequential)

time:
	gcc -c -std=gnu99 -Wall -pedantic -Wno-unused-variable -Wno-unused-but-set-variable ./SEQ_main_TIMES.c ./SEQ_LTsamples_fun.c $(EXdir)/LTfun2.c $(EXdir)/LTfun.c $(EXdir)/LTsings2.c \
								$(TALBOTdir)/COM/COM_Talbot_pack.c $(TALBOTdir)/COM_DE/COM_Talbot_pack_DE.c $(TALBOTdir)/FUN/SEQ_Talbot_pack.c $(TALBOTdir)/FUN_DE/SEQ_Talbot_pack_DE.c
#	gcc -c -std=gnu99 -Wall -Wextra -pedantic ./SEQ_main_TIMES.c ./SEQ_LTsamples_fun.c $(EXdir)/LTfun2.c $(EXdir)/LTfun.c $(EXdir)/LTsings2.c $(TALBOTdir)/COM/COM_Talbot_pack.c $(TALBOTdir)/COM_DE/COM_Talbot_pack_DE.c $(TALBOTdir)/FUN/SEQ_Talbot_pack.c $(TALBOTdir)/FUN_DE/SEQ_Talbot_pack_DE.c
	gcc -o ex0_time.exe ./SEQ_main_TIMES.o ./SEQ_LTsamples_fun.o ./LTfun2.o ./LTfun.o ./LTsings2.o ./COM_Talbot_pack.o ./COM_Talbot_pack_DE.o ./SEQ_Talbot_pack.o ./SEQ_Talbot_pack_DE.o -lm -lrt

	


# ----------------"clean" rule-------------------------------------
clean: 
	/bin/rm -f  ./*.exe
	/bin/rm -f  ./*.o
