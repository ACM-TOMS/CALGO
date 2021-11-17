# Makefile for driver programs
# Makefile for GNU gcc compiler


# Talbot Suite folder
TalbotDIR = ../../TalbotSuiteDE
EXdir = ..


# ----------------compile and link---------------------------------
# 	without -fopenmp and -lgomp options (sequential)
SEQ:
	gcc -c -std=gnu99 -Wall -pedantic -Wno-unused-variable ./SEQ_main.c $(EXdir)/ILTfun.c ./SEQ_LTsamples_fun.c $(EXdir)/LTfun.c $(EXdir)/LTsings.c \
	          $(TalbotDIR)/COM/COM_Talbot_pack.c $(TalbotDIR)/COM_DE/COM_Talbot_pack_DE.c \
	                                             $(TalbotDIR)/FUN_DE/SEQ_Talbot_pack_DE.c \
	                                                                  ./SEQ_Talbot_pack.c   # with no correction to NOPTS
#	                                                   $(TalbotDIR)/FUN/SEQ_Talbot_pack.c   # with correction to NOPTS
#	gcc -c -std=gnu99 -Wall -pedantic -Wextra ...
	gcc -o ex0.exe *.o -lm -lrt


# ----------------"clean" rule-------------------------------------
clean: 
	/bin/rm -f  ./*.exe
	/bin/rm -f  ./*.o
