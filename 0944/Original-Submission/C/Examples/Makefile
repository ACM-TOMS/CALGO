# Makefile for simple main examples
# Makefile for GNU gcc compiler


# Talbot Suite folder
TalbotDIR = ../Src

# --------------- ..._Talbot1 function ----------------------------

SEQ1:
# without -fopenmp and -lgomp options
	gcc -c -std=c99 -Wall -pedantic ./OMP_main1.c $(TalbotDIR)/COM_Talbot_pack.c $(TalbotDIR)/OMP_Talbot_pack.c
	gcc -o SEQ_talbot.exe *.o -lm

OMP1:
# with -fopenmp and -lgomp options
	gcc -c -std=c99 -fopenmp -Wall -pedantic ./OMP_main1.c $(TalbotDIR)/COM_Talbot_pack.c $(TalbotDIR)/OMP_Talbot_pack.c
	gcc -o OMP_talbot.exe *.o -lm -lgomp

MPI1:
	mpicc -c -std=c99 -Wall -pedantic ./MPI_main1.c $(TalbotDIR)/COM_Talbot_pack.c $(TalbotDIR)/MPI_Talbot_pack.c
	mpicc -o MPI_talbot.exe *.o -lm


# --------------- ..._Talbot2 function ----------------------------

SEQ2:
# without -fopenmp and -lgomp options
	gcc -c -std=c99 -Wall -pedantic ./OMP_main2.c $(TalbotDIR)/COM_Talbot_pack.c $(TalbotDIR)/OMP_Talbot_pack.c
	gcc -o SEQ_talbot.exe *.o -lm

OMP2:
# with -fopenmp and -lgomp options
	gcc -c -std=c99 -fopenmp -Wall -pedantic ./OMP_main2.c $(TalbotDIR)/COM_Talbot_pack.c $(TalbotDIR)/OMP_Talbot_pack.c
	gcc -o OMP_talbot.exe *.o -lm -lgomp

MPI2:
	mpicc -c -std=c99 -Wall -pedantic ./MPI_main2.c $(TalbotDIR)/COM_Talbot_pack.c $(TalbotDIR)/MPI_Talbot_pack.c
	mpicc -o MPI_talbot.exe *.o -lm


# --------------- HYB_Talbot3 function ----------------------------
HYB3:
	mpicc -c -std=c99 -fopenmp -Wall -pedantic ./HYB_main.c $(TalbotDIR)/COM_Talbot_pack.c $(TalbotDIR)/HYB_Talbot_pack.c
	mpicc -o HYB_talbot.exe *.o -lm -lgomp



# ----------------"clean" rule-------------------------------------
clean: 
	/bin/rm -f  ./*.exe
	/bin/rm -f  ./*.o
