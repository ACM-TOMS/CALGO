# Makefile for driver programs
# Makefile for GNU gcc compiler


# Talbot Suite folder
TalbotDIR = ../Src


# ----------------compile and link---------------------------------
SEQ:
# without -fopenmp and -lgomp options (as sequential)
	gcc -c -std=gnu99 -Wall -pedantic ./OMP_main.c $(TalbotDIR)/COM_Talbot_pack.c $(TalbotDIR)/OMP_Talbot_pack.c
	gcc -o SEQ_talbot.exe *.o -lm -lrt

OMP:
# with -fopenmp and -lgomp options
	gcc -c -std=gnu99 -fopenmp -Wall -pedantic ./OMP_main.c $(TalbotDIR)/COM_Talbot_pack.c $(TalbotDIR)/OMP_Talbot_pack.c
	gcc -o OMP_talbot.exe *.o -lm -lgomp -lrt

MPI:
	mpicc -c -std=gnu99 -Wall -pedantic ./MPI_main.c $(TalbotDIR)/COM_Talbot_pack.c $(TalbotDIR)/MPI_Talbot_pack.c
	mpicc -o MPI_talbot.exe *.o -lm -lrt

HYB:
	mpicc -c -std=gnu99 -fopenmp -Wall -pedantic ./HYB_main.c $(TalbotDIR)/COM_Talbot_pack.c $(TalbotDIR)/HYB_Talbot_pack.c
	mpicc -o HYB_talbot.exe *.o -lm -lgomp -lrt



# ----------------"clean" rule-------------------------------------
clean: 
	/bin/rm -f  ./*.exe
	/bin/rm -f  ./*.o
