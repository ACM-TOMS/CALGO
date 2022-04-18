/*
 -------------------------------------------------------------------------
 File taputil3.h of ADOL-C version 1.6 as of January 1,   1995          
 Included in ---->
                  fos_reverse.c
                  fov_reverse.c 
                  hos_forward.c
                  hos_reverse.c
                  hov_reverse.c
		  taputil2.c
                  utils.c


 ------------------------------------------------------------------------- 

 taputil3.h contains the prototypes for the functions 
 defined in taputil3.c.  The functions provide
 initialization and stopage of the taping process, as well as
 statistics gathering functions.

*/

/* Stats functions */

extern void get_op_stats(int,char **,int *,int *,unsigned char **);
extern void get_int_stats(int,char **,int *,int *,locint **);
extern void get_val_stats(int,char **,int *,int *,double **);
extern void tapestats(short,int *);

/* Tracing */

extern void start_trace(short,int);
extern void stop_trace(int,int);






