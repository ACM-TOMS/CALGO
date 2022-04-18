/*
  --------------------------------------------------------------------
  file tayutil.h of ADOL-C version 1.6 as of January 1,   1995
  Included in -->
                 fos_reverse.c
                 fov_reverse.c
                 hos_forward.c
                 hos_reverse.c
                 hov_reverse.c
                 taputil1.c
                 taputil3.c
                 tayutil.c
                 utils.c
   
   -------------------------------------------------------------------
   tayutil.h defines the prototypes for the functions from 
   tayutil.c.  See tayutil.c for an explanation of the functionality 
   of these routines.

*/

#ifndef __STDC__
   int unlink(char *); 
#endif
int taylor_access();
void close_taylor();
void taylor_begin(int, double**,int);
void taylor_close(int,int,int);
void write_taylor(locint, int);
void get_taylor(locint);
void get_taylors(locint, int);
void write_scaylor(revreal);
void write_scaylors(revreal*,int);
void taylor_back(revreal*, int*, int*, int*);
void taylor_back2(revreal**, int*, int*, int*);


