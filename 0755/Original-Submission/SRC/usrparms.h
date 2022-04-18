/*
 -------------------------------------------------------------------------
 File usrparms.h of ADOL-C version 1.6 as of January 1,   1995          
 Included in ---->
                  fos_reverse.c
                  fov_reverse.c 
                  hos_forward.c
                  hos_reverse.c
                  hov_reverse.c
		  taputil1.c
                  taputil2.c
                  taputil3.c
                  tayutil.c
		  utils.c


 ------------------------------------------------------------------------- 

 Usrparms.h contains the parameters which might affect 
 the performance of the ADOL-C system.  Intended to be tweeked by
 users and local maintainence personal.
                       
*/


#define bufsize	    65536 /*16384 or  524288  */
#define locint      unsigned short   
#define revreal     double
#define inf_num     1.0  /* don't undefine these;  on non-IEEE machines */
#define inf_den     0.0  /* change the values to get small fractions    */
#define non_num     0.0  /* (inf_num/inf_den) and (non_num/non_den)     */
#define non_den     0.0  /* respectively, see the documentation         */
/* #define DEBUG       DEBUG  */
#define store       dontusethisuglymessplease 
#define FNAME2       "_adol-rl_tape."
#define FNAME1       "_adol-in_tape."
#define FNAME        "_adol-op_tape."


