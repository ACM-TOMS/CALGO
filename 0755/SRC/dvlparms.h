/*
 -------------------------------------------------------------------------
 file dvlparms.h of ADOL-C version 1.6 as of January 1,   1995          
 Included in ---->
                  drivers.c
		  driversc.c
                  fos_reverse.c
                  fov_reverse.c 
                  hos_forward.c
                  hos_reverse.c
                  hov_reverse.c
		  taputil1.c
                  taputil2.c
                  taputil3.c
                  tayutil.c


 ------------------------------------------------------------------------- 

  The sole purpose of this file is to provide the developers and          
  maintainers of ADOL-C and include file that contains library wide       
  definitions, etc.                                                       

*/

/* Include the standard library, stdio routines, math functions */
/* and error number routines.                                   */
/* __GNUG__ is defined by g++, __GNUC__ by g++ and gcc          */
#ifdef __GNUG__ 
#include <std.h>
#else
#include <stdlib.h>
#endif

#include <stdio.h>
#include <math.h>
#include <errno.h>

/* Define Compsize */

#define compsize >

/* Possible ADOL-C options. */

#define overwrite   overwrite
#define conditional conditional

/*-------------------------------*/
/* Reserved for future expansion */
/*-------------------------------*/













