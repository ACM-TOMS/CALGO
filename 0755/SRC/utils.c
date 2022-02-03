/*
   --------------------------------------------------------------
   File utils.c of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
   Contains the definitions of trace_on/ trace_off and the 
   class clean_up.  The class clean_up makes sure the once the
   program leaves, any temporary taylor file is deleted.
   Note that this file must be compiled with C++.
*/


/* Basic Includes */

#include "dvlparms.h"
#include <stdio.h>
#include <stream.h>
#include "usrparms.h"

/* External routines from Adouble.c */

extern locint keep_stock();
extern void take_stock();   

extern "C" {
#include "tayutil.h"
#include "taputil3.h"
}

/***********************************************************************/
/* Added class clean-up, so that when the program leaves, it will clean*/
/* up the temporary file.                                              */
/***********************************************************************/

class cleanup{
    int valid;
public:
    cleanup();
    ~cleanup();
};
cleanup::cleanup()
{
  valid = 0;
}
cleanup::~cleanup()
{
      if (taylor_access())
	{
	  close_taylor();
          remove("adoltemp.xxx");     /*     Complies with ANSI standard */ 
      /*     unlink("adoltemp.xxx");   works on some UNIX systems */ 
	}
}

static cleanup at_end;

/***************************************************************************/
/* Trace_on:                                                               */
/* Initialization for the taping process.  Sets up the arrays op_tape,     */
/* int_tape, val_tape, and stats.  Op_tape, int_tape, val_tape are arrays  */
/* of pointers to individual buffers for operations, integers (locints),   */
/* and values (doubles).  Also initializes buffers for this tape, sets     */
/* files names, and calls appropriate setup routines.                      */
/***************************************************************************/

void trace_on(short tnum,int& revals)
{
  start_trace(tnum,revals);
  take_stock();   /* record all existing adoubles on the tape */
}

void trace_on(short tag)
{ 
  int dum = 0; 
  trace_on(tag,dum);
} 

/*************************************************************************/
/* Stop Tracing.  Clean up, and turn off trace_flag.                    **/
/*************************************************************************/
void trace_off(int flag)
{
  int locations;
  locations = keep_stock();     /* copy remaining live variables and turns */
                                /* off trace_flag  */
  stop_trace(locations,flag);   
  cout.flush();
}


/************************************************************************/
/* Trace_off() is essentially trace_off(0).                             */
/************************************************************************/
void trace_off()
{
  int lofl = 0;
  trace_off(lofl);
}

