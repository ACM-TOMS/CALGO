/****************************************************************************
 * RealPaver v. 0.4                                                         *
 *--------------------------------------------------------------------------*
 * Author: Laurent Granvilliers                                             *
 * Copyright (c) 1999-2003 Institut de Recherche en Informatique de Nantes  *
 * Copyright (c) 2004      Laboratoire d'Informatique de Nantes Atlantique  *
 *--------------------------------------------------------------------------*
 * RealPaver is distributed WITHOUT ANY WARRANTY. Read the associated       *
 * COPYRIGHT file for more details.                                         *
 *--------------------------------------------------------------------------*
 * profile.h                                                                *
 ****************************************************************************/

#ifndef __profile_h
#define __profile_h

#include <stdlib.h>
#include <time.h>
#include <sys/resource.h>
#include "config.h"


/****************************************************************************
 *                                  SPARC                                   *
 ****************************************************************************/
#if SYSTEM_SPARC

struct rusage IBRsrUsage;

#define IBGetTimes(usr_time)                                        \
    getrusage(RUSAGE_SELF,&IBRsrUsage);                             \
                                                                    \
    usr_time=IBRsrUsage.ru_utime.tv_sec*1000 +                      \
               IBRsrUsage.ru_utime.tv_usec/1000

/****************************************************************************
 *                             PC i386 & linux                              *
 ****************************************************************************/
#elif SYSTEM_LINUX_IX86

#define IBGetTimes(usr_time) \
       usr_time = (long)(((clock()/(double)CLOCKS_PER_SEC))*1000.0)


/****************************************************************************
 *                                MIPS SGI                                  *
 ****************************************************************************/
#elif SYSTEM_SGI

#define IBGetTimes(usr_time) \
       usr_time = (long)(((clock()/(double)CLOCKS_PER_SEC))*1000.0)

#endif
/****************************************************************************/


#define IBNbClock 10

/* One identifier per timed operation */
#define IBClockParse  0    /* parsing */
#define IBClockSelect 1    /* unused */
#define IBClockSolve  2    /* total solving time */
#define IBClockBC3    3    /* solving with algorithm BC3 */
#define IBClockHC4    4    /* solving with algorithm HC4 */
#define IBClockINwt   5    /* solving with algorithm Newton */
#define IBClockHC3    6    /* solving with algorithm HC3 */


/* clock type */
struct IBClockElem
{
  long btime;   /* starting time of operation which is timed */
  long accu;    /* total accumulated time */
};

/* clock functions */
void IBClockInit    ();               /* initialization for this module */
void IBClockFree    ();               /* desallocation of memory to be called after use */
void IBClockBegin   (int i);          /* starting time for i-th clock */
void IBClockEnd     (int i);          /* ending time for i-th clock */
long IBClockGet     (int i);          /* value of i-th clock */
void IBClockSet     (int i,long t);   /* value of i-th clock */
long IBClockObserve (int i);          /* elapsed time since the last IBClockBegin(i) */

#endif
