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
 * profile.c                                                                *
 ****************************************************************************/

#include "profile.h"


/* Main structure for clocks */
struct IBClockElem IBClock[IBNbClock];


void IBClockInit()
/***************************************************************************
*  Initialization of clocks before utilisation
*/
{
  int i;
  for( i=0; i<IBNbClock; i++ )
  {
    IBClock[i].accu = 0;
  }
}


inline void IBClockFree()
/***************************************************************************
*  Nothing to do !
*/
{
}


inline void IBClockBegin(int i)
/***************************************************************************
*  Beginning of timed operation number `i'
*/
{
  long usrtime;

  IBGetTimes(usrtime);
  IBClock[i].btime = usrtime;
}


inline void IBClockEnd(int i)
/***************************************************************************
*  End of timed operation number `i'
*  Amount of time between `btime' and now is summed in `accu'
*/
{
  long usrtime;

  IBGetTimes(usrtime);
  IBClock[i].accu += (usrtime - IBClock[i].btime);
}


inline long IBClockGet(int i)
/***************************************************************************
*  Returns the value in `accu'
*/
{
  return( IBClock[i].accu );
}


inline void IBClockSet(int i, long t)
/***************************************************************************
*  Set the time for the i-th operation
*/
{
  IBClock[i].accu = t;
}


inline long IBClockObserve (int i)
/***************************************************************************
*  Elapsed time since the last IBClockBegin(i)
*/
{
  long usrtime;

  IBGetTimes(usrtime);
  return( usrtime - IBClock[i].btime );
}
