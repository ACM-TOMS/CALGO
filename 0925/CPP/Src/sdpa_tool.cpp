/*-----------------------------------------
  sdpa_tool.cpp
-----------------------------------------*/

#include "sdpa_tool.h"
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>

#include <unistd.h>
#ifndef CLK_TCK
#define  CLK_TCK  sysconf(_SC_CLK_TCK)
#endif

namespace sdpa {

// These are constant.
// Do Not Change .
int IZERO =  0;
int IONE  =  1;
int IMONE = -1;
double DZERO =  0.0;
double DONE  =  1.0;
double DMONE = -1.0;

double Time::rGetUseTime()
{
  struct tms TIME;
  times(&TIME);
  return (double)TIME.tms_utime/(double)CLK_TCK; 
}

void Time::rSetTimeVal(struct timeval& targetVal)
{
  static struct timezone tz;
  gettimeofday(&targetVal,&tz);
}

double Time::rGetRealTime(const struct timeval& start,
			   const struct timeval& end)
{
  const long int second = end.tv_sec - start.tv_sec;
  const long int usecond = end.tv_usec - start.tv_usec;
  return ((double)second) + ((double)usecond)*(1.0e-6);
}

} // end of namespace 'sdpa'

