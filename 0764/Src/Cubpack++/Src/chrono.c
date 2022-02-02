/////////////////////////////////////////////////////////
//                                                     //
//    Cubpack++                                        //
//                                                     //
//        A Package For Automatic Cubature             //
//                                                     //
//        Authors : Ronald Cools                       //
//                  Dirk Laurie                        //
//                  Luc Pluym                          //
//                                                     //
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//File chrono.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   28 Mar 1996     V0.1h(long instead of int)
/////////////////////////////////////////////////////////

#include <chrono.h>
#ifdef GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#endif

/////////////////////////////////////////////////////////
Chrono::Chrono()
  {
  Reset();
  }
/////////////////////////////////////////////////////////
void
Chrono::Start()
  {
  if ( !Running)
    {
    Time= OldTime = times_();
    Running =True;
    };
  }
/////////////////////////////////////////////////////////
void
Chrono::Stop()
  {
  if (Running)
    {
    Time = times_();
    Running =False;
    };
  }
/////////////////////////////////////////////////////////
void
Chrono::Reset()
  {
  Time=0;
  OldTime=0;
  Running = False;
  }
/////////////////////////////////////////////////////////
unsigned long
Chrono::Read()
  {
  if (Running)
    {
    return (times_() - OldTime);
    }
  else
    {
    return (Time -OldTime);
    };
  }
/////////////////////////////////////////////////////////
long
Chrono::times_ ()
   {
   long t=0;

#ifdef GETRUSAGE
   struct rusage  info ;
   struct timeval *user_t = &(info.ru_utime), *sys_t  = &(info.ru_stime) ;

   getrusage (0, &info) ;

   t += user_t->tv_sec * 1000 + user_t->tv_usec / 1000 ;
   t += sys_t->tv_sec  * 1000 + sys_t->tv_usec  / 1000 ;
#endif
return t;
   }
/////////////////////////////////////////////////////////
