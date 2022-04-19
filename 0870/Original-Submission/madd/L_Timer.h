////////////////////////////////////////////////////
//
// Timer class 
//
// Author: Leonidas Linardakis
//
//
// (c) 2001 Leonidas Linardakis
//
////////////////////////////////////////////////////

#ifndef __L_Timer
#define __L_Timer


#include <time.h>

struct L_STRUC_time{
   double fraction;
   int sec;
   int min;
   int hour;
   int day;
   int week;
   double year; 
};

class L_Timer{

public:
   
        L_Timer(); // constructor 
   void printExpanded();
   void setSeconds(double secs);
   void start();
   long double stop();
   long double stopAndPrint(char* filename);
   long double clockInSecs();
   long double cycles();

private:

   void expandFromSeconds();
                               
   long double seconds;
   L_STRUC_time tmStruc;
   clock_t theClock, endClock;

};

#endif



