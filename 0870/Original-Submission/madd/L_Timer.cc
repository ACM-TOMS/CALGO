////////////////////////////////////////////////////
//
// Timer class 
//
// Author: Leonidas Linardakis
//
// (c) 2001 Leonidas Linardakis
//
//
////////////////////////////////////////////////////

using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include "math.h"
#include "L_Timer.h"

void
L_Timer::setSeconds(double secs)
{
   seconds = secs;
}


void 
L_Timer::printExpanded()
{
   expandFromSeconds();
   cout << endl
      << tmStruc.year << " years; "
      << tmStruc.week << " weeks; "
      << tmStruc.day  << " days; "
      << tmStruc.hour << " hours; "
      << tmStruc.min  << " mins; "
      << tmStruc.sec  << " secs and "
      << tmStruc.fraction << " of secs."
      << endl;
}

       
void 
L_Timer::expandFromSeconds()
{  
  double eval; 
  tmStruc.fraction = modf(seconds, &eval);
  tmStruc.sec = (long unsigned)eval % 60;
  eval -= tmStruc.sec;
  eval /= 60;
  tmStruc.min = (long unsigned)(eval) % 60;
  eval -= tmStruc.min;
  eval /= 60;
  tmStruc.hour = (long unsigned)(eval) % 24;
  eval -= tmStruc.hour;
  eval /= 24;
  
  // now we have days
  int tmpDay = (long unsigned)(eval) % 365;
  eval -= tmpDay;
  tmStruc.year = eval / 365;
  tmStruc.day = tmpDay % 7;
  tmStruc.week = (tmpDay - tmStruc.day) / 7;
}


   
void 
L_Timer::start()
{
   theClock = clock();
}


long  double 
L_Timer::stop()
{
   endClock = clock();
   return (long double)(endClock - theClock) / (long double)CLOCKS_PER_SEC; 
}



long double 
L_Timer::stopAndPrint(char* filename)
{
   endClock = clock();
   ofstream oFile;
   oFile.open(filename, ios::app);
   long double cyc = ((long double)(endClock) - (long double)(theClock));
   long double sec = cyc / (long double)(CLOCKS_PER_SEC);
   oFile << endl << cyc << "   " << sec;
   return sec; 
}


long double 
L_Timer::cycles()
{
   return ((long double)(endClock) - (long double)(theClock)); 
}


long  double 
L_Timer::clockInSecs()
{
   return (long double)((long long double)clock() / (long long double)(CLOCKS_PER_SEC));
}
   

L_Timer::L_Timer()
{
  
}
