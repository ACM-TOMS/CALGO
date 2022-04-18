//int getrusage(int, struct rusage *);
//int gettimeofday(struct timeval *, struct timezone *);

//*************************************************************************
// function: timer()     
//                                                                        
// Programmer: Karen Minser                  
// Date: September 24, 1993                                               
//                                                                        
// function:  timer()                                                 
//                                                                     
// User-supplied function returns elapsed cpu time (float)             
//*************************************************************************

#include <sys/time.h>
#include "resource.h"

double timer()
{
   double elasped_time;
   struct rusage mytime;
   getrusage(RUSAGE_SELF,&mytime);
   						// convert elapsed time to milliseconds 
   						// and return elapsed time in seconds 
   elasped_time =  ((double)mytime.ru_utime.tv_sec * 1000000 +
                  (double)mytime.ru_utime.tv_usec);
   return(elasped_time / 1000000);

}  // end of timer 


//*************************************************************************
// function: timer2()     
//                                                                        
// Programmer: Karen Minser                  
// Date: September 24, 1993                                               
//                                                                        
// function:  timer()                                                 
//                                                                     
// This returns a time stamp for wall clock time (double)             
//*************************************************************************

double timer2()
{
   struct timeval timval1;
   struct timezone timzon;
   double t;

   gettimeofday(&timval1,&timzon);
   t = timval1.tv_sec + timval1.tv_usec/1000000.;
   return t;

}     // end of timer2 
