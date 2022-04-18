#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "revolve.h"
#include <stdio.h>
#include <stdlib.h>

 int main()
{
  int check, capo, fine, steps, snaps;
  int info;
  enum action whatodo;
  capo = 0;
  snaps = 0;
  printf(" ENTER:   STEPS, SNAPS, INFO \n");
  scanf("%i",&steps);
  scanf("%i",&snaps);
  scanf("%i",&info);
  printf("  %d  %d  %d \n",steps,snaps,info);
  fine = steps + capo;
  check = -1;                  /* Neccessary for first call */
  do
   {
    whatodo = revolve(&check, &capo, &fine, snaps, &info);
    if ((whatodo == takeshot) && (info > 1))
      printf(" takeshot at %6d \n",capo);
    if ((whatodo == advance) && (info > 2))
      printf(" advance to %7d \n",capo);
    if ((whatodo == firsturn) && (info > 2)) 
      printf(" firsturn at %6d \n",capo);
    if ((whatodo == youturn) && (info > 2))
      printf(" youturn at %7d \n",capo);
    if ((whatodo == restore) && (info > 2)) 
      printf(" restore at %7d \n",capo);
    if (whatodo == error) 
     {
      printf(" irregular termination of revolve \n");
      switch(info)
       {
        case 10: printf(" number of checkpoints stored exceeds checkup, \n");
                 printf(" increase constant 'checkup' and recompile \n");
                 break;
        case 11: printf(" number of checkpoints stored = %d exceeds snaps = %d, \n"
                         ,check+1,snaps);
                 printf(" ensure 'snaps' > 0 and increase initial 'fine' \n");
                 break;
        case 12: printf(" error occurs in numforw \n");
                 break;
        case 13: printf(" enhancement of 'fine', 'snaps' checkpoints stored, \n");
                 printf(" increase 'snaps'\n");
                 break;
        case 14: printf(" number of snaps exceeds snapsup, ");
                 printf(" increase constant 'snapsup' and recompile \n");
                 break;
        case 15: printf(" number of reps exceeds repsup, ");
                 printf(" increase constant 'repsup' and recompile \n");       }
     }
   } while((whatodo != terminate) && (whatodo != error));
 return info;
}



