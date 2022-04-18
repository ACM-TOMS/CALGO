/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Authors:
!!!  Clay P. Breshears (clay.breshears@intel.com)
!!!  KAI Software, a division of Intel Americas, Inc.
!!!
!!!  Henry A. Gabb (henry.gabb@intel.com)
!!!  KAI Software, a division of Intel Americas, Inc.
!!! 
!!!  Richard J. Hanson (koolhans@rice.edu)
!!!  Rice University,  Rice Center for High Performance Software Research
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#include <pthread.h>
#include <limits.h>
#include <stdio.h>
#include<string.h>
#include <errno.h>
#include <unistd.h>

/* 
   Estimates of integer array sizes needed to hold data from C structures.  
   Configuration will use maximum of this value and actual system size
   to ensure enough room is available.

#define ISIZE_OF_C_INT               1
#define ISIZE_OF_fpthrd_T            2
#define ISIZE_OF_fpthrd_ATTR_T      44
#define ISIZE_OF_fpthrd_MUTEX_T     16
#define ISIZE_OF_fpthrd_MUTEXATTR_T 32
#define ISIZE_OF_fpthrd_COND_T      16
#define ISIZE_OF_fpthrd_CONDATTR_T  28
#define ISIZE_OF_FSCHED_PARAM        9
#define ISIZE_OF_FTIMESPEC           4
#define ISIZE_OF_FSIGSET_T           4
#define ISIZE_OF_FSIZE_T             2

#define max(x,y) (x>y)?x:y

*/

void usage()
{
	printf("\nPrint Usage example here.\n\n");
}

void fpthrd_type_decl(FILE *fp, 
                      char fsize_name[], 
                      int fsize, 
                      char ftype_name[], 
                      char in_name[])
{
   char int_param[30] = "      INTEGER, PARAMETER :: ";

   fprintf(fp,"%s %s = %d\n",int_param, fsize_name, fsize);

   fprintf(fp,"      TYPE %s\n", ftype_name);
   fprintf(fp,"         PRIVATE\n");
   fprintf(fp,"         INTEGER(KIND=IADDR) :: %s(%s)\n", in_name, fsize_name);
   fprintf(fp,"      END TYPE %s\n\n", ftype_name);
}

void fpthrd_stype_decl(FILE *fp, 
                       char fsize_name[], 
                       int fsize, 
                       char ftype_name[], 
                       char in_name[])
{
   char int_param[30] = "      INTEGER, PARAMETER :: ";

   fprintf(fp,"%s %s = %d\n",int_param, fsize_name, fsize);

   fprintf(fp,"      TYPE %s\n", ftype_name);
   fprintf(fp,"         INTEGER :: %s(%s)\n", in_name, fsize_name);
   fprintf(fp,"      END TYPE %s\n\n", ftype_name);
}

void get_init_values(FILE *fp, void *t, int n, const char *what)
{
   int *i;
   int j;

   fprintf(fp,"      INTEGER (KIND=IADDR), DIMENSION(%d), PARAMETER ::",n);
   fprintf(fp," %s_I = &\n",what);
   fprintf(fp,"      (/ ");

printf("fpthrd_config: Size of %s INIT is %d x 4 bytes\n",what,n);

   i = (int *)t;
	fprintf(fp, "%d",*i++);
   for (j = 1; j < n; j++) {
		fprintf(fp, ", %d",*i++);
	}
   fprintf(fp," /)\n\n");
}

void get_long_init_values(FILE *fp, void *t, int n, const char *what)
{
   long *i;
   int j;

   fprintf(fp,"      INTEGER (KIND=IADDR), DIMENSION(%d), PARAMETER ::",n);
   fprintf(fp," %s_I = &\n",what);
   fprintf(fp,"      (/ ");

printf("fpthrd_config: Size of %s INIT is %d x 8 bytes\n",what,n);

   i = (long *)t;
	fprintf(fp, "%ld",*i++);
   for (j = 1; j < n; j++) {
		fprintf(fp, ", %ld",*i++);
	}
   fprintf(fp," /)\n\n");
}

int findSize(int s, int n)
{
	if (s%n == 0) 
		return(s/n);
	else
		return((s/n)+1);
}

main(int argc, char *argv[])
{
   FILE *fp;
   char int_param[30] = "      INTEGER, PARAMETER :: ";
   int bit64 = 0;  /* Is this a 64-bit addressin platform? */
   int numbytes = 4;

   pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
   pthread_cond_t  c = PTHREAD_COND_INITIALIZER;

   if (argc!=2) {
		printf("fpthrd_config: Incorrect number of command line arguments.\n");
      usage();
		return 1;
	}

   if ((fp = fopen("constants.inc", "w")) == NULL) {
      printf("fpthrd_config: Unable to open file 'constants.inc'\n");
      return 1;
   }

/* FPTHRD Internal data and address sizes */

   if (strcmp(argv[1],"SUN") == 0){
		fprintf(fp,"%s RR=6\n",int_param);
	}
	else if (strcmp(argv[1],"IBM") == 0) {
		fprintf(fp,"%s RR=6\n",int_param);
	}
	else if (strcmp(argv[1],"SGI") == 0) {
		fprintf(fp,"%s RR=6\n",int_param);
	}
	else if (strcmp(argv[1],"SGI64") == 0) {
      bit64 = 1;
      numbytes = 8;
		fprintf(fp,"%s RR=12\n",int_param);
	}
 	else if (strcmp(argv[1],"CPQ") == 0) {
      bit64 = 1;
      numbytes = 8;
		fprintf(fp,"%s RR=12\n",int_param);
	}
	else {
		printf("fpthrd_config: Invalid command line argument.\n");
		fclose(fp);
     	usage();
		return 1;
	}

   fprintf(fp,"%s IADDR = SELECTED_INT_KIND(R=RR)\n",int_param);
   fprintf(fp,"      INTEGER (KIND=IADDR) gauge_IADDR\n\n");

	if (sizeof(int) == 8) {
		printf("fpthrd_config: C integer size is 8 bytes.\n");
      fprintf(fp,"%s CINT = SELECTED_INT_KIND(R=12)\n",int_param);
	}
	else {
		if (sizeof(int) == 4)
			printf("fpthrd_config: C integer size is 4 bytes.\n");
		else {
			printf("fpthrd_config: C integer size is %ld bytes.\n",sizeof(int));
			printf("fpthrd_config: *Setting size to be 4 bytes in FPTHRD.*\n");
		}
      fprintf(fp,"%s CINT = SELECTED_INT_KIND(R=6)\n",int_param);
	}
   fprintf(fp,"      INTEGER (KIND=CINT) gauge_CINT\n\n");
	

/* POSIX thread constants */

   fprintf(fp,"%s FPTHRD_CREATE_DETACHED = ",int_param);
   fprintf(fp,"%d\n",PTHREAD_CREATE_DETACHED);

   fprintf(fp,"%s FPTHRD_CREATE_JOINABLE = ",int_param);
   fprintf(fp,"%d\n\n",PTHREAD_CREATE_JOINABLE);

#ifdef _POSIX_THREAD_PROCESS_SHARED
   fprintf(fp,"%s FPTHRD_PROCESS_PRIVATE = ",int_param);
   fprintf(fp,"%d\n",PTHREAD_PROCESS_PRIVATE);

   fprintf(fp,"%s FPTHRD_PROCESS_SHARED = ",int_param);
   fprintf(fp,"%d\n\n",PTHREAD_PROCESS_SHARED);
#else
   fprintf(fp,"%s FPTHRD_PROCESS_PRIVATE = 0\n",int_param);
   fprintf(fp,"%s FPTHRD_PROCESS_SHARED = 0\n\n",int_param);
#endif 

#ifdef _POSIX_THREAD_PRIO_PROTECT
   fprintf(fp,"%s FPTHRD_PRIO_PROTECT = ",int_param);
   fprintf(fp,"%d\n",PTHREAD_PRIO_PROTECT);

   fprintf(fp,"%s FPTHRD_PRIO_INHERIT = ",int_param);
   fprintf(fp,"%d\n",PTHREAD_PRIO_INHERIT);

   fprintf(fp,"%s FPTHRD_PRIO_NONE = ",int_param);
   fprintf(fp,"%d\n\n",PTHREAD_PRIO_NONE);
#else
   fprintf(fp,"%s FPTHRD_PRIO_PROTECT = 0\n",int_param);
   fprintf(fp,"%s FPTHRD_PRIO_INHERIT = 0\n",int_param);
   fprintf(fp,"%s FPTHRD_PRIO_NONE = 0\n\n",int_param);
#endif

   fprintf(fp,"%s FPTHRD_CANCEL_ENABLE = ",int_param);
   fprintf(fp,"%d\n",PTHREAD_CANCEL_ENABLE);

   fprintf(fp,"%s FPTHRD_CANCEL_DISABLE = ",int_param);
   fprintf(fp,"%d\n\n",PTHREAD_CANCEL_DISABLE);

   fprintf(fp,"%s FPTHRD_CANCEL_DEFERRED = ",int_param);
   fprintf(fp,"%d\n",PTHREAD_CANCEL_DEFERRED);

   fprintf(fp,"%s FPTHRD_CANCEL_ASYNCHRONOUS = ",int_param);
   fprintf(fp,"%d\n\n",PTHREAD_CANCEL_ASYNCHRONOUS);

   fprintf(fp,"      INTEGER (KIND=CINT), PARAMETER :: ");
   fprintf(fp,"FPTHRD_CANCELED = huge(1)-1\n\n");

#ifdef _POSIX_THREAD_PRIORITY_SCHEDULING
   fprintf(fp,"%s FPTHRD_SCOPE_SYSTEM = ",int_param);
   fprintf(fp,"%d\n",PTHREAD_SCOPE_SYSTEM);

   fprintf(fp,"%s FPTHRD_SCOPE_PROCESS = ",int_param);
   fprintf(fp,"%d\n\n",PTHREAD_SCOPE_PROCESS);
#else
   fprintf(fp,"%s FPTHRD_SCOPE_SYSTEM = 0\n",int_param);
   fprintf(fp,"%s FPTHRD_SCOPE_PROCESS = 0\n\n",int_param);
#endif

#ifdef _POSIX_THREAD_PRIORITY_SCHEDULING
   fprintf(fp,"%s FPTHRD_INHERIT_SCHED = ",int_param);
   fprintf(fp,"%d\n",PTHREAD_INHERIT_SCHED);

   fprintf(fp,"%s FPTHRD_EXPLICIT_SCHED = ",int_param);
   fprintf(fp,"%d\n\n",PTHREAD_EXPLICIT_SCHED);
#else
   fprintf(fp,"%s FPTHRD_INHERIT_SCHED = 0\n",int_param);
   fprintf(fp,"%s FPTHRD_EXPLICIT_SCHED = 0\n\n",int_param);
#endif

#ifdef PTHREAD_KEYS_MAX
   fprintf(fp,"%s FPTHRD_KEYS_MAX = ",int_param);
   fprintf(fp,"%d\n",PTHREAD_KEYS_MAX);
#else
   printf("fpthrd_config: Using sysconf() value for FPTHRD_KEYS_MAX\n");
   fprintf(fp,"%s FPTHRD_KEYS_MAX = ",int_param);
   fprintf(fp,"%d\n",(int)sysconf(_SC_THREAD_KEYS_MAX));
#endif

#ifdef PTHREAD_STACK_MIN
   fprintf(fp,"%s FPTHRD_STACK_MIN = ",int_param);
   fprintf(fp,"%d\n",PTHREAD_STACK_MIN);
#else
   printf("fpthrd_config: Using sysconf() value for FPTHRD_STACK_MIN\n");
   fprintf(fp,"%s FPTHRD_STACK_MIN = ",int_param);
   fprintf(fp,"%d\n",(int)sysconf(_SC_THREAD_STACK_MIN));
#endif

#ifdef PTHREAD_THREADS_MAX
   fprintf(fp,"%s FPTHRD_THREADS_MAX = ",int_param);
   fprintf(fp,"%d\n",PTHREAD_THREADS_MAX);
#else
   printf("fpthrd_config: Using sysconf() value for FPTHRD_THREADS_MAX\n");
   fprintf(fp,"%s FPTHRD_THREADS_MAX = ",int_param);
   fprintf(fp,"%d\n",(int)sysconf(_SC_THREAD_THREADS_MAX));
#endif

#ifdef PTHREAD_DESTRUCTOR_ITERATIONS
   fprintf(fp,"%s FPTHRD_DESTRUCTOR_ITERATIONS = ",int_param);
   fprintf(fp,"%d\n\n",PTHREAD_DESTRUCTOR_ITERATIONS);
#else
   printf("fpthrd_config: Using sysconf() value for FPTHRD_DESTURCTOR_ITERATIONS\n");
   fprintf(fp,"%s FPTHRD_DESTRUCTOR_ITERATIONS = ",int_param);
   fprintf(fp,"%d\n\n",(int)sysconf(_SC_THREAD_DESTRUCTOR_ITERATIONS));
#endif

#ifdef _POSIX_THREAD_PRIORITY_SCHEDULING
   fprintf(fp,"%s FSCHED_FIFO = %d\n",int_param,SCHED_FIFO);

   fprintf(fp,"%s FSCHED_RR = %d\n",int_param,SCHED_RR);

   fprintf(fp,"%s FSCHED_OTHER = %d\n\n",int_param,SCHED_OTHER);
#else
   fprintf(fp,"%s FSCHED_FIFO = 0\n",int_param);
   fprintf(fp,"%s FSCHED_RR = 0\n",int_param);
   fprintf(fp,"%s FSCHED_OTHER = 0\n\n",int_param);
#endif


/* POSIX thread error codes */

   fprintf(fp,"%s ESRCH = %d\n",int_param, ESRCH);
   fprintf(fp,"%s EINVAL = %d\n",int_param, EINVAL);
   fprintf(fp,"%s EFAULT = %d\n",int_param, EFAULT);
   fprintf(fp,"%s ENOTSUP = %d\n",int_param, ENOTSUP);
   fprintf(fp,"%s EAGAIN = %d\n",int_param, EAGAIN);
   fprintf(fp,"%s EDEADLK = %d\n",int_param, EDEADLK);
   fprintf(fp,"%s ENOSYS = %d\n",int_param, ENOSYS);
   fprintf(fp,"%s EPERM = %d\n",int_param, EPERM);
   fprintf(fp,"%s EBUSY = %d\n",int_param, EBUSY);
   fprintf(fp,"%s ENOMEM = %d\n",int_param, ENOMEM);
   fprintf(fp,"%s ETIMEDOUT = %d\n",int_param, ETIMEDOUT);
   fprintf(fp,"%s EINTR = %d\n",int_param, EINTR);
   fprintf(fp,"%s ENOSPC = %d\n",int_param, ENOSPC);

   fprintf(fp,"\n%s NUMBER_OF_ERRORS = 13\n",int_param);
   fprintf(fp,"      integer(KIND=CINT), private ::");
   fprintf(fp," fpthrd_errors(NUMBER_OF_ERRORS) = &\n");
   fprintf(fp,"      (/ %d, ",ESRCH);
   fprintf(fp,"%d, ",EINVAL);
   fprintf(fp,"%d, ",EFAULT);
   fprintf(fp,"%d, ",ENOTSUP);
   fprintf(fp,"%d, ",EAGAIN);
   fprintf(fp,"%d, ",EDEADLK);
   fprintf(fp,"%d, ",ENOSYS);
   fprintf(fp,"%d, ",EPERM);
   fprintf(fp,"%d, ",EBUSY);
   fprintf(fp,"%d, ",ENOMEM);
   fprintf(fp,"%d, ",ETIMEDOUT);
   fprintf(fp,"%d, ",EINTR);
   fprintf(fp,"%d /)\n\n",ENOSPC);


   fprintf(fp,"%s ISIZE_OF_C_INT = %d\n\n",int_param,1);

/*
        INTEGER               ISIZE_OF_ADDRESS
*/

/* POSIX thread data types.  Some contain pointers to C structures. */


   fpthrd_type_decl( fp,
                    "ISIZE_OF_fpthrd_T", 
                     findSize(sizeof(pthread_t),numbytes),
                    "fpthrd_t",
                    "thread");

   fpthrd_type_decl( fp,
                    "ISIZE_OF_fpthrd_ATTR_T",
                     findSize(sizeof(pthread_attr_t),numbytes),
                    "fpthrd_attr_t",
                    "threadattr");

   fpthrd_type_decl( fp,
                    "ISIZE_OF_fpthrd_MUTEX_T",
                     findSize(sizeof(pthread_mutex_t),numbytes),
                    "fpthrd_mutex_t",
                    "mutex");

   if (!bit64) 
		get_init_values(fp,&m,sizeof(m)/4,"MUTEX");
	else
		get_long_init_values(fp,&m,sizeof(m)/8,"MUTEX");
	

   fpthrd_type_decl( fp,
                    "ISIZE_OF_fpthrd_MUTEXATTR_T",
                     findSize(sizeof(pthread_mutexattr_t),numbytes),
                    "fpthrd_mutexattr_t",
                    "mutexattr");

   fpthrd_type_decl( fp,
                    "ISIZE_OF_fpthrd_COND_T",
                     findSize(sizeof(pthread_cond_t),numbytes),
                    "fpthrd_cond_t",
                    "conditional");

   if (!bit64) 
		get_init_values(fp,&c,sizeof(c)/4,"COND");
	else
		get_long_init_values(fp,&c,sizeof(c)/8,"COND");
	
   fpthrd_type_decl( fp,
                    "ISIZE_OF_fpthrd_CONDATTR_T",
                     findSize(sizeof(pthread_condattr_t),numbytes),
                    "fpthrd_condattr_t",
                    "condattr");

   fpthrd_stype_decl( fp,
                     "ISIZE_OF_FSCHED_PARAM",
                      findSize(sizeof(struct sched_param),numbytes),
                     "fsched_param",
                     "sched_priority");

   fpthrd_stype_decl( fp,
                     "ISIZE_OF_FTIMESPEC",
                      findSize(sizeof(struct timespec),numbytes),
                     "ftimespec",
                     "timespec");

   fpthrd_stype_decl( fp,
                     "ISIZE_OF_FSIGSET_T",
                      findSize(sizeof(sigset_t),numbytes),
                     "fsigset_t",
                     "mask");

   fprintf(fp,"      INTEGER, PARAMETER :: ISIZE_OF_FSIZE_T = ");
   fprintf(fp,"%d\n",findSize(sizeof(size_t),numbytes));
   fprintf(fp,"      TYPE fsize_t\n");
   fprintf(fp,"         INTEGER (kind=IADDR) :: size(ISIZE_OF_FSIZE_T)\n");
   fprintf(fp,"      END TYPE\n");

/* Close file */

   fclose(fp);

   printf("FPTHRD Configuration complete\n");
}

