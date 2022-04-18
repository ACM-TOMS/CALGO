/*******************************************************************************
***                                                                          ***
***  This program is furnished by the U.S. Army Engineer Research and        ***
***  Development Center, Major Shared Resource Center (ERDC MSRC) "as is"    ***
***  and is accepted and used by the recipient with the express              ***
***  understanding that the Government makes no warranties, expressed or     ***
***  implied, concerning the accuracy, completeness, reliability, usability  ***
***  or suitability for any particular purpose of the information and data   ***
***  within this program or furnished in connection therewith, and the       ***
***  Government shall be under no liability whatsoever to any person by      ***
***  reason of any use made thereof. This program belongs to the U.S.        ***
***  Government; therefore, the recipient further agrees not to assert any   ***
***  proprietary rights therein or to represent the source code to anyone    ***
***  as other than a Government program.                                     ***
***                                                                          ***
*******************************************************************************/
 
/*---------------------------------------------------------------------------*\
 |                                                                            |
 |  Authors:                                                                  |
 |  Richard J. Hanson (koolhans@rice.edu)                                     |
 |  Rice University, Center for High Performance Software Research            |
 |                                                                            |
 |  Clay P. Breshears (clay.breshears@intel.com)                              |
 |  KAI Software, a division of Intel Americas, Inc.                          |
 |                                                                            |
 |  Henry A. Gabb (henry.gabb@intel.com)                                      |
 |  KAI Software, a division of Intel Americas, Inc.                          |
\*---------------------------------------------------------------------------*/
 
 
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
/* This code is part of the package "A Fortran Interface to Posix Threads,"  to be
   published in ACM-TOMS.  Authors: R. Hanson, C. Breshears, and H. Gabb.
   This is summary.h. Last change on 4 December 2000. */
 
#ifdef SGI
#define PTHREAD_STACK_MIN    _sysconf(_SC_THREAD_STACK_MIN)
#define PTHREAD_THREADS_MAX  1024 /*_sysconf(_SC_THREADS_MAX)*/
#endif
 
#ifdef SUN
/* #include <thread_db.h> /* Needed to get defn of sigset_t */
#include <sys/signal.h> /* Needed to get defn of sigset_t */
#define PTHREAD_STACK_MIN  _sysconf(_SC_THREAD_STACK_MIN)
#define PTHREAD_THREADS_MAX _POSIX_THREAD_THREADS_MAX
#define PTHREAD_KEYS_MAX    _POSIX_THREAD_KEYS_MAX
#define  PTHREAD_PRIO_INHERIT 2
#define  PTHREAD_PRIO_PROTECT 1
#define  PTHREAD_PRIO_NONE    0
#define PTHREAD_DESTRUCTOR_ITERATIONS 4
#endif
 
#ifdef IBM
#define  PTHREAD_CREATE_JOINABLE 0
#define  PTHREAD_PROCESS_PRIVATE 0
#define  PTHREAD_PROCESS_SHARED 0
#define  PTHREAD_PRIO_INHERIT 2
#define  PTHREAD_PRIO_PROTECT 1
#define  PTHREAD_PRIO_NONE    0
#define PTHREAD_DESTRUCTOR_ITERATIONS 4
#endif
 
#ifdef CPQ
#define PTHREAD_THREADS_MAX _POSIX_THREAD_THREADS_MAX
#define  PTHREAD_PRIO_INHERIT 2
#define  PTHREAD_PRIO_PROTECT 1
#define  PTHREAD_PRIO_NONE    0
#endif
 
#ifdef LINUX
#define  PTHREAD_PRIO_INHERIT 2
#define  PTHREAD_PRIO_PROTECT 1
#define  PTHREAD_PRIO_NONE    0
#endif
/*
 * NOTE: the "do {" ... "} while (0);" bracketing around the macros
 * allows the various printing and abort macros to be used as if they
 * were function calls, even in contexts where a trailing ";" would
 * generate a null statement. For example,
 *
 *      if (status != 0)
 *          err_abort (status, "message");
 *      else
 *          return status;
 *
 * will not compile if err_abort is a macro ending with "}", because
 * C does not expect a ";" to follow the "}". Because C does expect
 * a ";" following the ")" in the do...while construct, err_abort and
 * errno_abort can be used as if they were function calls.
 */
#define test_comment(test_number,text)\
do {fprintf (stdout, "Testing %3d%s at \"%s\", line %d\n",test_number, text, __FILE__, __LINE__);} while(0)
#define err_abort(test_number,code,text) do { \
    if(code != 0)fprintf (stdout, "Failed %4d %s at \"%s\":%d: Status= %3d\n", \
        test_number, text, __FILE__, __LINE__, code); \
    if(code != 0) abort ();\
    } while (0)
 
#define skip do{fprintf (stdout,"\n");} while(0)
