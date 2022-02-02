/* ******************************************************************************
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
****************************************************************************** */

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


/* This code is part of the package "A Fortran Interface to Posix Threads,"
  to be published in ACM-TOMS.  Authors: R. Hanson, C. Breshears, and H. Gabb.
   This is ptf90.c.
   Last change:  CPB   Tue Jan 22 11:40:46 CST 2002
*/
#ifdef LINUX
#  define _REENTRANT
#  define _POSIX_SOURCE
#endif

#include <pthread.h>
#include <limits.h>
#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include "summary.h"

/* Long integers are required for 64-bit memory addresses: */
#ifdef SGI
#if BIT64 == SGI
#define INT_CAST long long
#else
#define INT_CAST int
#endif
#else
#define INT_CAST int
#endif

#define PTHRD_INT  int

/* 8-byte integer paramters are required for -i8 compilation */

#ifdef SGI
#if SGI == I8
#define PARAM_INT long long
#else
#define PARAM_INT int
#endif
#else
#define PARAM_INT int
#endif

/* This value forces an abort if the Fortran code has not called fpthrd_data_exchange(). */

#define FORTRAN_NULL       0x7FFFFFFF
#define FPTHREAD_CANCELED  0x7FFFFFFE



/* IBM : The IBM Fortran compiler folds symbol names to lower case and does not add
         extra characters to the end.   */
#ifdef NO_CHANGE
#define fpthrd_attr_init                 fpthr_attr_init
#define fpthrd_attr_destroy              fpthr_attr_destroy
#define fpthrd_attr_setstacksize         fpthr_attr_setstacksize
#define fpthrd_attr_getstacksize         fpthr_attr_getstacksize
#define fpthrd_attr_setdetachstate       fpthr_attr_setdetachstate
#define fpthrd_attr_getdetachstate       fpthr_attr_getdetachstate
#define fpthrd_attr_setscope             fpthr_attr_setscope
#define fpthrd_attr_getscope             fpthr_attr_getscope
#define fpthrd_attr_setinheritsched      fpthr_attr_setinheritsched
#define fpthrd_attr_getinheritsched      fpthr_attr_getinheritsched
#define fpthrd_attr_setschedpolicy       fpthr_attr_setschedpolicy
#define fpthrd_attr_getschedpolicy       fpthr_attr_getschedpolicy
#define fpthrd_attr_setschedparam        fpthr_attr_setschedparam
#define fpthrd_attr_getschedparam        fpthr_attr_getschedparam
#define fpthrd_create                    fpthrd_create
#define fpthrd_join                      fpthr_join
#define fpthrd_exit                      fpthrd_exit
#define fpthrd_detach                    fpthr_detach
#define fpthr_self                       fpthr_self
#define fpthrd_self                      fpthrd_self
#define fpthr_equal                      fpthr_equal
#define fpthrd_equal                     fpthrd_equal
#define fpthrd_getschedparam             fpthr_getschedparam
#define fpthrd_setschedparam             fpthr_setschedparam
#define fpthrd_cancel                    fpthr_cancel
#define fpthrd_setcancelstate            fpthr_setcancelstate
#define fpthrd_setcanceltype             fpthr_setcanceltype
#define fpthr_testcancel                 fpthr_testcancel
#define fpthrd_testcancel                fpthrd_testcancel
#define fpthrd_cleanup_push              fpthrd_cleanup_push
#define fpthrd_cleanup_pop               fpthrd_cleanup_pop
#define fpthrd_mutexattr_init            fpthr_mutexattr_init
#define fpthrd_mutexattr_destroy         fpthr_mutexattr_destroy
#define fpthrd_mutexattr_getpshared      fpthr_mutexattr_getpshared
#define fpthrd_mutexattr_setpshared      fpthr_mutexattr_setpshared
#define fpthrd_mutexattr_setprotocol     fpthr_mutexattr_setprotocol
#define fpthrd_mutexattr_getprotocol     fpthr_mutexattr_getprotocol
#define fpthrd_mutexattr_setprioceiling  fpthr_mutexattr_setprioceiling
#define fpthrd_mutexattr_getprioceiling  fpthr_mutexattr_getprioceiling
#define fpthrd_mutex_init                fpthr_mutex_init
#define fpthrd_mutex_destroy             fpthr_mutex_destroy
#define fpthrd_mutex_lock                fpthr_mutex_lock
#define fpthrd_mutex_trylock             fpthr_mutex_trylock
#define fpthrd_mutex_unlock              fpthr_mutex_unlock
#define fpthrd_mutex_setprioceiling      fpthr_mutex_setprioceiling
#define fpthrd_mutex_getprioceiling      fpthr_mutex_getprioceiling
#define fpthrd_condattr_init             fpthr_condattr_init
#define fpthrd_condattr_destroy          fpthr_condattr_destroy
#define fpthrd_condattr_getpshared       fpthr_condattr_getpshared
#define fpthrd_condattr_setpshared       fpthr_condattr_setpshared
#define fpthrd_cond_init                 fpthr_cond_init
#define fpthrd_cond_destroy              fpthr_cond_destroy
#define fpthrd_cond_signal               fpthr_cond_signal
#define fpthrd_cond_broadcast            fpthr_cond_broadcast
#define fpthrd_cond_wait                 fpthr_cond_wait
#define fpthrd_cond_timedwait            fpthr_cond_timedwait
#define fpthrd_setconcurrency            fpthr_setconcurrency
#define fpthrd_getconcurrency            fpthrd_getconcurrency
#define fpthrd_strerror                  fpthrd_strerror
#define fpthrd_set_ftimespec             fpthrd_set_ftimespec
#define fpthrd_get_fsched_param          fpthrd_get_fsched_param
#define fpthrd_set_fsched_param          fpthrd_set_fsched_param
#define fpthrd_get_fsize                 fpthrd_get_fsize
#define fpthrd_set_fsize                 fpthrd_set_fsize
#endif

#ifdef APPEND_UNDERSCORE
#define fpthrd_attr_init                 fpthr_attr_init_
#define fpthrd_attr_destroy              fpthr_attr_destroy_
#define fpthrd_attr_setstacksize         fpthr_attr_setstacksize_
#define fpthrd_attr_getstacksize         fpthr_attr_getstacksize_
#define fpthrd_attr_setdetachstate       fpthr_attr_setdetachstate_
#define fpthrd_attr_getdetachstate       fpthr_attr_getdetachstate_
#define fpthrd_attr_setscope             fpthr_attr_setscope_
#define fpthrd_attr_getscope             fpthr_attr_getscope_
#define fpthrd_attr_setinheritsched      fpthr_attr_setinheritsched_
#define fpthrd_attr_getinheritsched      fpthr_attr_getinheritsched_
#define fpthrd_attr_setschedpolicy       fpthr_attr_setschedpolicy_
#define fpthrd_attr_getschedpolicy       fpthr_attr_getschedpolicy_
#define fpthrd_attr_setschedparam        fpthr_attr_setschedparam_
#define fpthrd_attr_getschedparam        fpthr_attr_getschedparam_
#define fpthrd_create                    fpthrd_create_
#define fpthrd_join                      fpthr_join_
#define fpthrd_exit                      fpthrd_exit_
#define fpthrd_detach                    fpthr_detach_
#define fpthr_self                       fpthr_self_
#define fpthrd_self                      fpthrd_self_
#define fpthr_equal                      fpthr_equal_
#define fpthrd_equal                     fpthrd_equal_
#define fpthrd_getschedparam             fpthr_getschedparam_
#define fpthrd_setschedparam             fpthr_setschedparam_
#define fpthrd_cancel                    fpthr_cancel_
#define fpthrd_setcancelstate            fpthr_setcancelstate_
#define fpthrd_setcanceltype             fpthr_setcanceltype_
#define fpthr_testcancel                 fpthr_testcancel_
#define fpthrd_testcancel                fpthrd_testcancel_
#define fpthrd_cleanup_push              fpthrd_cleanup_push_
#define fpthrd_cleanup_pop               fpthrd_cleanup_pop_
#define fpthrd_mutexattr_init            fpthr_mutexattr_init_
#define fpthrd_mutexattr_destroy         fpthr_mutexattr_destroy_
#define fpthrd_mutexattr_getpshared      fpthr_mutexattr_getpshared_
#define fpthrd_mutexattr_setpshared      fpthr_mutexattr_setpshared_
#define fpthrd_mutexattr_setprotocol     fpthr_mutexattr_setprotocol_
#define fpthrd_mutexattr_getprotocol     fpthr_mutexattr_getprotocol_
#define fpthrd_mutexattr_setprioceiling  fpthr_mutexattr_setprioceiling_
#define fpthrd_mutexattr_getprioceiling  fpthr_mutexattr_getprioceiling_
#define fpthrd_mutex_init                fpthr_mutex_init_
#define fpthrd_mutex_destroy             fpthr_mutex_destroy_
#define fpthrd_mutex_lock                fpthr_mutex_lock_
#define fpthrd_mutex_trylock             fpthr_mutex_trylock_
#define fpthrd_mutex_unlock              fpthr_mutex_unlock_
#define fpthrd_mutex_setprioceiling      fpthr_mutex_setprioceiling_
#define fpthrd_mutex_getprioceiling      fpthr_mutex_getprioceiling_
#define fpthrd_condattr_init             fpthr_condattr_init_
#define fpthrd_condattr_destroy          fpthr_condattr_destroy_
#define fpthrd_condattr_getpshared       fpthr_condattr_getpshared_
#define fpthrd_condattr_setpshared       fpthr_condattr_setpshared_
#define fpthrd_cond_init                 fpthr_cond_init_
#define fpthrd_cond_destroy              fpthr_cond_destroy_
#define fpthrd_cond_signal               fpthr_cond_signal_
#define fpthrd_cond_broadcast            fpthr_cond_broadcast_
#define fpthrd_cond_wait                 fpthr_cond_wait_
#define fpthrd_cond_timedwait            fpthr_cond_timedwait_
#define fpthrd_setconcurrency            fpthr_setconcurrency_
#define fpthrd_getconcurrency            fpthrd_getconcurrency_
#define fpthrd_strerror                  fpthrd_strerror_
#define fpthrd_set_ftimespec             fpthrd_set_ftimespec_
#define fpthrd_get_fsched_param          fpthrd_get_fsched_param_
#define fpthrd_set_fsched_param          fpthrd_set_fsched_param_
#define fpthrd_get_fsize                 fpthrd_get_fsize_
#define fpthrd_set_fsize                 fpthrd_set_fsize_
#endif

#ifdef APPEND_TWO_UNDERSCORES
#define fpthrd_attr_init                 fpthr_attr_init__
#define fpthrd_attr_destroy              fpthr_attr_destroy__
#define fpthrd_attr_setstacksize         fpthr_attr_setstacksize__
#define fpthrd_attr_getstacksize         fpthr_attr_getstacksize__
#define fpthrd_attr_setdetachstate       fpthr_attr_setdetachstate__
#define fpthrd_attr_getdetachstate       fpthr_attr_getdetachstate__
#define fpthrd_attr_setscope             fpthr_attr_setscope__
#define fpthrd_attr_getscope             fpthr_attr_getscope__
#define fpthrd_attr_setinheritsched      fpthr_attr_setinheritsched__
#define fpthrd_attr_getinheritsched      fpthr_attr_getinheritsched__
#define fpthrd_attr_setschedpolicy       fpthr_attr_setschedpolicy__
#define fpthrd_attr_getschedpolicy       fpthr_attr_getschedpolicy__
#define fpthrd_attr_setschedparam        fpthr_attr_setschedparam__
#define fpthrd_attr_getschedparam        fpthr_attr_getschedparam__
#define fpthrd_create                    fpthrd_create__
#define fpthrd_join                      fpthr_join__
#define fpthrd_exit                      fpthrd_exit__
#define fpthrd_detach                    fpthr_detach__
#define fpthr_self                       fpthr_self__
#define fpthrd_self                      fpthrd_self__
#define fpthr_equal                      fpthr_equal__
#define fpthrd_equal                     fpthrd_equal__
#define fpthrd_getschedparam             fpthr_getschedparam__
#define fpthrd_setschedparam             fpthr_setschedparam__
#define fpthrd_cancel                    fpthr_cancel__
#define fpthrd_setcancelstate            fpthr_setcancelstate__
#define fpthrd_setcanceltype             fpthr_setcanceltype__
#define fpthr_testcancel                 fpthr_testcancel__
#define fpthrd_testcancel                fpthrd_testcancel__
#define fpthrd_cleanup_push              fpthrd_cleanup_push__
#define fpthrd_cleanup_pop               fpthrd_cleanup_pop__
#define fpthrd_mutexattr_init            fpthr_mutexattr_init__
#define fpthrd_mutexattr_destroy         fpthr_mutexattr_destroy__
#define fpthrd_mutexattr_getpshared      fpthr_mutexattr_getpshared__
#define fpthrd_mutexattr_setpshared      fpthr_mutexattr_setpshared__
#define fpthrd_mutexattr_setprotocol     fpthr_mutexattr_setprotocol__
#define fpthrd_mutexattr_getprotocol     fpthr_mutexattr_getprotocol__
#define fpthrd_mutexattr_setprioceiling  fpthr_mutexattr_setprioceiling__
#define fpthrd_mutexattr_getprioceiling  fpthr_mutexattr_getprioceiling__
#define fpthrd_mutex_init                fpthr_mutex_init__
#define fpthrd_mutex_destroy             fpthr_mutex_destroy__
#define fpthrd_mutex_lock                fpthr_mutex_lock__
#define fpthrd_mutex_trylock             fpthr_mutex_trylock__
#define fpthrd_mutex_unlock              fpthr_mutex_unlock__
#define fpthrd_mutex_setprioceiling      fpthr_mutex_setprioceiling__
#define fpthrd_mutex_getprioceiling      fpthr_mutex_getprioceiling__
#define fpthrd_condattr_init             fpthr_condattr_init__
#define fpthrd_condattr_destroy          fpthr_condattr_destroy__
#define fpthrd_condattr_getpshared       fpthr_condattr_getpshared__
#define fpthrd_condattr_setpshared       fpthr_condattr_setpshared__
#define fpthrd_cond_init                 fpthr_cond_init__
#define fpthrd_cond_destroy              fpthr_cond_destroy__
#define fpthrd_cond_signal               fpthr_cond_signal__
#define fpthrd_cond_broadcast            fpthr_cond_broadcast__
#define fpthrd_cond_wait                 fpthr_cond_wait__
#define fpthrd_cond_timedwait            fpthr_cond_timedwait__
#define fpthrd_setconcurrency            fpthr_setconcurrency__
#define fpthrd_getconcurrency            fpthrd_getconcurrency__
#define fpthrd_strerror                  fpthrd_strerror__
#define fpthrd_set_ftimespec             fpthrd_set_ftimespec__
#define fpthrd_get_fsched_param          fpthrd_get_fsched_param__
#define fpthrd_set_fsched_param          fpthrd_set_fsched_param__
#define fpthrd_get_fsize                 fpthrd_get_fsize__
#define fpthrd_set_fsize                 fpthrd_set_fsize__
#endif

/*******************************************************************************
***********                                                          ***********
***********                                                          ***********
***********        Bindings to Pthreads library functions            ***********
***********                                                          ***********
***********                                                          ***********
*******************************************************************************/

/*----------------------------------------------------------------------------*\
 |                                                                            |
 |                       Dope passed from C to Fortran                        |
 |                                                                            |
\*----------------------------------------------------------------------------*/

  void fpthrd_set_ftimespec(PARAM_INT *change_sec, PARAM_INT *change_nanosec, struct timespec *waittime)
{
     struct timespec t;
     t.tv_sec=time(NULL)+ (int)*change_sec;
     t.tv_nsec=(int)*change_nanosec;
     *waittime=t;
}
  void fpthrd_set_fsched_param(PARAM_INT *schedule_value, struct sched_param *sched)
{
  struct sched_param Lsched;
  Lsched=*sched;
  Lsched.sched_priority=(int)*schedule_value;
  *sched=Lsched;
}
  void fpthrd_get_fsched_param(PARAM_INT *schedule_value, struct sched_param *sched)
{
    struct sched_param Lsched;
    Lsched=*sched;
   *schedule_value=(PARAM_INT) Lsched.sched_priority;
}
  void fpthrd_set_fsize(PARAM_INT *size_value, size_t *size)
{
  size_t Lsize;

  Lsize=(size_t) *size_value;
  *size=Lsize;
}
  void fpthrd_get_fsize(PARAM_INT *size_value, size_t *size)
{
  size_t Lsize;
  Lsize=*size;
  *size_value=(PARAM_INT)Lsize;
}

void fpthrd_strerror(PARAM_INT *errnum, int *MESSAGE, int* MESSAGE_SIZE)
{
 char *message; int i; int message_size;
 int status;

 status = (int) *errnum;
 message=strerror(status);
 if(message == NULL)
  {
    MESSAGE[0]=0; return;
  }
 /* The error message is moved as integers (ASCII sequence) to Fortran.
    Fortran will pack the integers into a character string. */
  message_size=*MESSAGE_SIZE;
  for(i=0;i < message_size;i++)
                              {
/* Copy until a zero byte is  passed. */
                              MESSAGE[i]=message[i];
      if(message[i] == 0)     break;
                              }
}
/*----------------------------------------------------------------------------*\
 |                                                                            |
 |                       Thread creation attributes                           |
 |                                                                            |
\*----------------------------------------------------------------------------*/
void fpthrd_attr_init(pthread_attr_t *attr, PARAM_INT *ierr)
{
      *ierr = (PARAM_INT)pthread_attr_init(attr);
       return;
}

void fpthrd_attr_destroy(pthread_attr_t *attr, PARAM_INT *ierr)
{
      *ierr = (PARAM_INT)pthread_attr_destroy(attr);
       return;
}

void fpthrd_attr_setstacksize(pthread_attr_t *attr, size_t *stacksize, PARAM_INT *ierr)
{
   size_t temp;
   temp= *stacksize;
   *ierr = (PARAM_INT)pthread_attr_setstacksize(attr, temp);
   return;
}

 void fpthrd_attr_getstacksize(const pthread_attr_t *attr, size_t *stacksize, PARAM_INT *ierr)
{
   size_t temp;
   *ierr = (PARAM_INT)pthread_attr_getstacksize(attr, &temp);
   *stacksize= temp;
   return;
}


void fpthrd_attr_setdetachstate(pthread_attr_t *attr, PARAM_INT *detachstate, PARAM_INT *ierr)
{
   int ldetachstate;
   ldetachstate = (int) *detachstate;
   *ierr = (PARAM_INT)pthread_attr_setdetachstate(attr, ldetachstate);
    return;
}

void fpthrd_attr_getdetachstate(const pthread_attr_t *attr, PARAM_INT *detachstate, PARAM_INT *ierr)
{
   int ldetachstate;

   *ierr = (PARAM_INT)pthread_attr_getdetachstate(attr, &ldetachstate);
    *detachstate = (PARAM_INT) ldetachstate;
    return;
}



/*----------------------------------------------------------------------------*\
 |                                                                            |
 |                       Thread scheduling attributes                         |
 |                                                                            |
\*----------------------------------------------------------------------------*/

void fpthrd_attr_setscope(pthread_attr_t *attr, PARAM_INT *scope, PARAM_INT *ierr)
{
   int lscope;
   lscope = (int) *scope;
   *ierr = (PARAM_INT)pthread_attr_setscope(attr, lscope);
    return;
}

void fpthrd_attr_getscope(const pthread_attr_t *attr, PARAM_INT *scope, PARAM_INT *ierr)
{
      int lscope;
   *ierr = (PARAM_INT)pthread_attr_getscope(attr, &lscope);
   *scope = (PARAM_INT)lscope;
    return;
}

void fpthrd_attr_setinheritsched(pthread_attr_t *attr, PARAM_INT *inherit, PARAM_INT *ierr)
{
      int linherit;
   linherit = (int) *inherit;
   *ierr = (PARAM_INT)pthread_attr_setinheritsched(attr, linherit);
    return;
}

void fpthrd_attr_getinheritsched(pthread_attr_t *attr, PARAM_INT *inheritsched, PARAM_INT *ierr)
{
      int linherit;
   *ierr = (PARAM_INT)pthread_attr_getinheritsched(attr, &linherit);
      *inheritsched = (PARAM_INT)linherit;
   return;
}

void fpthrd_attr_setschedpolicy(pthread_attr_t *attr, PARAM_INT *policy, PARAM_INT *ierr)
{
      int lpolicy;
   lpolicy = (int) *policy;
   *ierr = (PARAM_INT)pthread_attr_setschedpolicy(attr, lpolicy);
    return;
}

void fpthrd_attr_getschedpolicy(pthread_attr_t *attr, PARAM_INT *policy, PARAM_INT *ierr)
{
      int lpolicy;
   *ierr = (PARAM_INT)pthread_attr_getschedpolicy(attr, &lpolicy);
   *policy = (PARAM_INT)lpolicy;
    return;
}

void fpthrd_attr_setschedparam(pthread_attr_t *attr, const struct sched_param *param, PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_attr_setschedparam(attr, param);
    return;
}

void fpthrd_attr_getschedparam(pthread_attr_t *attr, struct sched_param *param, PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_attr_getschedparam(attr, param);
    return;
}



/*----------------------------------------------------------------------------*\
 |                                                                            |
 |                       Thread creation and control                          |
 |                                                                            |
\*----------------------------------------------------------------------------*/

void fpthrd_create(pthread_t *thread_id, pthread_attr_t *attr,
                     void *(*start_routine)(void *), void *arg, PARAM_INT *ierr)
{
  pthread_t *tid;
  pthread_attr_t *lattr;
  PARAM_INT arg_in;


/*  This argument should never be NULL.*/
/*  Butenhof uses NULL in one of the examples in his book.      */

/*
  if(*(PARAM_INT *) thread_id == FORTRAN_NULL)
*/
  if(*(int *)thread_id == FORTRAN_NULL)
    {
     tid=NULL;
printf("Thread id is NULL\n");
    }
  else
    {
     tid=thread_id;
printf("Thread id is not empty\n");
    }

  tid=thread_id;

  if(*(int *)attr == FORTRAN_NULL)
    {
      lattr=NULL;
printf("Thread attributes is NULL\n");
    }
  else
  {
   lattr=attr;
printf("Thread attributes is non-empty\n");
  }
arg_in = *(PARAM_INT *)arg;
if(arg_in == FORTRAN_NULL)
{
printf("Function argument is NULL\n");
  *ierr = (PARAM_INT)pthread_create(tid, lattr, start_routine,  NULL);
}
else
{
printf("Function argument is not empty\n");
  *ierr = (PARAM_INT)pthread_create(tid, lattr, start_routine, arg);
}
  return;
}


void fpthrd_join(pthread_t *thread_id, PARAM_INT *exitcode, PARAM_INT *ierr)
{
   void *exit;
   exit =exitcode;
   if(*exitcode == (PARAM_INT)FORTRAN_NULL)
/*
   if(*(int *)exitcode == FORTRAN_NULL)
*/
   {
     exit=NULL;
   }
   *ierr = (PARAM_INT)pthread_join(*thread_id, &exit);

/* #if SGI == I8 */
#ifdef SGI
/*
 *  If exitcode is long long, examine high bytes for cancel code.
 */
      if ( exit == PTHREAD_CANCELED) {
      *exitcode=(PARAM_INT)FPTHREAD_CANCELED;
      }
/*#endif*/

/* Use a reserved value to signal Fortran that a thread was cancelled. */
/*#ifdef SUN*/
#else
   if(exit == PTHREAD_CANCELED) *exitcode=FPTHREAD_CANCELED;
   else
   {
     if(*exitcode == FORTRAN_NULL) return;
     *exitcode = (INT_CAST *)exit;
   }
#endif

}


void fpthrd_exit(PARAM_INT *val)
{
   void *lval;

   lval = (int *)(*val);   /* exit codes restricted to integer type */
   pthread_exit(lval);
}

void fpthrd_detach(pthread_t *thread_id, PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_detach(*thread_id);
}


void fpthrd_self(pthread_t *self)
{
   *self = pthread_self();   /* always succeeds, no error code */
}

void fpthr_self(pthread_t *self)
{
  fpthrd_self(self);
}

void fpthrd_equal(pthread_t *thread_id1, pthread_t *thread_id2, PARAM_INT *flag)
{
     *flag = (PARAM_INT)pthread_equal(*thread_id1, *thread_id2);   /* non-zero if equal */
}

void fpthr_equal(pthread_t *thread_id1, pthread_t *thread_id2, int *flag)
{
  fpthrd_equal(thread_id1, thread_id2, flag);
}

/*
void fpthrd_once(pthread_once_t *once_block, void (*init_routine)(void), PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_once(once_block, init_routine);
}
*/



/*----------------------------------------------------------------------------*\
 |                                                                            |
 |                       Thread scheduling control                            |
 |                                                                            |
\*----------------------------------------------------------------------------*/

void fpthrd_getschedparam(pthread_t *thread_id, PARAM_INT *policy,
struct sched_param *param, PARAM_INT *ierr)
{
      int lpolicy;
   *ierr = (PARAM_INT)pthread_getschedparam(*thread_id, &lpolicy, param);
      *policy = (PARAM_INT)lpolicy;
   return;
}

void fpthrd_setschedparam(pthread_t *thread_id, int *policy,
                            const struct sched_param *param, int *ierr)
{
      int lpolicy;
      lpolicy = *(int *)policy;
   *ierr = pthread_setschedparam(*thread_id, lpolicy, param);
    return;
}



/*----------------------------------------------------------------------------*\
 |                                                                            |
 |                          Thread cancellation                               |
 |                                                                            |
\*----------------------------------------------------------------------------*/

void fpthrd_cancel(pthread_t *thread_id, PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_cancel(*thread_id);
}


void fpthrd_setcancelstate(PARAM_INT *state, PARAM_INT *oldstate, PARAM_INT *ierr)
{
   int lstate, loldstate;
   lstate = (int) *state;
   *ierr = (PARAM_INT)pthread_setcancelstate(lstate, &loldstate);
   *oldstate = (PARAM_INT)loldstate;
}

void fpthrd_setcanceltype(PARAM_INT *type, PARAM_INT *oldtype, PARAM_INT *ierr)
{
   int ltype, loldtype;
   ltype = (int) *type;
   *ierr = (PARAM_INT)pthread_setcanceltype(ltype, &loldtype);
   *oldtype = (PARAM_INT)loldtype;
}


void fpthrd_testcancel()
{
   pthread_testcancel();   /* always succeeds; no error code */
}

void fpthr_testcancel()
{
   pthread_testcancel();   /* always succeeds; no error code */
}





/*----------------------------------------------------------------------------*\
 |                                                                            |
 |                   Mutex attribute manipulation and initialization          |
 |                                                                            |
\*----------------------------------------------------------------------------*/

void fpthrd_mutexattr_init(pthread_mutexattr_t *attr, PARAM_INT *ierr)
{
     *ierr = (PARAM_INT)pthread_mutexattr_init(attr);
      return;
}


void fpthrd_mutexattr_destroy(pthread_mutexattr_t *attr, PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_mutexattr_destroy(attr);
    return;
}
void fpthrd_mutexattr_getpshared(pthread_mutexattr_t *attr, PARAM_INT *pshared, PARAM_INT *ierr)
#ifdef LINUX
{
        *ierr=ENOTSUP;
        return;
}
#else
{
      int lpshared;
     *ierr = (PARAM_INT)pthread_mutexattr_getpshared(attr, &lpshared);
      *pshared = (PARAM_INT)lpshared;
      return;
}
#endif

void fpthrd_mutexattr_setpshared(pthread_mutexattr_t *attr, PARAM_INT *pshared, PARAM_INT *ierr)
#ifdef LINUX
{
        *ierr=ENOTSUP;
        return;
}
#else
{
   int lpshared;
   lpshared = (int) *pshared;
   *ierr = (PARAM_INT)pthread_mutexattr_setpshared(attr, lpshared);
    return;
}
#endif

void fpthrd_mutex_getprioceiling(pthread_mutex_t *mutex, PARAM_INT *ceiling, PARAM_INT *ierr)
#ifdef LINUX
{
        *ierr=ENOTSUP;
        return;
}
#elif CPQ
{
        *ierr=ENOTSUP;
        return;
}
#else
{
      int lceiling;
     *ierr = (PARAM_INT)pthread_mutex_getprioceiling(mutex, &lceiling);
      *ceiling = (PARAM_INT)lceiling;
      return;
}
#endif

void fpthrd_mutex_setprioceiling(pthread_mutex_t *mutex, PARAM_INT *ceiling,
                                 PARAM_INT *old_ceiling, PARAM_INT *ierr)
#ifdef LINUX
{
        *ierr=ENOTSUP;
        return;
}
#elif CPQ
{
        *ierr=ENOTSUP;
        return;
}
#else
{
      int lceiling, loldceiling;
   lceiling = (int) *ceiling;
   *ierr = (PARAM_INT)pthread_mutex_setprioceiling(mutex, lceiling, &loldceiling);
   *old_ceiling = (PARAM_INT)loldceiling;
    return;
}
#endif

/*----------------------------------------------------------------------------*\
 |                                                                            |
 |                     Mutex scheduling attributes                            |
 |                                                                            |
\*----------------------------------------------------------------------------*/
void fpthrd_mutexattr_setprotocol(pthread_mutexattr_t *attr,
                                    PARAM_INT *protocol, PARAM_INT *ierr)
#ifdef LINUX
{
        *ierr=ENOTSUP;
        return;
}
#elif CPQ
{
        *ierr=ENOTSUP;
        return;
}
#else
{
   int lprotocol;
   lprotocol = (int) *protocol;
   *ierr = (PARAM_INT)pthread_mutexattr_setprotocol(attr, lprotocol);
    return;
}
#endif


void fpthrd_mutexattr_getprotocol(pthread_mutexattr_t *attr,
                                    PARAM_INT *protocol, PARAM_INT *ierr)
#ifdef LINUX
{
        *ierr=ENOTSUP;
        return;
}
#elif CPQ
{
        *ierr=ENOTSUP;
        return;
}
#else
{
   int lprotocol;
   *ierr = (PARAM_INT)pthread_mutexattr_getprotocol(attr, &lprotocol);
   *protocol = (PARAM_INT)lprotocol;
    return;
}
#endif
/*----------------------------------------------------------------------------*\
 |                                                                            |
 |                     Mutex creation and control                             |
 |                                                                            |
\*----------------------------------------------------------------------------*/

void fpthrd_mutex_init(pthread_mutex_t *mutex, pthread_mutexattr_t *attr, PARAM_INT *ierr)
{
   pthread_mutexattr_t *lattr;


      if(*(PARAM_INT *) attr == FORTRAN_NULL)
      {
        lattr=NULL;
      }
      else
      {
        lattr = attr;
      }

      *ierr = (PARAM_INT)pthread_mutex_init(mutex, lattr);
       return;
}


void fpthrd_mutex_destroy(pthread_mutex_t *mutex, PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_mutex_destroy(mutex);
    return;
}


void fpthrd_mutex_lock(pthread_mutex_t *mutex, PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_mutex_lock(mutex);
    return;
}


void fpthrd_mutex_trylock(pthread_mutex_t *mutex, PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_mutex_trylock(mutex);
    return;
}


void fpthrd_mutex_unlock(pthread_mutex_t *mutex, PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_mutex_unlock(mutex);
    return;
}


/*----------------------------------------------------------------------------*\
 |                                                                            |
 |                       Mutex scheduling control                             |
 |                                                                            |
\*----------------------------------------------------------------------------*/
void fpthrd_mutexattr_setprioceiling(pthread_mutexattr_t *attr,
                                       PARAM_INT *prioceiling, PARAM_INT *ierr)
#ifdef LINUX
{
        *ierr=ENOTSUP;
        return;
}
#elif CPQ
{
        *ierr=ENOTSUP;
        return;
}
#else
{
      int lprioceiling;
   lprioceiling = (int) *prioceiling;
   *ierr = (PARAM_INT)pthread_mutexattr_setprioceiling(attr, lprioceiling);
    return;
}
#endif

void fpthrd_mutexattr_getprioceiling(pthread_mutexattr_t *attr,
                                       PARAM_INT *prioceiling, PARAM_INT *ierr)
#ifdef LINUX
{
        *ierr=ENOTSUP;
        return;
}
#elif CPQ
{
        *ierr=ENOTSUP;
        return;
}
#else
{
      int lprioceiling;
   *ierr = (PARAM_INT)pthread_mutexattr_getprioceiling(attr, &lprioceiling);
   *prioceiling = (PARAM_INT)lprioceiling;
    return;
}
#endif

/*----------------------------------------------------------------------------*\
 |                                                                            |
 |            Condition variable initialization attributes                    |
 |                                                                            |
\*----------------------------------------------------------------------------*/

void fpthrd_condattr_init(pthread_condattr_t *attr, PARAM_INT *ierr)
{
     *ierr = (PARAM_INT)pthread_condattr_init(attr);
      return;
}


void fpthrd_condattr_destroy(pthread_condattr_t *attr, PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_condattr_destroy(attr);
    return;
}
void fpthrd_condattr_getpshared(pthread_condattr_t *attr, PARAM_INT *pshared, PARAM_INT *ierr)
#ifndef LINUX
{
   int lpshared;
   *ierr = (PARAM_INT)pthread_condattr_getpshared(attr, &lpshared);
   *pshared = (PARAM_INT)lpshared;
    return;
}
#else
{
        *ierr=ENOTSUP;
        return;
}
#endif

void fpthrd_condattr_setpshared(pthread_condattr_t *attr, PARAM_INT *pshared, PARAM_INT *ierr)
#ifndef LINUX
{
   int lpshared;
   lpshared = (int) *pshared;
   *ierr = (PARAM_INT)pthread_condattr_setpshared(attr, lpshared);
    return;
}
#else
{
        *ierr=ENOTSUP;
        return;
}
#endif

/*----------------------------------------------------------------------------*\
 |                                                                            |
 |                     Condition variable operations                          |
 |                                                                            |
\*----------------------------------------------------------------------------*/

void fpthrd_cond_init(pthread_cond_t *cond, pthread_condattr_t *attr, PARAM_INT *ierr)
{
   pthread_condattr_t *lattr;
      if(*(PARAM_INT *)attr == FORTRAN_NULL)
      {
        lattr=NULL;
      }
      else
      {
        lattr = attr;
      }

      *ierr = (PARAM_INT)pthread_cond_init(cond, lattr);
       return;
}

void fpthrd_cond_destroy(pthread_cond_t *cond, PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_cond_destroy(cond);
   return;
}


void fpthrd_cond_signal(pthread_cond_t *cond, PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_cond_signal(cond);
    return;
}


void fpthrd_cond_broadcast(pthread_cond_t *cond, PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_cond_broadcast(cond);
    return;
}


void fpthrd_cond_wait(pthread_cond_t *cond, pthread_mutex_t *mutex, PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_cond_wait(cond, mutex);
    return;
}

void fpthrd_cond_timedwait(pthread_cond_t *cond, pthread_mutex_t *mutex,
                             const struct timespec *abstime, PARAM_INT *ierr)
{
   *ierr = (PARAM_INT)pthread_cond_timedwait(cond, mutex, abstime);
    return;
}


/*----------------------------------------------------------------------------*\
 |                                                                            |
 |   Bindings to nonportable functions                                        |
 |                                                                            |
\*----------------------------------------------------------------------------*/
  void fpthrd_setconcurrency(PARAM_INT *level, PARAM_INT *ierr)
   {
#ifdef SGI
      int llevel;

      llevel = (int) *level;
      *ierr = (PARAM_INT)pthread_setconcurrency(llevel);
#endif
#ifdef LINUX
      int llevel;

      llevel = (int) *level;
      *ierr = (PARAM_INT)pthread_setconcurrency(llevel);
#endif

#ifdef          SUN
      int       llevel;

      llevel =  *(int *)level;
      if(llevel > 0) {thr_setconcurrency(llevel); *ierr = 0; return;}
                *ierr=EINVAL;
#endif

#ifdef          IBM
                *ierr=ENOTSUP;
#endif
   }

void                             fpthrd_getconcurrency(PARAM_INT *level)
                                 {
#ifdef                           SGI
      *level =                   (PARAM_INT)pthread_getconcurrency();
#endif
#ifdef                           LINUX
      *level =                   (PARAM_INT)pthread_getconcurrency();
#endif
#ifdef                           SUN
      *level =                   (PARAM_INT)thr_getconcurrency();
#endif

#ifdef                           IBM
      *level =                   0;
#endif
  }

