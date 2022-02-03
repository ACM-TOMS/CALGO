#ifndef _UTIL_H
#define _UTIL_H

#include "common.h"

/* Tracing with EXTRAE */
#ifdef SEQTRACE
# include <mpitrace_user_events.h>
# define TRACE_INIT() SEQtrace_init();
# define TRACE_FINI() SEQtrace_fini();
# define TRACE_START() SEQtrace_eventandcounters( 1000, 1 );
# define TRACE_STOP() SEQtrace_eventandcounters( 1000, 0 );

/* Tracing with PAPI */
#elif PAPITRACE
# define TRACE_INIT()                                                           \
  int i, code, xEventSet, nEvents;                                              \
  char *env, *event[16];                                                        \
  char envCopy[256];                                                            \
  long long PAPI_Counters[256];                                                 \
                                                                                \
  xEventSet = PAPI_NULL;                                                        \
  if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) handle_error(1); \
  if (PAPI_create_eventset(&xEventSet) != PAPI_OK) handle_error(1);             \
                                                                                \
  /* Add events to our EventSet */                                              \
  nEvents = 0;                                                                  \
  if (env = getenv("PAPI_COUNTERS")) {                                          \
    strcpy((char*)&envCopy, env);                                               \
    if (event[nEvents] = strtok((char*)&envCopy, ",")) {                        \
      PAPI_event_name_to_code(event[nEvents], &code);                           \
      PAPI_Counters[nEvents] = 0;                                               \
      if (PAPI_add_event(xEventSet, code) != PAPI_OK) handle_error(1);          \
      nEvents++;                                                                \
    }                                                                           \
                                                                                \
    while (event[nEvents] = strtok(NULL, ",")) {                                \
      PAPI_event_name_to_code(event[nEvents], &code);                           \
      PAPI_Counters[nEvents] = 0;                                               \
      if (PAPI_add_event(xEventSet, code) != PAPI_OK) handle_error(1);          \
      nEvents++;                                                                \
    }                                                                           \
  }                                                                             \
  else {                                                                        \
    /* Add default counter if env. variable not set - PAPI_TOT_INS */           \
    PAPI_Counters[nEvents] = 0;                                                 \
    if (PAPI_add_event(xEventSet, PAPI_TOT_INS) != PAPI_OK) handle_error(1);    \
    strcpy(envCopy, "PAPI_TOT_INS");                                            \
    event[0] = &envCopy[0];                                                     \
    nEvents++;                                                                  \
  }

# define TRACE_FINI()
# define TRACE_START()                                                          \
  /* Start counting */                                                          \
  if (PAPI_start(xEventSet) != PAPI_OK) handle_error(1);

# define TRACE_STOP()                                                           \
  /* Stop counting */                                                           \
  if (PAPI_stop(xEventSet, PAPI_Counters) != PAPI_OK) handle_error(1);          \
                                                                                \
  /* Print HWC results */                                                       \
  for (i= 0; i< nEvents; i++) {                                                 \
    printf("%lld\t%s\n", PAPI_Counters[i], event[i]);                           \
  }


/* Tracing with likwid */
#elif LIKWID
# include <likwid.h>
# define TRACE_INIT()  likwid_markerInit();
# define TRACE_FINI()  likwid_markerClose();
# define TRACE_START() likwid_markerStartRegion("benchmark");
# define TRACE_STOP()  likwid_markerStopRegion("benchmark");

/* Tracing with PAPIEX */
#else
# define TRACE_INIT()
# define TRACE_FINI()
# define TRACE_START()
# define TRACE_STOP()
#endif


/* Maximum tolerance (double precision) for relative error check */
#define TOLERANCE 1.0E-15


/* Error messages for PAPI calls */
static const char *ErrMsg[] = {
  "Example program",
  "Couldn't find any hardware counters",
  "Couldn't start counting events",
  "Couldn't stop counting events",
  "Couldn't translate event code into event name",
  "Couldn't initialize overflow handler",
  "Couldn't reset event set"
};


/*
  This initializes the array A to be all 1's.  
  This is nearly superfluous (could use memset), but
  provides convenience and consistency nonetheless...
 */
void StencilInit(int nx, int ny, int nz, /* size of the array */
                 double *A);             /* the array to initialize to 1's */

/*
  This function determines ticks per second.
  Inspired by OSKI function (bebop.cs.berkeley.edu)
*/
double seconds_per_tick();

/*
  Function to clear the cache, preventing data items in cache
  from making subsequent trials run faster.
*/
double clear_cache();

/*
  Handles error for PAPI routines
*/
void handle_error(int errcode);

/*
  Compare results between two executions using TOLERANCE as relative error
*/
void check_vals(double* A, double* B, int nx, int ny, int nz);

/*
  Standard prototype for any StencilProbe algorithm
*/
void StencilProbe(double* A0, double* Anext, int nx, int ny, int nz,
                  int tx, int ty, int tz, int timesteps, int length);

#endif
