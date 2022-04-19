/*********************************************************************/
/*                                                                   */
/* Copyright (C) 2006  Andrey Nikolayevich Chernikov                 */
/*                                                                   */
/* defs.h                                                            */
/*                                                                   */
/* `#define' statements for the PCDM                                 */
/*                                                                   */
/*********************************************************************/

#ifndef __DEFS_H__
#define __DEFS_H__

// ** TYPES **
#define UINT        unsigned int
#define REGION      Region
#define SIZE_INT    (sizeof (int))
#define BITS_INT    (8 * sizeof (int))

// ** TIMING **
#define TIME_DIFF(t2, t1)                                                      \
    (t2.tv_sec - t1.tv_sec) + 1e-6 * (t2.tv_usec - t1.tv_usec) 

#define BEGIN_TIMING(i) { gettimeofday (&(timer1[i]), 0); }

#define END_TIMING(i)   {                                                      \
                            gettimeofday (&(timer2), 0);                       \
                            timer0[i] += TIME_DIFF (timer2, timer1[i]);        \
                        }
                        
#define TIME_TOTAL               0
#define TIME_INIT                1
#define TIME_METIS               2
#define TIME_CDT                 3
#define TIME_REFINE_PCDM         4
#define TIME_IO_READ             5
#define TIME_IO_WRITE            6
#define TIME_COMM                7
#define TIME_SYNC                8
#define TIME_MPI_PROBE           9
#define TIME_FREE               10
#define TIME_STATS              11
#define NUM_TIMES               12
                        
// ** MEMORY MANAGEMENT **
#define MALLOC(p, t, n) { p = (t) malloc (n);     assert((p != 0)); }
#define CALLOC(p, t, n) { p = (t) calloc (1, n);  assert((p != 0)); }
#define REALLOC(p, t, n){ p = (t) realloc (p, n); assert((p != 0)); }
#define FREE(p)         { if (p != 0) free (p); p = 0; }
#define NEW(p, smth)    { p = new smth; assert((p != 0)); }

// ** STL WRAPPERS **
#define PUSH_BACK(A, x)  { (A).push_back(x); }
#define PUSH_FRONT(A, x) { (A).push_front(x); }
#define FOR(i, A)        for (i = (A).begin(); i != (A).end(); i++) 
#define EXISTS(S, x)     ((S).find(x) != (S).end())
                        
// ** ARITHMETIC **
#define MIN(a, b)       (a < b ? a : b)
#define MAX(a, b)       (a > b ? a : b)
#define ABS(x)          (x > 0 ? x : -x)
#define ONE             ((int) 1)
#define PI              3.141592653589793238462643383279502884197169399375105820974944592308
#define DEGREES(r)      (r * 180.0 / PI)
#define RADIANS(d)      (d * PI / 180.0)

// ** MESSAGE TYPES **
#define MSG_REQ_WORK     0
#define MSG_OBJ_DATA     1
#define MSG_OBJ_DONE     2
#define MSG_TERMINATE    3
#define MSG_REQ_INFO     4
#define MSG_OBJ_INFO     5
#define MSG_NEED_LB      6
#define MSG_DENY_LB      7
#define MSG_STAT         8
#define MSG_SPLIT        9
#define MSG_FINISHED    10 
#define MSG_UPDATE      11
#define MSG_POINT_NUM   12

#define MSG_HEADER_INTS  4
#define MSG_HEADER_SIZE (MSG_HEADER_INTS * sizeof(int))


// ** MESHING **
#define MAX_CAVITY_SIZE 512
#define Q_COUNT         (int)32
#define LOG_Q_COUNT     (int) 5

#endif
