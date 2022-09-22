#ifndef CALS_DEFINITIONS_H
#define CALS_DEFINITIONS_H

#ifndef NDEBUG
#define DEBUG(exp) exp
#else
#define DEBUG(exp) do {} while(0);
#endif

#if WITH_TIME
#define TIME(exp) exp
#else
#define TIME(exp) ;
#endif

typedef size_t dim_t;

#endif // CALS_DEFINITIONS_H
