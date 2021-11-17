#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#ifdef __GNUC__
#elif defined(__xlC__)
#define CLOCKTICKS CLK_TCK;
#else
#endif
