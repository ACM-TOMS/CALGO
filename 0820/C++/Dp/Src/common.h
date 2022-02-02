//***************************************************************
// common.h
//****************************************************************

#ifndef COMMON_H
#define COMMON_H

#include <stdlib.h>

typedef double real;
typedef long integer;

#define maximum(x,y)   ((x)>(y) ? (x):(y))
#define minimum(x,y)   ((x)<(y) ? (x):(y))
#define absval(x)  ((x)<0 ? (-(x)):(x))

integer Log2( integer r );                 // log base 2 of integer r
bool is_pow_of_2( integer r );             // check if r is power of 2

//integer NumOfData(const char *filename);   
     // count total=number of reals in file, return total
#endif
