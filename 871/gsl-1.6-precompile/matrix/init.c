#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

/*#define BASE_GSL_COMPLEX_LONG
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT
*/
#define BASE_LONG_DOUBLE
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef  BASE_FLOAT


/*
//this is not needed for init_source but we do it because we also need to disable
//these types in the other files like matrix.c and matrix_source.c
#define BASE_ULONG
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
#include "init_source.c"
#include "templates_off.h"
#undef  BASE_CHAR
*/

