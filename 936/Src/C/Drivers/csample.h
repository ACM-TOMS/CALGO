//Time stamp: December 9, 2013, 11:25
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include "cmessy.h"
#define fp CTYPE_ // Defined as float, double, or long double
enum sample_enums{
  setup_sample=0, partial_message, finish_message
};
struct SAMPLE_TY_{
  struct CMESSY_TY_ se;
  int what;
};
void SAMPLE_(struct SAMPLE_TY_ *s);
