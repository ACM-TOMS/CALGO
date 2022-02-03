#include "pampac.h"
/**********************************************************************/
/* This routine prints debug messages for the user's output benefit.  */
/* The verbose field of the options data structure restricts output   */
/* and the calling function provides its name in the string fname.    */
/* Variable amounts of data are passed after the format string.       */
/**********************************************************************/
void debug_print (int thresh, options_struct *opts,
                  const char *fname, const char *format, ...) {
    va_list args;

    if (opts->verbose < thresh)
       return;

    printf ("%s: ", fname);
    va_start (args, format);
    vprintf(format, args);
    fflush (stdout);
    va_end(args);
    return;
}
