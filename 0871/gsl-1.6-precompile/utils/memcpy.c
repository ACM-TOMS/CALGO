#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* Copy LEN bytes starting at SRCADDR to DESTADDR.  Result undefined
   if the source overlaps with the destination.
   Return DESTADDR. */

#if HAVE_CONFIG_H
#include <config.h>
#endif

char *
memcpy (destaddr, srcaddr, len)
     char *destaddr;
     const char *srcaddr;
     int  len;
{
  char *dest = destaddr;

  while (len-- > 0)
    *destaddr++ = *srcaddr++;
  return dest;
}
