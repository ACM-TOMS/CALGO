#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* system.h: System-dependent declarations.  Include this first.
   $Id: system.h,v 1.2 2003/06/17 21:06:10 bjg Exp $

   Copyright (C) 1997 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  */

#ifndef TEXINFO_SYSTEM_H
#define TEXINFO_SYSTEM_H

#define _GNU_SOURCE

#include <config.h>

/* <unistd.h> should be included before any preprocessor test
   of _POSIX_VERSION.  */
#if HAVE_UNISTD_H
#include <unistd.h>
#endif /* HAVE_UNISTD_H */

#include <stdio.h>
#include <sys/types.h>
#include <ctype.h>

#if HAVE_LOCALE_H
#include <locale.h>
#endif
#include <libintl.h>

/* Don't use bcopy!  Use memmove if source and destination may overlap,
   memcpy otherwise.  */
#if HAVE_STRING_H
# if !STDC_HEADERS && HAVE_MEMORY_H
#  include <memory.h>
# endif
# include <string.h>
#else
# include <strings.h>
char *memchr ();
#endif

#ifdef STDC_HEADERS
#define getopt system_getopt
#include <stdlib.h>
#undef getopt
#else
extern char *getenv ();
#endif

#ifndef HAVE_STRERROR
extern char *strerror ();
#endif

#include <errno.h>
#ifndef errno
extern int  errno;
#endif
#ifdef VMS
#include <perror.h>
#endif

#include <sys/stat.h>

#if HAVE_SYS_FILE_H
#include <sys/file.h>
#endif /* HAVE_SYS_FILE_H */

#ifndef O_RDONLY
/* Since <fcntl.h> is POSIX, prefer that to <sys/fcntl.h>.
   This also avoids some useless warnings on (at least) Linux.  */
#if HAVE_FCNTL_H
#include <fcntl.h>
#else /* not HAVE_FCNTL_H */
#if HAVE_SYS_FCNTL_H
#include <sys/fcntl.h>
#endif /* not HAVE_SYS_FCNTL_H */
#endif /* not HAVE_FCNTL_H */
#endif /* not O_RDONLY */

#if HAVE_PWD_H
#include <pwd.h>
#endif
/* Some systems don't declare this function in pwd.h. */
struct passwd *getpwnam ();

/* Our library routines not included in any system library.  */
extern void *xmalloc (), *xrealloc ();
extern char *xstrdup ();

#endif /* TEXINFO_SYSTEM_H */
