#we're building bertini, version 1.5.1, and the corresponding email is jon's
AC_INIT([bertini], [1.5.1], [hauenstein@nd.edu])

#
AC_CONFIG_AUX_DIR([config])

# turn on the keeping of produced objects in their folders.  this is for non-recursive make
# and autotools
# see Recursive Make Considered Harmful, and any number of demos.
AM_INIT_AUTOMAKE([subdir-objects])

#the language for Bertini is C
AC_LANG([C])


#another directive to use the m4 folder
AC_CONFIG_MACRO_DIR([m4])

#the only produced file will be a single Makefile.
AC_CONFIG_FILES([Makefile])

#find the CC compiler
AC_PROG_CC

#set up for building libraries
AC_PROG_LIBTOOL

#find flex / lex
AC_PROG_LEX

#find bison / yacc
AC_PROG_YACC

#find the linker
AC_PROG_LN_S

#find the command for making directories with their parents
AC_PROG_MKDIR_P

#enable the creation of shared libraries
AC_ENABLE_SHARED

#enable the creation of static libraries
AC_ENABLE_STATIC

#fire up libtool
LT_INIT

# the form of the following commands --
# AC_SEARCH_LIBS(function, libraries-list, action-if-found, action-if-not-found, extra-libraries)



AC_ARG_WITH([gmp],
            [AS_HELP_STRING([--without-gmp],
              [disable searching support for gmp])],
            [],
            [with_gmp=yes])

          AS_IF([test "x$with_gmp" != xno],
            AC_SEARCH_LIBS([__gmpz_init],[gmp], [],[
				AC_MSG_ERROR([unable to find gmp])
				]),[]
 				)


AC_ARG_WITH([mpfr],
            [AS_HELP_STRING([--without-mpfr],
              [disable searching support for mpfr])],
            [],
            [with_mpfr=yes])

          AS_IF([test "x$with_mpfr" != xno],
            AC_SEARCH_LIBS([mpfr_get_version],[mpfr], [],[
				AC_MSG_ERROR([unable to find mpfr])
				]),[]
 				)


AC_ARG_WITH([m],
            [AS_HELP_STRING([--without-m],
              [disable searching support for m])],
            [],
            [with_m=yes])

          AS_IF([test "x$with_m" != xno],
            AC_SEARCH_LIBS([cos],[m], [],[
				AC_MSG_ERROR([unable to find math library with cos()])
				]),[]
 				)




#AC_SEARCH_LIBS([mpfr_get_version],[mpfr], [],[
#				AC_MSG_ERROR([unable to find mpfr])
#				])






#this calls a file in the m4/ directory, which sets up the MPI wrapper stuffs
LX_FIND_MPI

AM_CONDITIONAL(BUILD_PARALLEL,[test "$have_C_mpi" = "yes"])


#wrap it up.
AC_OUTPUT



