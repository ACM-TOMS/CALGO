




#we're building bertini_real, version 1.4.0, and the corresponding email is dani brake's
AC_INIT([bertini_real], [1.4.0], [danielthebrake@gmail.com],[bertini_real], [http://www.danielthebrake.org/])


# Force autoconf to be at least this version number:
AC_PREREQ([2.68])


#
AC_CONFIG_AUX_DIR([config])

# turn on the keeping of produced objects in their folders.  this is for non-recursive make
# and autotools
# see Recursive Make Considered Harmful, and any number of demos.
AM_INIT_AUTOMAKE([1.13 subdir-objects])

#the language for BertiniReal is C++
AC_LANG([C++])

#another directive to use the m4 folder
AC_CONFIG_MACRO_DIR([m4])

#the only produced file will be a single Makefile.
AC_CONFIG_FILES([Makefile])

#find the CC compiler
AC_PROG_CXX


#enable the creation of shared libraries
AC_ENABLE_SHARED

#enable the creation of static libraries
AC_ENABLE_STATIC

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

#fire up libtool
LT_INIT


AX_CXX_COMPILE_STDCXX_11


#this calls a file in the m4/ directory, which sets up the MPI wrapper stuffs
LX_FIND_MPI

AS_IF([test "$have_CXX_mpi" = "yes"],
	AC_MSG_NOTICE([mpi c++ appears to work correctly]),
	AC_MSG_ERROR([it appears you do not have mpi installed])
	)



# the form of the following commands --
# AC_SEARCH_LIBS(function, libraries-list, action-if-found, action-if-not-found, extra-libraries)

#find gmp
AC_SEARCH_LIBS([__gmpz_init],[gmp], [],[
	AC_MSG_ERROR([unable to find gmp])
])

# find mpfr
AC_SEARCH_LIBS([mpfr_get_version],[mpfr], [],[
	AC_MSG_ERROR([unable to find mpfr])
])



#find bertini
AC_SEARCH_LIBS([head_zero_dim_track_d],[bertini-parallel], [],
	[AC_MSG_ERROR([unable to find libbertini-parallel.  please build and install bertini1 from source.  bertini.nd.edu])],
	$MPI_CXXLDFLAGS)

#find the math library
AC_SEARCH_LIBS([cos], [m], [], [
  AC_MSG_ERROR([unable to find the cos() function])
  ])



AX_BOOST_BASE([1.50],, [AC_MSG_ERROR([bertini_real needs Boost, but it was not found in your system])])
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM

AX_BOOST_REGEX
AX_BOOST_TIMER

AX_BOOST_CHRONO

AM_CONFIG_HEADER(config.h)


#wrap it up.
AC_OUTPUT







