# This macro checks for the Intel Compiler ICC
# that otherwise is mistakenly held for GNU gcc

     AC_DEFUN([AC_PATH_X],
     [AC_REQUIRE_CPP()[]dnl
     AC_MSG_CHECKING([for X])
     # ...omitted...
       AC_MSG_RESULT([libraries $x_libraries, headers $x_includes])
     fi[]dnl
     ])# AC_PATH_X


     # _AC_EMXOS2
     # ----------
     # Check for EMX on OS/2.
m4_define([AC_CXX_ICC],
[AC_CACHE_CHECK([for Intel ICC compiler], [ac_cv_cxx_icc],
[AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [return __EMX__;])],
                        [ac_cv_emxos2=yes],
                        [ac_cv_emxos2=no])])
     test "$ac_cv_emxos2" = yes && EMXOS2=yes[]dnl
     ])# _AC_EMXOS2

