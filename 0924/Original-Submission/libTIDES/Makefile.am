AM_CPPFLAGS = -I/usr/local/include

lib_LIBRARIES = libTIDES.a

libTIDES_a_SOURCES = commonITER.c commonITER.h \
                     doubITER.c doubITER.h \
                     doubNUM.c doubNUM.h doubNUMdef.h \
                     taylorODE.c taylorODE.h

if HAVE_MP
  libTIDES_a_SOURCES += mpfrITER.c mpfrITER.h \
                        mpfrNUM.c mpfrNUM.h mpfrNUMdef.h \
                        mpfrtaylorODE.c mpfrtaylorODE.h
endif
