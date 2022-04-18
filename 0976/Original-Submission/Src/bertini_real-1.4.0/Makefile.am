

AM_YFLAGS = -d -p `basename $* | sed 's,y$$,,'`
AM_LFLAGS = -s -P`basename $* | sed 's,l$$,,'` -olex.yy.c

AM_CPPFLAGS = -I$(top_srcdir)/include -D_HAVE_MPI $(BOOST_CPPFLAGS) $(MPI_CXXFLAGS)

ACLOCAL_AMFLAGS = -I m4

#set up the empty variables
bin_PROGRAMS =
BUILT_SOURCES =
CLEANFILES =
noinst_LIBRARIES =
noinst_LTLIBRARIES =


#installed headers
include_HEADERS =
rootinclude_HEADERS =

SUFFIXES = .c .h .y .l
#see https://www.gnu.org/software/automake/manual/html_node/Suffixes.html


rootincludedir = $(includedir)/bertini_real

#now include the Makemodule files from the subdirectories.  they will add on to the variables which are declared above.



# some sublibraries
include src/Makemodule.am  					#unclassified

include src/bertini1/Makemodule.am 			#bertini1
include src/cells/Makemodule.am 			#cells
include src/containers/Makemodule.am 		#containers
include src/decompositions/Makemodule.am 	#decompositions
include src/io/Makemodule.am 				#io
include src/nag/Makemodule.am 				#nag
include src/symbolics/Makemodule.am 		#symbolics


noinst_LTLIBRARIES += libbr_common.la
libbr_common_ladir = $(includedir)/bertini_real

libbr_common_la_SOURCES = \
	$(unclassified_sources) \
	$(bertini1_sources) \
	$(cells_sources) \
	$(containers_sources) \
	$(decompositions_sources) \
	$(io_sources) \
	$(nag_sources) \
	$(symbolics_sources)

libbr_common_la_SOURCES += \
	$(unclassified_headers) \
	$(bertini1_headers) \
	$(cells_headers) \
	$(containers_headers) \
	$(decompositions_headers) \
	$(io_headers) \
	$(nag_headers) \
	$(symbolics_headers)

libbr_common_la_LIBADD = libpartitionParse.la $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB) $(BOOST_CHRONO_LIB) $(BOOST_REGEX_LIB) $(BOOST_TIMER_LIB) $(MPI_CXXLDFLAGS)


# the two programs built from this source code
include src/bertini_real/Makemodule.am
include src/sampler/Makemodule.am








