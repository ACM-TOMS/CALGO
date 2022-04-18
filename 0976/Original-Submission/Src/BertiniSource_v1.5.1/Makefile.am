


AM_CFLAGS = -I include

AM_YFLAGS = -d -p `basename $* | sed 's,y$$,,'`
AM_LFLAGS = -s -P`basename $* | sed 's,l$$,,'` -olex.yy.c


ACLOCAL_AMFLAGS = -I m4
#set up the empty variables
bin_PROGRAMS =
BUILT_SOURCES =
lib_LTLIBRARIES =
CLEANFILES =


INSTALL_EXEC_HOOKS =


include_HEADERS =


SUFFIXES = .c .h .y .l
#see https://www.gnu.org/software/automake/manual/html_node/Suffixes.html



#now include the Makemodule files from the subdirectories.  they will add on to the variables whi

include src/Makemodule.am




# formatting note: do NOT put tabs in front of your lines in the if's.
if BUILD_PARALLEL
#set extension of default bertini version to be linked
EXT = parallel
else
#not parallel, but serial
#set extension of default bertini version to be linked
EXT = serial
endif




#set up an install hook to link against the proper version of bertini's executable
install-exec-hook:
	cd $(DESTDIR)$(bindir) && \
	rm -f bertini$(EXEEXT) && \
	$(LN_S) bertini-$(EXT)$(EXEEXT) bertini$(EXEEXT)


#an uninstall hook to remove the symlink (or whatever got created)
uninstall-hook:
	rm $(DESTDIR)$(bindir)/bertini









