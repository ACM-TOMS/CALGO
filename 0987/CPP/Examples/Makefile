CXXFLAGS=-Wall 
#-Weffc++
CXX = g++
LDFLAGS=

SUBDIRS = Example1 Example2 Fun Fun100 Fun1000

all: example1 example2 fun fun100 fun1000

.PHONY: example1 example2 fun

example1: 
	cd Example1 && $(MAKE) example1

example2:
	cd Example2 && $(MAKE) example2

fun:
	cd Fun && $(MAKE) fun

fun100:
	cd Fun100 && $(MAKE) fun100

fun1000:
	cd Fun1000 && $(MAKE) fun1000

clean: 
	for subdir in $(SUBDIRS) ; do \
	  (cd $$subdir && $(MAKE) clean); \
	done
