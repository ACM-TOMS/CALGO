#
# Available options are as follows (everything is disabled by default):
#
# ENABLE_LOG    : enables log messages (too verbose for big examples)
# ENABLE_INFO   : enables info messages
# ENABLE_STAT   : enables operation count statistics
# ENABLE_TIME   : enables time measurements
# ENABLE_LONG   : use long integers instead of basic integers
#
CXXFLAGS	?= -DENABLE_TIME -DENABLE_STAT #-DENABLE_INFO

#CXXFLAGS	+= -DENABLE_LONG

CXXFLAGS	+= -O3
#CXXFLAGS	+= -Og -g -DENABLE_LOG

CXXFLAGS	+= -std=c++11 -fopenmp -Iinclude -Wall -Wextra -pedantic
LDFLAGS		+=
LDLIBS		+=

TOOLS		:= $(patsubst tool/%.cpp,bin/%,$(wildcard tool/*.cpp))
TESTS		:= $(patsubst test/%.cpp,build/test/%,$(wildcard test/*.cpp))

DEPS		:= $(patsubst tool/%.cpp,build/dep/%.d,$(wildcard tool/*.cpp))
DEPS		+= $(patsubst test/%.cpp,build/dep/%.d,$(wildcard test/*.cpp))

all: dirs $(TOOLS)

dirs:
	$(shell mkdir -p build/dep build/test bin)

-include $(DEPS)

bin/%: tool/%.cpp Makefile
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@ -MMD -MF $(@:bin/%=build/dep/%.d) -MT "$@ $(@:bin/%=build/dep/%.d)"

build/test/%: test/%.cpp Makefile
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@ -MMD -MF $(@:build/test/%=build/dep/%.d) -MT "$@ $(@:build/test/%=build/dep/%.d)"

clean:
	rm -rf build bin

run: dirs $(TOOLS)
	bin/cppr test/sample/inp.max

test: dirs $(TESTS)
	@for t in $(TESTS); do echo $$t; $$t; done

.PHONY: all dirs clean run test
