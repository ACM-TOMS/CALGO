CXX = mpic++

CXXFLAGS = -O2 -std=c++0x -I$(BOOST_ROOT) -I../include

LDFLAGS = 

RM = rm -rf

PROGRAMS = paraperm

SOURCES = $(patsubst %,%.cpp,$(PROGRAMS))
OBJECTS = $(patsubst %,%.o,$(PROGRAMS))

.PHONY: all
all: $(PROGRAMS)

$(PROGRAMS): %: %.o
	$(CXX) -o $@ $< $(LDFLAGS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	$(RM) $(PROGRAMS) $(OBJECTS)
