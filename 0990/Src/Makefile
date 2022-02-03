# C++ compiler
CC = g++

# C flags
CFLAGS= -w -g -c
CLINKFLAGS = -Wall -g
CSTANDARD = -std=c++11

# Paths
backEndSRC = Source/backend_TOMS_Submission

BIN = bin
backendBIN = bin/backend

EXE = bin/EASAL

SRC = Source
BIN = bin
EXE = $(BIN)/EASAL

# Check if windows
ifeq ($(OS), Windows_NT)
    uname_S := Windows
else
    uname_S := $(shell uname -s)
endif

# Include directories
INC_DIRS = -I ./include/ 

# Libraries
LIBS = -lFOX-1.6 -lglut -lGL -lGLU -lpthread

# Binary objects (add any new *.cpp file)
_OBJ = ActiveConstraintGraph ActiveConstraintRegion Atlas AtlasBuilder \
	AtlasNode Point CartesianRealizer \
	CayleyParameterization CayleyPoint CgMarker \
	ConstraintCheck ConvexChart main PointSet Orientation \
	PredefinedInteractions SaveLoader Settings Stree \
	Utils
OBJ = $(patsubst %, $(backendBIN)/%.o, $(_OBJ))



# Make instructions

backend: $(backendBIN) $(OBJ) $(EXE)
	@echo "Finished."

$(backendBIN): 
	mkdir -p bin/backend 

$(backendBIN)/%.o: $(backEndSRC)/%.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(CSTANDARD) $(INC_DIRS)

$(EXE): $(OBJ)
	$(CC) -o $@ $^ $(CLINKFLAGS) $(LIB_DIRS) 

.PHONY: clean
clean:
	rm -rf $(BIN) 
