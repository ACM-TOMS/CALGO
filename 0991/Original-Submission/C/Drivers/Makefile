executables =  driver2d stability timing

all: $(executables)

CC = gcc
CFLAGS = -std=gnu99 -g -Wall -O3 -I ../Src 
LDLIBS = -lm

vpath %.c ../Src
vpath %.h ../Src

OBJS = tswdft2d.o utils.o 
ALLOBJS = $(OBJS) swdft2d.o swfft2d.o

driver2d: driver2d.o $(OBJS)
	$(CC) $(CFLAGS) -I ../Src -o $@ $^ $(LDLIBS)

stability: stability.o $(ALLOBJS)
	$(CC) $(CFLAGS) -I ../Src -o $@ $^ $(LDLIBS)

timing: timing.o $(ALLOBJS)
	$(CC) $(CFLAGS) -I ../Src -o $@ $^ $(LDLIBS)

clean: 
	rm -f $(ALLOBJS) 
	rm -f driver2d.o stability.o timing.o $(executables)
