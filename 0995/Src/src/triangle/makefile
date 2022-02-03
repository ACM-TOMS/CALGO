SRC = ./
BIN = ./

# CC should be set to the name of your favorite C compiler.
CC = gcc

CSWITCHES = -O -std=gnu11 -I/usr/X11R6/include -L/usr/X11R6/lib

# RM should be set to the name of your favorite rm (file deletion program).
RM = /bin/rm

$(BIN)triangle.o: $(SRC)triangle.c $(SRC)triangle.h
	$(CC) $(CSWITCHES) -DTRILIBRARY -c -o $(BIN)triangle.o \
		$(SRC)triangle.c
distclean:
	$(RM) $(BIN)triangle.o
