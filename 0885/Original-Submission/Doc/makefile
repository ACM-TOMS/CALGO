CC = gcc
CFLAGS = -lm


main: main.c lnnorm.c normal.c log1px.c anorm.c enorm.c lnanorm.c lnenorm.c
	$(CC) $(CFLAGS) -o main main.c lnnorm.c normal.c log1px.c anorm.c enorm.c lnanorm.c lnenorm.c

clean:
	rm -f *.o  *.itex
	rm -f main
