# Example of a UNIX makefile provided to
# UNIX users for a simple compilation and linkage of separate files of
# the application project.
#
# Authors:                                                                
#
# M.Gaviano, D.E.Kvasov, D.Lera, and Ya.D.Sergeyev                     
#
gkls : example.o gkls.o rnd_gen.o
	$(CC) $(CFLAGS) -lm example.o gkls.o rnd_gen.o -o gkls
example.o : ../Src/gkls.h ../Src/rnd_gen.h example.c
	$(CC) -I. $(CFLAGS) -c example.c
gkls.o : ../Src/gkls.h ../Src/rnd_gen.h ../Src/gkls.c
	$(CC) -I. $(CFLAGS) -c ../Src/gkls.c
rnd_gen.o : ../Src/rnd_gen.h ../Src/rnd_gen.c
	$(CC) -I. $(CFLAGS) -c ../Src/rnd_gen.c
clean:
	-rm *.o
