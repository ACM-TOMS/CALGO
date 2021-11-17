# Public-domain work without warranty — see README:LICENSING for details.

vpath %.c ../Src

project=409
cc=gcc
o=-O3

$(project): main.c $(project).c
	$(cc) -Wall -Wextra -Wconversion $o -o $@ $^ -lm

cov: $(project).c.gcov
$(project).c.gcov: main.c $(project).c
	make o='-fprofile-arcs -ftest-coverage' clean $(project)
	./$(project)
	gcov $(project).c
	rm -f *.gcda *.gcno

clean:
	rm -f $(project) $(project).c.gcov
