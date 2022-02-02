make
g++ -shared -fPIC *.o ../block/*.o -o libgsllinalg.so -L /home/wschrep/gsl-MpIeee/lib  -lMpIeee -lSpecialValue -L . libm.a libc.a

