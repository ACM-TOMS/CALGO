make
g++ -Wall -O2 -shared -fpic *.o -o libgsllinalg.so -L /home/wschrep/gsl-MpIeee/lib  -lMpIeee -lSpecialValue -lgmp 

