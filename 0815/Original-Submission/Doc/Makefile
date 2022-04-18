FC = /opt/NAGWare/bin/f95
FFLAGS = -g -gline -u -dcfuns -nan -save -C -C=undefined

OBJ1 = driver1.o src1.o
OBJ2 = driver2.o src2.o

all: res1.1 res1.2 res2.1 res2.2

res1.1: obj1
	obj1 < data1 > res1.1

res1.2: obj2
	obj2 < data1 > res1.2

res2.1: obj1
	obj1 < data2 > res2.1

res2.2: obj2
	obj2 < data2 > res2.2
obj2 : $(OBJ2)
	$(FC) $(FFLAGS) $(OBJ2) -o obj2

obj1 : $(OBJ1)
	$(FC) $(FFLAGS) $(OBJ1) -o obj1
