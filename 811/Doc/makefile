OBJS = src.o subs.o
all: res1 res2 res3 res4 res5 res6 res7 res8
res2: driver2.o ${OBJS} 
	${FC} ${FFLAGS} tbunu.o ${OBJS}
	a.out < data1 > $@
res1: driver1.o ${OBJS} 
	${FC} ${FFLAGS} tbunl.o ${OBJS}
	a.out > $@
res4: driver4.o ${OBJS} 
	${FC} ${FFLAGS} tminu.o ${OBJS}
	a.out < data1 > $@
res3: driver3.o ${OBJS} 
	${FC} ${FFLAGS} tminl.o ${OBJS}
	a.out > $@
res6: driver6.o ${OBJS} 
	${FC} ${FFLAGS} tnewu.o ${OBJS}
	a.out < data1 > $@
res5: driver5.o ${OBJS} 
	${FC} ${FFLAGS} tnewl.o ${OBJS}
	a.out > $@
res8: driver8.o ${OBJS} 
	${FC} ${FFLAGS} tnewu.o ${OBJS}
	a.out < data1 > $@
res7: driver7.o ${OBJS} 
	${FC} ${FFLAGS} tnewl.o ${OBJS}
	a.out > $@
