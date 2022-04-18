include ../make.inc

all:  testlsame $(CRETMG)

testlsame:  lsame.o lsametst.o
	$(LOADER) $(LOADOPT) -o testlsame lsame.o lsametst.o
	testlsame

$(CRETMG): second.o secondtst.o dsecnd.o dsecndtst.o
	$(LOADER) $(LOADOPT) -o testsecond second.o secondtst.o
	$(LOADER) $(LOADOPT) -o testdsecnd dsecnd.o dsecndtst.o
	testsecond; testdsecnd

.f.o: ; $(FORTRAN) $(TSTOPT) -c $<

clean:
	rm -f *.o testlsame testsecond testdsecnd

