include ../make.inc

all:      TDGESQB TDGESQR TDGESUB TDGESUV TDHABAL \
          TDHAESU TDHAORD TDHAPVB TDHAPVL TDHASUB \
          TDHGPQR TDSHBAL TDSHES  TDSHPVB TDSHPVL

clean:
	rm -f *.o

TDGESQB:  tdgesqb.o ; $(LOADER) -o tdgesqb tdgesqb.o $(LOADOPTS)
TDGESQR:  tdgesqr.o ; $(LOADER) -o tdgesqr tdgesqr.o $(LOADOPTS)
TDGESUB:  tdgesub.o ; $(LOADER) -o tdgesub tdgesub.o $(LOADOPTS)
TDGESUV:  tdgesuv.o ; $(LOADER) -o tdgesuv tdgesuv.o $(LOADOPTS)
TDHABAL:  tdhabal.o ; $(LOADER) -o tdhabal tdhabal.o $(LOADOPTS)
TDHAESU:  tdhaesu.o ; $(LOADER) -o tdhaesu tdhaesu.o $(LOADOPTS)
TDHAORD:  tdhaord.o ; $(LOADER) -o tdhaord tdhaord.o $(LOADOPTS)
TDHAPVB:  tdhapvb.o ; $(LOADER) -o tdhapvb tdhapvb.o $(LOADOPTS)
TDHAPVL:  tdhapvl.o ; $(LOADER) -o tdhapvl tdhapvl.o $(LOADOPTS)
TDHASUB:  tdhasub.o ; $(LOADER) -o tdhasub tdhasub.o $(LOADOPTS)
TDHGPQR:  tdhgpqr.o ; $(LOADER) -o tdhgpqr tdhgpqr.o $(LOADOPTS)
TDSHBAL:  tdshbal.o ; $(LOADER) -o tdshbal tdshbal.o $(LOADOPTS)
TDSHES:   tdshes.o  ; $(LOADER) -o tdshes  tdshes.o  $(LOADOPTS)
TDSHPVB:  tdshpvb.o ; $(LOADER) -o tdshpvb tdshpvb.o $(LOADOPTS)
TDSHPVL:  tdshpvl.o ; $(LOADER) -o tdshpvl tdshpvl.o $(LOADOPTS)


.f.o: ; $(FORTRAN) $(OPTS) -c $<
