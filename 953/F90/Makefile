all: lib examples testing

lib:
	( cd Src; make; cd ..; cd Tools; make; cd .. )

examples:
	( cd Examples; make; cd .. )

testing:
	( cd Testing; make; cd .. )

tuning:
	( cd Tuning; make; cd .. )

clean:
	( cd SRC; make clean; \
	cd ../Tools; make clean; \
	cd ../Examples; make clean; \
	cd ../Testing; make clean; \
	cd ../Tuning; make clean; \
	cd .. )
