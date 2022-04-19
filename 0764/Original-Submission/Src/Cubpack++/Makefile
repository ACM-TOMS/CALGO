SRC_DIR = Src

install:
	./quick_install

run:
	./quick_run
	mv $(SRC_DIR)/tst .
	echo 'executable "tst" is available.'

clean:
	cd $(SRC_DIR) ; /bin/rm -f *.o  tst core makefile
	cd $(SRC_DIR) ; /bin/rm -f main.c libcubpack.a
	cd $(SRC_DIR) ; /bin/rm -r -f ptrepository tempinc Templates.DB

