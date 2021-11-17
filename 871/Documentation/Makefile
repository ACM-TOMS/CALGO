BML	= /home/wschrep/bin/baseml
MAKE	= make
DELETE	= rm
PAGES	= examples.html download.html index.html paper.html gsl-conversion.html usage.html mpl.html arprec.html
all: $(PAGES)

%.html: %.bml
	$(BML) $?  > $@ 

clean:
	@$(DELETE) -vf $(PAGES) *~ DEADJOE


upload: all
	scp *.html wschrep@hmacs.cmi.ua.ac.be:public_html/precompiler 
	scp css/*.css wschrep@hmacs.cmi.ua.ac.be:public_html/precompiler/css/
	scp -r files/* wschrep@hmacs.cmi.ua.ac.be:public_html/precompiler/files
