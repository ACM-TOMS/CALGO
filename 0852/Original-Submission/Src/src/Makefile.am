###########################################################################
# RealPaver v. 0.4                                                        #
###########################################################################
# Author: Laurent Granvilliers                                            #
# Copyright (c) 1999-2003 Institut de Recherche en Informatique de Nantes # 
# Copyright (c) 2004      Laboratoire d'Informatique de Nantes Atlantique #
###########################################################################
# RealPaver is distributed WITHOUT ANY WARRANTY. Read the associated      #
# COPYRIGHT file for more details.                                        #
###########################################################################
# Makefile.am                                                             #
###########################################################################


bin_PROGRAMS=realpaver

realpaver_SOURCES=parser.y lexer.l \
	profile.c interval.c union_interval.c domain.c \
	default_interval.c constant.c variable.c constraint.c evaluator.c \
	narrowing_hullbox.c narrowing_newton.c \
	propagation_local.c propagation_strong.c \
	list_domains.c search.c main.c

realpaver_LDADD=@MATH_LIB@ @LEXLIB@

INCLUDES=-I..

realpaver_DEPENDENCIES=parser.y lexer.l

parser.c:
	$(YACC) -d -oparser.c parser.y

lexer.c:
	$(LEX) -olexer.c lexer.l

clean:
	rm -f *~ *.o realpaver parser.c parser.h lexer.c
