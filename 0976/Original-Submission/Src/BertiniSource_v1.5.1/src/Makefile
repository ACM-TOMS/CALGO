# This file is deprecated
INC=-I../include 
OPT=-O3 -funroll-loops -fexpensive-optimizations -Wall #-D_HAVE_MPI
ARGS=$(INC) $(OPT) 
OBJFILES=misc.o track.o prog_eval.o eval_functions.o endgame.o post_processing.o matrixSolve.o parse_input.o startpoint_maker.o setup_functions.o rank_finder.o zero_dim_track.o zero_dim_setup.o zero_dim_main.o clear_functions.o QR.o svd_from_QR.o PSEG.o sharpen.o svd_analysis.o ode.o parallel_datatypes.o copy_functions.o parallel_send_functions.o parallel_main.o min_svd.o adaptiveMP2.o PSEG_amp.o CauchyEG.o Cauchy_rank.o corank.o regen_eval.o regen_setup.o regen_track.o regen_parallel.o eqbyeq_eval.o eqbyeq_setup.o eqbyeq_output.o eqbyeq_track.o eqbyeq_parallel.o pos_dim.o pureDecomp.o numericalIrredDecomp.o membership.o membership_eval.o membership_junk.o deflation.o dimbydim_eval.o dimbydim_output.o dimbydim_setup.o dimbydim_track.o dimbydim_parallel.o cascade_eval.o cascade_output.o cascade_setup.o cascade_track.o cascade_parallel.o QR_block.o corank_rrv.o ldt_eval.o ldt_mult_mat.o ldt_main.o ldt_rank.o ldt_junk_removal.o ldt_dual_basis.o QLP.o NID_parallel.o CauchyEG_new.o parallel_endgame.o trackBackEG.o regen_pos_dim_setup.o regen_pos_dim_track.o regen_pos_dim_eval.o regen_pos_dim_copy.o regen_pos_dim_parallel.o ldt_diff.o diff_deflatable.o track_error.o adaptiveMP_error.o setupSLP.o setupSLP_defSF.o parallel_diff.o function_main.o newton_main.o witnessGeneration.o isosingular.o witnessProjection.o printWitness.o generalSliceMoving.o regen_pos_dim_extend.o paramHomotopy.o 
POBJFILES=lex.splitParse.o lex.ppParse.o ppParse.tab.o lex.pParse.o pParse.tab.o
COMP=gcc
FLEX=flex
BISON=bison

all : $(OBJFILES) $(POBJFILES) bertini libbertini;

libbertini: $(OBJFILES) $(POBJFILES)
	$(COMP) $(INC) -Dmain=bertini_main $(OPT) -c bertini.c -o bertini.o 
	ar rcs libbertini.a bertini.o $(OBJFILES) $(POBJFILES)

bertini : bertini.c $(OBJFILES) $(POBJFILES) ;
	$(COMP) $(ARGS) -o bertini bertini.c $(OBJFILES) $(POBJFILES) -static -lm -lmpfr -lgmp

.c.o : 
	$(COMP) $(ARGS) -c $*.c 

lex.pParse.o pParse.tab.o: lex.pParse.c pParse.tab.c
	$(COMP) $(ARGS) -c lex.pParse.c pParse.tab.c

pParse.tab.c :  pParse.y
	$(BISON) -p pParse -b pParse -d pParse.y

lex.pParse.c:   pParse.l
	$(FLEX) -PpParse pParse.l

lex.ppParse.o ppParse.tab.o: lex.ppParse.c ppParse.tab.c
	$(COMP) $(ARGS) -c lex.ppParse.c ppParse.tab.c

ppParse.tab.c :  ppParse.y
	$(BISON) -p ppParse -b ppParse -d ppParse.y

lex.ppParse.c:   ppParse.l
	$(FLEX) -PppParse ppParse.l

lex.splitParse.o:	lex.splitParse.c
	$(COMP) $(ARGS) -c lex.splitParse.c

lex.splitParse.c:	splitParse.l
	$(FLEX) -PsplitParse splitParse.l

clean : ;
	rm -f $(OBJFILES) $(POBJFILES) bertini SV_cmp_d SV_cmp_mp gmon.out deg.out names.out ODE_data
	rm -f par.out paramDerivs.out arr.out eval.out eval2.out const.out finalFile.out jacV.out jacP.out num.out
	rm -f lex.parser.c mpf.out parser.tab.c parser.tab.h compareR compareI odeOUT sys.out tempBertiniFile
	rm -f lex.parse_and_hom.c parse_and_hom.tab.c parse_and_hom.tab.h lex.splitParse.c lex.ppParse.c ppParse.tab.c ppParse.tab.h lex.pParse.c pParse.tab.c pParse.tab.h
	rm -f libbertini.* bertini.o
