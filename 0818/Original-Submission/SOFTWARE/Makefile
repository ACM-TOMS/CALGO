#
FFLAGS = -g
CFLAFS =
#
LDFLAGS = 
#
FC = f90
CC = cc
LD = $(FC)
AR = ar -r -v
RANLIB =
RM = rm -f
#
MOD_SUF = mod
###############################################################################
OBJS = conv_tools.o dense.o hash.o info.o link.o\
	lmbv_coo.o lmbv_csc.o lmbv_csr.o lmbv_dia.o\
	lmbv_bco.o lmbv_bsc.o lmbv_bsr.o lmbv_bdi.o lmbv_vbr.o\
	lsbv_coo.o lsbv_csc.o lsbv_csr.o lsbv_dia.o\
	lsbv_bco.o lsbv_bsc.o lsbv_bsr.o lsbv_bdi.o lsbv_vbr.o\
	mbv.o blas_sparse_namedconstants.o properties.o\
	rmbv_coo.o rmbv_csc.o rmbv_csr.o rmbv_dia.o\
	rmbv_bco.o rmbv_bsc.o rmbv_bsr.o rmbv_bdi.o rmbv_vbr.o\
	rsbv_coo.o rsbv_csc.o rsbv_csr.o rsbv_dia.o\
	rsbv_bco.o rsbv_bsc.o rsbv_bsr.o rsbv_bdi.o rsbv_vbr.o\
	sbv.o types.o\
	usconv_coo2csr.o usconv_csr2coo.o usconv_csc2coo.o usconv_coo2csc.o\
        usconv_coo2dia.o usconv_dia2coo.o usconv_bco2bsr.o usconv_bco2bsc.o\
        usconv_bsr2bco.o usconv_bsc2bco.o usconv_bco2bdi.o usconv_bdi2bco.o\
        uscr_coo.o uscr_csc.o uscr_csr.o uscr_dia.o  \
        uscr_bco.o uscr_bsc.o uscr_bsr.o uscr_bdi.o uscr_vbr.o\
        uscr_begin.o  uscr_block_begin.o uscr_variable_block_begin.o\
        ussp.o usgp.o\
	uscr.o usds.o usmv.o usmm.o ussv.o ussm.o\
        usdot.o usaxpy.o usga.o usgz.o ussc.o\
        SparseBLAS.o  SparseBLAS1.o\
        Entry.o INS_ROUTINER.o  INSERTING.o uscr_insert_block.o\
        uscr_insert_col.o uscr_insert_row.o  uscr_insert_entries.o\
        uscr_insert_clique.o  uscr_insert_entry.o uscr_end.o\
        blas_sparse_proto.o blas_sparse.o
COBJS = 
LIBS = 
###############################################################################
libSparseBLAS_$(SBLAS_ARCH).a:$(OBJS) 
	$(AR) $@ $(OBJS)
###############################################################################
blas_sparse.o : blas_sparse_namedconstants.o blas_sparse_proto.o
blas_sparse_proto.o: SparseBLAS1.o Entry.o
conv_tools.o : blas_sparse_namedconstants.o
dense.o : properties.o
Entry.o: uscr_begin.o uscr_insert_entry.o  uscr_insert_col.o uscr_insert_row.o uscr_insert_entries.o uscr_end.o uscr_block_begin.o  uscr_insert_block.o  uscr_variable_block_begin.o uscr_insert_clique.o INS_ROUTINER.o INSERTING.o ussp.o  usgp.o
hash.o : blas_sparse_namedconstants.o
info.o : link.o properties.o types.o
INS_ROUTINER.o: INSERTING.o  SparseBLAS1.o properties.o
INSERTING.o : properties.o blas_sparse_namedconstants.o
link.o : properties.o types.o
mbv.o : lmbv_coo.o rmbv_coo.o lmbv_csc.o rmbv_csc.o lmbv_csr.o rmbv_csr.o\
        lmbv_dia.o rmbv_dia.o rmbv_bsr.o lmbv_bsr.o lmbv_bsc.o rmbv_bsc.o\
        lmbv_bdi.o rmbv_bdi.o lmbv_vbr.o rmbv_vbr.o lmbv_bco.o rmbv_bco.o
lmbv_coo.o : properties.o link.o
rmbv_coo.o : properties.o link.o
lmbv_csc.o : properties.o link.o 
rmbv_csc.o : properties.o link.o 
lmbv_csr.o : properties.o link.o 
rmbv_csr.o : properties.o link.o 
lmbv_dia.o : properties.o link.o 
rmbv_dia.o : properties.o link.o 
lmbv_bco.o : properties.o link.o dense.o
rmbv_bco.o : properties.o link.o dense.o
lmbv_bsr.o : properties.o link.o dense.o
rmbv_bsr.o : properties.o link.o dense.o
lmbv_bsc.o : properties.o link.o dense.o
rmbv_bsc.o : properties.o link.o dense.o
lmbv_bdi.o : properties.o link.o dense.o
rmbv_bdi.o : properties.o link.o dense.o
lmbv_vbr.o : properties.o link.o dense.o
rmbv_vbr.o : properties.o link.o dense.o
lsbv_coo.o : properties.o link.o hash.o
rsbv_coo.o : properties.o link.o hash.o
lsbv_csc.o : properties.o link.o
rsbv_csc.o : properties.o link.o
lsbv_csr.o : properties.o link.o
rsbv_csr.o : properties.o link.o
lsbv_dia.o : properties.o link.o
rsbv_dia.o : properties.o link.o
lsbv_bco.o : properties.o link.o hash.o dense.o
rsbv_bco.o : properties.o link.o hash.o dense.o
lsbv_bsr.o : properties.o link.o dense.o
rsbv_bsr.o : properties.o link.o dense.o
lsbv_bsc.o : properties.o link.o dense.o
rsbv_bsc.o : properties.o link.o dense.o
lsbv_bdi.o : properties.o link.o dense.o
rsbv_bdi.o : properties.o link.o dense.o
lsbv_vbr.o : properties.o link.o dense.o
rsbv_vbr.o : properties.o link.o dense.o
properties.o : blas_sparse_namedconstants.o
sbv.o : lsbv_coo.o rsbv_coo.o lsbv_csc.o rsbv_csc.o lsbv_csr.o rsbv_csr.o\
        lsbv_dia.o rsbv_dia.o lsbv_bsr.o rsbv_bsr.o lsbv_bsc.o rsbv_bsc.o\
        lsbv_bdi.o rsbv_bdi.o lsbv_vbr.o rsbv_vbr.o lsbv_bco.o rsbv_bco.o
SparseBLAS.o: uscr.o usds.o usmv.o ussv.o usmm.o ussm.o info.o
SparseBLAS1.o: SparseBLAS.o usconv_coo2csr.o usconv_coo2csc.o usconv_csr2coo.o usconv_csc2coo.o usconv_coo2dia.o usconv_dia2coo.o usconv_bco2bsr.o usconv_bco2bsc.o usconv_bsr2bco.o usconv_bsc2bco.o usconv_bco2bdi.o usconv_bdi2bco.o usdot.o usaxpy.o  usga.o usgz.o ussc.o properties.o
types.o : blas_sparse_namedconstants.o
usaxpy.o : blas_sparse_namedconstants.o
usconv_coo2csr.o :conv_tools.o properties.o link.o 
usconv_coo2csc.o :conv_tools.o properties.o link.o
usconv_csr2coo.o :conv_tools.o properties.o link.o
usconv_csc2coo.o :conv_tools.o properties.o link.o
usconv_coo2dia.o :conv_tools.o properties.o link.o
usconv_dia2coo.o :conv_tools.o properties.o link.o
usconv_bco2bsr.o :conv_tools.o properties.o link.o
usconv_bco2bsc.o :conv_tools.o properties.o link.o
usconv_bsr2bco.o :conv_tools.o properties.o link.o
usconv_bsc2bco.o :conv_tools.o properties.o link.o
usconv_bco2bdi.o :conv_tools.o properties.o link.o
usconv_bdi2bco.o :conv_tools.o properties.o link.o
uscr.o : uscr_coo.o uscr_csc.o uscr_csr.o uscr_dia.o\
         uscr_bco.o uscr_bsc.o uscr_bsr.o uscr_bdi.o uscr_vbr.o
uscr_coo.o : properties.o link.o usds.o
uscr_csc.o : properties.o link.o usds.o
uscr_csr.o : properties.o link.o usds.o
uscr_dia.o : properties.o link.o usds.o
uscr_bco.o : properties.o link.o usds.o
uscr_bsr.o : properties.o link.o usds.o
uscr_bsc.o : properties.o link.o usds.o
uscr_bdi.o : properties.o link.o usds.o
uscr_vbr.o : properties.o link.o usds.o
uscr_end.o : INS_ROUTINER.o INSERTING.o properties.o
uscr_begin.o: INSERTING.o properties.o blas_sparse_namedconstants.o
uscr_insert_block.o: INS_ROUTINER.o INSERTING.o blas_sparse_namedconstants.o
uscr_block_begin.o: INSERTING.o properties.o blas_sparse_namedconstants.o
uscr_variable_block_begin.o:INSERTING.o properties.o blas_sparse_namedconstants.o
uscr_insert_entry.o: INS_ROUTINER.o INSERTING.o blas_sparse_namedconstants.o
uscr_insert_col.o : uscr_insert_entry.o blas_sparse_namedconstants.o
uscr_insert_row.o : uscr_insert_entry.o blas_sparse_namedconstants.o
uscr_insert_clique.o: uscr_insert_entry.o blas_sparse_namedconstants.o
uscr_insert_entries.o: uscr_insert_entry.o blas_sparse_namedconstants.o
usdot.o : blas_sparse_namedconstants.o
usga.o : blas_sparse_namedconstants.o
ussc.o : blas_sparse_namedconstants.o
usgz.o : blas_sparse_namedconstants.o
ussp.o:INSERTING.o properties.o
usgp.o:INSERTING.o properties.o
usds.o : link.o
usmv.o : properties.o link.o mbv.o
usmm.o : properties.o link.o mbv.o
ussv.o : properties.o link.o sbv.o
ussm.o : properties.o link.o sbv.o
###############################################################################
.SUFFIXES: .o .F .c .f .f90

.f.o :
	$(FC) $(FFLAGS) -c $*.f

.F.o :
	$(FC) $(FFLAGS) -c $*.F
.f90.o :
	$(FC) $(FFLAGS) -c $*.f90

.c.o :
	$(CC) $(CFLAGS) -c $*.c

clean :
	$(RM) *.o *.$(MOD_SUF) libSparseBLAS*.a
