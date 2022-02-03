#include "pdsdp5plapack.h"
#include "src/solver/dsdpschurmat_impl.h"
static int printtimes=0;
static int printmore=0;

static void PPDSDPPrintTime(int rank, char* label, double itertime, double cumtime){
  if (printtimes==0) return;
  if (rank==0){
    printf(" %10s:  This iteration: %4.4e, Cumulative Time %4.4e \n",label,itertime,cumtime);
  }
  return;
}

static void PPDSDPPrint(char *label, double dd){
  int rank;
  if (printmore==0) return;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
  if (rank==0){
    printf(" %d: %s %4.4e\n",rank,label,dd);
  }
  return;
}

#include <sys/time.h>
static void wallclock(double * ttime) {
  /* Linux/Unix */
  static struct timeval _tp;
  gettimeofday(&_tp,(struct timezone *)0);
  (*ttime)=((double)_tp.tv_sec)+(1.0e-6)*(_tp.tv_usec);
}

/* PDSDP simply provides a different implementation of the Schur operations. */

typedef struct {
  MPI_Comm mpi_comm, plapack_comm;
  PLA_Obj AMat, vVec, wVec;
  PLA_Obj one,zero,dxerror;
  PLA_Template templ;
  int global_size;
  int nb_distr;
  int rank,nprocs;
  int rowrank,numrownodes;
  int colrank,numcolnodes;
  double ratio;
  double t0,t1,t2,thessian,tsolve;
  DSDP dsdp;
} plapackM;

static struct DSDPSchurMat_Ops plapackdsdpops;
static int DSDPPlapackOpsInit(struct DSDPSchurMat_Ops*);

extern int DSDPSetSchurMatOps(DSDP,struct DSDPSchurMat_Ops*,void*);

#undef __FUNCT__
#define __FUNCT__ "PDSDPUsePLAPACKLinearSolver"
int PDSDPUsePLAPACKLinearSolver(DSDP dsdp, MPI_Comm comm, double ratio, int nb_distr){
  int info;
  plapackM* ctx;

  DSDPFunctionBegin;
  DSDPCALLOC1(&ctx,plapackM,&info);DSDPCHKERR(info);
  info=DSDPPlapackOpsInit(&plapackdsdpops);DSDPCHKERR(info);
  info=DSDPSetSchurMatOps(dsdp,&plapackdsdpops,(void*)ctx);DSDPCHKERR(info);

  ctx->dsdp=dsdp;
  ctx->ratio=ratio;
  ctx->nb_distr=nb_distr;
  ctx->mpi_comm=comm;

  ctx->AMat=NULL;
  ctx->vVec=NULL;
  ctx->wVec=NULL;
  ctx->zero=NULL;
  ctx->one=NULL;
  ctx->dxerror=NULL;
  ctx->templ=NULL;

  DSDPFunctionReturn(0);
}

static int onrow(plapackM* ctx,int row){
  int block,it1,it2,it3;
  block=ctx->nb_distr;
  it1= row / block;  /* Which block */
  it2= (it1) % ctx->numcolnodes;   /* Which colrank does this block belong to */
  it3= (row%block) + block*(it1 / ctx->numcolnodes) ;  /* Of all col */
  if (it2 == ctx->colrank && it3%ctx->numrownodes  == ctx->rowrank ){
    return 1;
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "pmatsetup"
int pmatsetup(void *MM, int m){
  plapackM* ctx=(plapackM*)MM;
  MPI_Comm rowcomm,colcomm;
  int itmp,nprocs,info;

  DSDPFunctionBegin;
  ctx->global_size=m;

  info = MPI_Comm_size(ctx->mpi_comm,&nprocs); DSDPCHKERR(info);
  itmp=(m-nprocs+1)/nprocs;
  itmp=DSDPMax(2,itmp);
  ctx->nb_distr=DSDPMin(ctx->nb_distr,itmp);

  info = PLA_Comm_1D_to_2D_ratio(ctx->mpi_comm,ctx->ratio,&ctx->plapack_comm); DSDPCHKERR(info);
  info = PLA_Init(ctx->plapack_comm); DSDPCHKERR(info);
  info = PLA_Temp_create(ctx->nb_distr, 0, &ctx->templ); DSDPCHKERR(info);

  info=PLA_Matrix_create(MPI_DOUBLE, m, m, ctx->templ,
			 PLA_ALIGN_FIRST, PLA_ALIGN_FIRST, &ctx->AMat);DSDPCHKERR(info);  
  info=PLA_Mvector_create(MPI_DOUBLE, m, 1, ctx->templ, PLA_ALIGN_FIRST, &ctx->vVec);DSDPCHKERR(info);  
  info=PLA_Mvector_create(MPI_DOUBLE, m, 1, ctx->templ, PLA_ALIGN_FIRST, &ctx->wVec);DSDPCHKERR(info);  
  info=PLA_Mscalar_create( MPI_DOUBLE, PLA_ALL_ROWS, PLA_ALL_COLS, 1, 1, ctx->templ, &ctx->dxerror );DSDPCHKERR(info);  
  info=PLA_Mscalar_create( MPI_DOUBLE, PLA_ALL_ROWS, PLA_ALL_COLS, 1, 1, ctx->templ, &ctx->one );DSDPCHKERR(info);
  info=PLA_Mscalar_create( MPI_DOUBLE, PLA_ALL_ROWS, PLA_ALL_COLS, 1, 1, ctx->templ, &ctx->zero );DSDPCHKERR(info);
  info=PLA_Obj_set_to_one(ctx->one);DSDPCHKERR(info);
  info=PLA_Obj_set_to_zero(ctx->zero);DSDPCHKERR(info);

  info = MPI_Comm_rank(ctx->plapack_comm,&ctx->rank); DSDPCHKERR(info);
  info = MPI_Comm_size(ctx->plapack_comm,&ctx->nprocs); DSDPCHKERR(info);

  info = PLA_Temp_comm_col_info(ctx->templ, &rowcomm, &ctx->rowrank, &ctx->numrownodes); DSDPCHKERR(info);
  info = PLA_Temp_comm_row_info(ctx->templ, &colcomm, &ctx->colrank, &ctx->numcolnodes); DSDPCHKERR(info);

  ctx->t0=0;ctx->t1=0;ctx->t2=0;
  ctx->thessian=0;ctx->tsolve=0;
  wallclock(&ctx->t0);
  DSDPFunctionReturn(0);
}

static int pmatzero(void*MM){
  plapackM* ctx=(plapackM*)MM;
  DSDPFunctionBegin;
  wallclock(&ctx->t1);
  PLA_Obj_set_to_zero(ctx->AMat);
  PLA_API_begin();
  PLA_Obj_API_open(ctx->AMat);
  DSDPFunctionReturn(0);
}


static int pmatonprocessor( void* MM,int row ,int*yesorno){
  plapackM* ctx=(plapackM*)MM;
  if (onrow(ctx,row)){
    *yesorno=1;
  } else {
    *yesorno=0;
  }
  return 0;
}

static int pmatlocalvariables( void* MM,double vars[] ,int m){
  plapackM* ctx=(plapackM*)MM;
  int row;
  for (row=0;row<m;row++){
    if (onrow(ctx,row)){
      vars[row]=1;
    } else {
      vars[row]=0;
    }
  }
  return 0;
}

static int pmatmult(void *MM, double x[], double y[], int n){
  plapackM* ctx=(plapackM*)MM;
  double d_one=1.0,drank=1.0/ctx->nprocs;
  int i,info;

  DSDPFunctionBegin;
  info=PLA_Obj_set_to_zero(ctx->vVec);DSDPCHKERR(info);
  info=PLA_Obj_set_to_zero(ctx->wVec);DSDPCHKERR(info);
  info=PLA_API_begin();DSDPCHKERR(info);
  info=PLA_Obj_API_open(ctx->vVec);DSDPCHKERR(info);
  info=PLA_API_axpy_vector_to_global(n, &d_one, x, 1, 
				     ctx->vVec, 0); DSDPCHKERR(info);
  /* Copy solution from PLAPACK vector to DSDPVector */
  info=PLA_Obj_API_close(ctx->vVec); DSDPCHKERR(info);
  info=PLA_API_end(); DSDPCHKERR(info);

  PLA_Symv( PLA_LOWER_TRIANGULAR, ctx->one, ctx->AMat, ctx->vVec, ctx->zero, ctx->wVec ); 
  /* Copy solution from PLAPACK vector to DSDPVector */

  memset((void*)y,0,n*sizeof(double));
  info=PLA_API_begin();
  info=PLA_Obj_API_open(ctx->wVec);

  info=PLA_API_axpy_global_to_vector(n, &d_one, ctx->wVec, 0,
				     y, 1); DSDPCHKERR(info);
  info=PLA_Obj_API_close(ctx->wVec); DSDPCHKERR(info);
  info=PLA_API_end(); DSDPCHKERR(info);
  for (i=0;i<n;i++){ y[i]*=drank;}  /* Should be in PLA_API_axpy_vector_to_global */

  DSDPFunctionReturn(0);
}

static int pmatrowcolumns(void* MM,int row ,double cols[],int*ncols, int nrows){
  plapackM* ctx=(plapackM*)MM;
  int i,nncols=nrows-row;
  DSDPFunctionBegin;
  if (onrow(ctx,row)){
    *ncols=nncols;
    for (i=0;i<nncols;i++){ cols[row+i]=1.0;}
  } else {
    *ncols=0;
  }
  DSDPFunctionReturn(0);
}

static int pmatrowonprocessor(void* MM,int row ,int*flag){
  plapackM* ctx=(plapackM*)MM;
  DSDPFunctionBegin;
  if (onrow(ctx,row)){
    *flag=1;
  } else {
    *flag=0;
  }
  DSDPFunctionReturn(0);
}

static int pmatshiftdiagonal(void* MM,double shift){
  plapackM* ctx=(plapackM*)MM;
  int info,row,m=ctx->global_size;
  double d_one=1.0,shift2;
  DSDPFunctionBegin;
  if (shift==0) return 0;
  info=PLA_API_begin();DSDPCHKERR(info);
  info=PLA_Obj_API_open(ctx->AMat);DSDPCHKERR(info);
  for (row=0;row<m;row++){
    if (onrow(ctx,row)){
      shift2=shift*10;d_one=1.0;
      info = PLA_API_axpy_matrix_to_global(1,1,&d_one, &shift2,
					   ctx->global_size, ctx->AMat, 
					   row, row); DSDPCHKERR(info);
    }
  }
  info=PLA_Obj_API_close(ctx->AMat); DSDPCHKERR(info);
  info=PLA_API_end(); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

static int pmataddelement(void* MM,int row, double shift){
  plapackM* ctx=(plapackM*)MM;
  int info;
  double d_one=1.0;
  DSDPFunctionBegin;
  if (onrow(ctx,row)){
    info = PLA_API_axpy_matrix_to_global(1,1,&d_one,&shift,
					 ctx->global_size, ctx->AMat, 
					 row, row); DSDPCHKERR(info);
  }
  DSDPFunctionReturn(0);
}


static int pmatadddiagonal(void*MM, double val[],int n){
  plapackM* ctx=(plapackM*)MM;
  int row,info,m=ctx->global_size;
  double d_one=1.0;
  DSDPFunctionBegin;
  for (row=0;row<m;row++){
    if (onrow(ctx,row)){
      d_one=1.0;
      info = PLA_API_axpy_matrix_to_global(1,1,&d_one,val+row,
					   ctx->global_size, ctx->AMat, 
					   row, row); DSDPCHKERR(info);
    }
  }
  DSDPFunctionReturn(0);
}

/* Should I cache some lines before inserting them into the matrix ? */
static int pmataddline(void*MM, int row, double dd, double vals[], int m){
  plapackM* ctx=(plapackM*)MM;
  int info;
  double d_one=1.0*dd;
  DSDPFunctionBegin;
  vals[row]+=1.0e-11;  /*  vals[row]*=dd*(1.0 + 1.0e-10); */
  info = PLA_API_axpy_matrix_to_global(m-row,1,&d_one, vals+row,
				       m, ctx->AMat, row, row); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

static int pmatassemble(void*MM){
  plapackM* ctx=(plapackM*)MM;
  int info;
  DSDPFunctionBegin;
  /*  info=pmatshiftdiagonal(MM,1.0e-13); */
  info=PLA_Obj_API_close(ctx->AMat);DSDPCHKERR(info);
  info=PLA_API_end();DSDPCHKERR(info);
  wallclock(&ctx->t2);
  ctx->thessian+=ctx->t2-ctx->t1;
  PPDSDPPrintTime(ctx->rank,"Compute M",ctx->t2-ctx->t1,ctx->thessian);
  DSDPFunctionReturn(0);
}

static int pmatdistributed(void*MM, int*flag){
  DSDPFunctionBegin;
  *flag=1;
  DSDPFunctionReturn(0);
}

static int pmatfactor(void*MM, int *flag){
  plapackM* ctx=(plapackM*)MM;
  int info,dummy;
  double ddxerror;
  DSDPFunctionBegin;
  wallclock(&ctx->t1);
  info=PLA_Obj_set_to_one(ctx->wVec);DSDPCHKERR(info);
  info=PLA_Obj_set_to_zero(ctx->vVec);DSDPCHKERR(info);
  info=PLA_Symv( PLA_LOWER_TRIANGULAR, ctx->one, ctx->AMat, ctx->wVec, ctx->zero, ctx->vVec ); DSDPCHKERR(info);
  *flag=0;
  info = PLA_Chol(PLA_LOWER_TRIANGULAR, ctx->AMat); DSDPCHKERR(info);
  if (info!=0) {
    *flag=1;
    printf("PLAPACK WARNING: Non positive-definite Matrix M : Row: %d\n",info);
  }
  info = PLA_Trsv(PLA_LOWER_TRIANGULAR, PLA_NO_TRANSPOSE, PLA_NONUNIT_DIAG, ctx->AMat, ctx->vVec);DSDPCHKERR(info);
  info = PLA_Trsv(PLA_LOWER_TRIANGULAR, PLA_TRANSPOSE, PLA_NONUNIT_DIAG, ctx->AMat,ctx->vVec); DSDPCHKERR(info);  

  info=PLA_Obj_set_to_minus_one(ctx->wVec);DSDPCHKERR(info); 
  info=PLA_Axpy( ctx->one, ctx->vVec, ctx->wVec );DSDPCHKERR(info); 
  info=PLA_Nrm2( ctx->wVec, ctx->dxerror );DSDPCHKERR(info); 
  PLA_Obj_get_local_contents( ctx->dxerror, PLA_NO_TRANS, &dummy, &dummy,
			      &ddxerror, 1, 1 );
  if (ddxerror/sqrt(1.0*ctx->global_size) > 0.1){
    *flag=1;
    if (ctx->rank==-1){
      printf("PDSDPPLAPACK: Non positive-definite Matrix. %4.2e\n",ddxerror);
    }
  }
  wallclock(&ctx->t2);
  ctx->tsolve+=ctx->t2-ctx->t1;
  PPDSDPPrintTime(ctx->rank,"PLAPACK: Factor M",ctx->t2-ctx->t1,ctx->tsolve);
  PPDSDPPrintTime(ctx->rank,"Subtotal Time",0,ctx->t2-ctx->t1);
  DSDPFunctionReturn(0);
}

static int pmatreduce(void*MM,double *v, int m){
  plapackM* ctx=(plapackM*)MM;
  double d_one=1.0;
  int info;
  int i;
  DSDPFunctionBegin;

  /* 
     Copy vec from DSDPVector to PLAPACK vector. This assumes the entries in the local
     DSDP Vectors are not duplicated on multiple processors. It adds the first element of
     each vector together, second element, ... 
  */
  info=PLA_Obj_set_to_zero(ctx->vVec);DSDPCHKERR(info);
  info=PLA_API_begin();DSDPCHKERR(info);
  info=PLA_Obj_API_open(ctx->vVec);DSDPCHKERR(info);
  info=PLA_API_axpy_vector_to_global(m, &d_one, v , 1, 
                                       ctx->vVec, 0); DSDPCHKERR(info);
  info=PLA_Obj_API_close(ctx->vVec); DSDPCHKERR(info);
  /* Copy solution from PLAPACK vector to DSDPVector */
  memset((void*)v,0,m*sizeof(double));

  info=PLA_Obj_API_open(ctx->vVec);DSDPCHKERR(info);
  info=PLA_API_axpy_global_to_vector(m, &d_one, ctx->vVec, 0,
				     v, 1); DSDPCHKERR(info);
  info=PLA_Obj_API_close(ctx->vVec); DSDPCHKERR(info);
  info=PLA_API_end(); DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}

static int pmatsolve(void* MM, double bb[], double xx[], int n){
  plapackM* ctx=(plapackM*)MM;
  double d_one=1.0,drank=1.0/ctx->nprocs;
  int i,info;

  DSDPFunctionBegin;
  wallclock(&ctx->t1);
  /* 
     Copy RHS from DSDPVector to PLAPACK vector. This assumes the entries in the local
     DSDP Vectors are not duplicated on multiple processors. It adds the first element of
     each vector together, second element, ... 
  */
  info=PLA_Obj_set_to_zero(ctx->vVec);DSDPCHKERR(info);
  info=PLA_API_begin();DSDPCHKERR(info);
  info=PLA_Obj_API_open(ctx->vVec);DSDPCHKERR(info);
  info=PLA_API_axpy_vector_to_global(n, &d_one, bb , 1, 
				     ctx->vVec, 0); DSDPCHKERR(info);
  info=PLA_Obj_API_close(ctx->vVec); DSDPCHKERR(info);
  info=PLA_API_end(); DSDPCHKERR(info);


  /* Assuming the matrix is already factored, solve the equations. */
  info = PLA_Trsv(PLA_LOWER_TRIANGULAR, PLA_NO_TRANSPOSE, PLA_NONUNIT_DIAG, ctx->AMat, ctx->vVec);DSDPCHKERR(info);
  info = PLA_Trsv(PLA_LOWER_TRIANGULAR, PLA_TRANSPOSE, PLA_NONUNIT_DIAG, ctx->AMat,ctx->vVec); DSDPCHKERR(info);  

  /* Copy solution from PLAPACK vector to DSDPVector */
  memset((void*)xx,0,n*sizeof(double));
  info=PLA_API_begin();
  info=PLA_Obj_API_open(ctx->vVec);
  info=PLA_API_axpy_global_to_vector(n, &d_one, ctx->vVec, 0,
				     xx, 1); DSDPCHKERR(info);
  info=PLA_Obj_API_close(ctx->vVec); DSDPCHKERR(info);
  info=PLA_API_end(); DSDPCHKERR(info);
  for (i=0;i<n;i++){xx[i]*=drank;}

  wallclock(&ctx->t2);
  ctx->tsolve+=ctx->t2-ctx->t1;
  /* PPDSDPPrintTime(ctx->rank,"Solve M",ctx->t2-ctx->t1,ctx->tsolve);*/
  DSDPFunctionReturn(0);
}

static int pmatdestroy(void*MM){
  plapackM* ctx=(plapackM*)MM;
  int info;
  DSDPFunctionBegin;
  PPDSDPPrint("PDSDP: Compute M: %5.5e \n",ctx->thessian);
  PPDSDPPrint("PDSDP: Solve M %4.5e seconds\n ",ctx->tsolve);
  info=PLA_Obj_free(&ctx->AMat); DSDPCHKERR(info);
  info=PLA_Obj_free(&ctx->vVec); DSDPCHKERR(info);
  info=PLA_Obj_free(&ctx->wVec); DSDPCHKERR(info);
  info=PLA_Obj_free(&ctx->one); DSDPCHKERR(info);
  info=PLA_Obj_free(&ctx->zero); DSDPCHKERR(info);
  info=PLA_Obj_free(&ctx->dxerror); DSDPCHKERR(info);
  info=PLA_Temp_free(&ctx->templ); DSDPCHKERR(info);
  info=PLA_Finalize(); DSDPCHKERR(info);
  DSDPFREE(&ctx,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


static const char *plapackmatname="PLAPACK matrix";

static int DSDPPlapackOpsInit(struct DSDPSchurMat_Ops* sops){
  int info;
  if (!sops) return 0;
  info=DSDPSchurMatOpsInitialize(sops); DSDPCHKERR(info);
  sops->matzero=pmatzero;
  sops->matrownonzeros=pmatrowcolumns;
  sops->mataddrow=pmataddline;
  sops->matadddiagonal=pmatadddiagonal;
  sops->mataddelement=pmataddelement ;
  sops->matshiftdiagonal=pmatshiftdiagonal;
  sops->matassemble=pmatassemble;
  sops->matscaledmultiply=0;
  sops->matfactor=pmatfactor;
  sops->matsolve=pmatsolve;
  sops->pmatwhichdiag=0;
  sops->matscaledmultiply=pmatmult;
  sops->pmatonprocessor=pmatonprocessor;
  sops->pmatlocalvariables=pmatlocalvariables;
  sops->pmatdistributed=pmatdistributed;
  sops->pmatreduction=pmatreduce;
  sops->matsetup=pmatsetup;
  sops->matdestroy=pmatdestroy;
  sops->id=25;
  sops->matname=plapackmatname;
  return 0;
}
