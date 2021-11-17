#include "pdsdp5scalapack.h"
#include "dsdpsys.h"
#include "src/solver/dsdpschurmat_impl.h"

#include "pblas.h"
#include "PBpblas.h"
#include "PBtools.h"
#include "PBblas.h"
#include "PBblacs.h"


typedef long int ftnlen;
/* External Scalapack and BLACS routines from Fortran routines */
extern int pdpotrf_(char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info, ftnlen uplo_len);
extern int pdpotrs_(char *uplo, int *n, int *nrhs, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb, int *info, ftnlen uplo_len);
extern int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern int descinit_(int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld, int *info);
extern void Cpdgemr2d(int m, int n, double* ptrmyblock, int ia, int ja, int* desca, double* ptrmynewblock, int  ib, int jb, int* descb, int globcontext);

extern void Cpdtrmr2d(char* uplo, char* diag,int m, int n, double* ptrmyblock, int ia, int ja, int* desca,double* ptrmynewblock, int  ib, int jb, int* descb, int globcontext);

extern int pdlaprnt_(int *m, int *n, double *a, int *ia, int *ja, int *desca, int *irprnt, int *icprnt, char *cmatnm, int *nout, double *work, ftnlen cmatnm_len);

static int pvecreduce(void*,double *, int);

static int printtimes=0;
static int printmore=1;
static void PPDSDPPrintTime(int rank, char* label, double itertime, double cumtime){
  if (printtimes==0) return;
  if (rank==0){
    printf(" %10s:  This iteration: %4.4e, Cumulative Time %4.4e \n",label,itertime,cumtime);
  }
  return;
}

static void PPDSDPPrint(char *label, double dd){
  int rank=0;
  if (printmore==0) return;

  MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
  if (rank==0){
    printf(" %d: ",rank);
    printf(label,dd);
  }
  return;
}


/* PDSDP simply provides a different implementation of the Schur operations. */

typedef struct {

  int m;
  int rank;
  double ratio;
  int fflag;

  double t0,t1,t2,thessian,tfactor,tsolve;

  /*
    mb defines how many continuous rows are processed by one processor.
    mb2 defines block size of two cyclic distribution of ScaLAPACK
  */
  int nprocs;

  int ictxt;
  int myrow;
  int mycol;
  int nprow;
  int npcol;
  int mb;

  int ictxt2;
  int myrow2;
  int mycol2;
  int nprow2;
  int npcol2;
  int mb2;

  double* pB;
  double* pg;
  double* pB2;
  double* pg2;
  double* rr;
  double *diag;

  int descB[DLEN1_];
  int descg[DLEN1_];
  int descB2[DLEN1_];
  int descg2[DLEN1_];
  int descDy[DLEN1_];
 
  int maxmp;
  int maxnp;
  int maxmp2;
  int maxnp2;

  char UPLQ;

} scalapackM;

static struct DSDPSchurMat_Ops scalapackdsdpops;
static int DSDPScalapackOpsInit(struct DSDPSchurMat_Ops*);

extern int DSDPSetSchurMatOps(DSDP,struct DSDPSchurMat_Ops*,void*);

#undef __FUNCT__
#define __FUNCT__ "PDSDPUseSCALAPACKLinearSolver"
int PDSDPUseSCALAPACKLinearSolver(DSDP dsdp){
  int info;
  scalapackM* ctx;

  DSDPFunctionBegin;
  DSDPCALLOC1(&ctx,scalapackM,&info);
  info=DSDPScalapackOpsInit(&scalapackdsdpops);DSDPCHKERR(info);
  info=DSDPSetSchurMatOps(dsdp,&scalapackdsdpops,(void*)ctx);DSDPCHKERR(info);

  ctx->fflag=0;
  ctx->mb=1;
  ctx->mb2=1;
  ctx->mb2=40;

  ctx->UPLQ='U';
  ctx->rank=0;
  ctx->nprocs=0;
  DSDPFunctionReturn(0);
}

static int rowonp(scalapackM* ctx,int row,int *mrow){
  if (row % ctx->npcol==ctx->mycol){
    if (mrow) *mrow=(row/ctx->npcol)*ctx->maxmp;
    return 1;
  } else {
    if (mrow) *mrow=-1;
    return 0;
  }
}

#undef __FUNCT__
#define __FUNCT__ "pmatsetup"
int pmatsetup(void *MM, int m){
  scalapackM* ctx=(scalapackM*)MM;
  int maxmp,maxnp,maxmp2,maxnp2;
  int info;
  int IZERO = 0, IONE = 1;

  DSDPFunctionBegin;
  ctx->m=m;

  /*
  ctx->nprow=1;
  ctx->npcol=1;
  sl_init__(&ctx->ictxt,&ctx->nprow,&ctx->npcol);
  */

  Cblacs_pinfo(&ctx->rank, &ctx->nprocs);

  Cblacs_get(-1,0,&ctx->ictxt);

  Cblacs_gridinit(&ctx->ictxt,"Row-major",1,ctx->nprocs);
  Cblacs_gridinfo(ctx->ictxt,&ctx->nprow,&ctx->npcol,&ctx->myrow,&ctx->mycol);

  Cblacs_get(-1,0,&ctx->ictxt2);
  ctx->nprow2 = (int) sqrt((double)ctx->nprocs);
  ctx->npcol2 = ctx->nprocs/ctx->nprow2;
  while (ctx->nprow2*ctx->npcol2 != ctx->nprocs) {
    ctx->nprow2--;
    ctx->npcol2 = ctx->nprocs/ctx->nprow2;
  }
  Cblacs_gridinit(&ctx->ictxt2,"Row-major",ctx->nprow2,ctx->npcol2);
  Cblacs_gridinfo(ctx->ictxt2,&ctx->nprow2,&ctx->npcol2,&ctx->myrow2,&ctx->mycol2);
  maxmp = DSDPMax(1,numroc_(&m,&ctx->mb,&ctx->myrow,&IZERO,&ctx->nprow));
  maxnp = DSDPMax(1,numroc_(&m,&ctx->mb,&ctx->mycol,&IZERO,&ctx->npcol));
  maxmp2 = DSDPMax(1,numroc_(&m,&ctx->mb2,&ctx->myrow2,&IZERO,&ctx->nprow2));
  maxnp2 = DSDPMax(1,numroc_(&m,&ctx->mb2,&ctx->mycol2,&IZERO,&ctx->npcol2));

  ctx->maxmp=maxmp;
  ctx->maxnp=maxnp;
  ctx->maxmp2=maxmp2;
  ctx->maxnp2=maxnp2;
  if (0){
    printf("RANK: %-3d, NPROCS: %-3d, M: %-3d\n",ctx->rank,ctx->nprocs,m);
    printf("RANK: %-3d, ICTXT : %-3d,MYROW: %-3d of NPROW: %-3d, MYCOL: %-3d, NPCOL: %-3d\n",
	   ctx->rank,ctx->ictxt,ctx->myrow,ctx->nprow,ctx->mycol,ctx->npcol);
    printf("RANK: %-3d, ICTXT2: %-3d,MYROW: %-3d of NPROW: %-3d, MYCOL: %-3d, NPCOL: %-3d\n",
	   ctx->rank,ctx->ictxt2,ctx->myrow2,ctx->nprow2,ctx->mycol2,ctx->npcol2);
    
    printf("RANK: %-3d, block : %-3d, maxmp : %-3d, maxnp : %-3d\n",ctx->rank,ctx->mb,maxmp,maxnp);
    printf("RANK: %-3d, block2: %-3d, maxmp2: %-3d, maxnp2: %-3d\n",ctx->rank,ctx->mb2,maxmp2,maxnp2);
  }
  DSDPCALLOC2(&ctx->pB,double,maxnp*maxmp,&info);
  if (info){printf("OUT of Memory\n"); return 1;}

  DSDPCALLOC2(&ctx->pg,double,maxmp,&info);
  if (info){printf("OUT of Memory\n"); return 1;}

  DSDPCALLOC2(&ctx->pB2,double,maxmp2*maxnp2,&info);
  if (info){printf("OUT of Memory\n"); return 1;}

  DSDPCALLOC2(&ctx->pg2,double,maxmp2,&info);
  if (info){printf("OUT of Memory\n"); return 1;}

  DSDPCALLOC2(&ctx->rr,double,m,&info);
  if (info){printf("OUT of Memory\n"); return 1;}

  DSDPCALLOC2(&ctx->diag,double,m,&info);
  if (info){printf("OUT of Memory\n"); return 1;}

  descinit_(ctx->descB,&m,&m,&m,&ctx->mb,&IZERO,&IZERO,
	    &ctx->ictxt,&m,&info);
  if (info!=0) { printf("DESCINIT ERROR1"); return info; }

  descinit_(ctx->descg,&m,&IONE,&ctx->mb,&IONE,&IZERO,&IZERO,
	    &ctx->ictxt,&maxmp,&info);
  if (info!=0) { printf("DESCINIT ERROR2"); return info; }

  descinit_(ctx->descB2,&m,&m,&ctx->mb2,&ctx->mb2,&IZERO,&IZERO,
	    &ctx->ictxt2,&maxmp2,&info);
  if (info!=0) { printf("DESCINIT ERROR3"); return info; }

  descinit_(ctx->descg2,&m,&IONE,&ctx->mb2,&IONE,&IZERO,&IZERO,
	    &ctx->ictxt2,&maxmp2,&info);
  if (info!=0) {printf("DESCINIT ERROR4"); return info;}

  descinit_(ctx->descDy,&m,&IONE,&m,&IONE,&IZERO,&IZERO,
	    &ctx->ictxt,&m,&info);
  if (info!=0) { printf("DESCINIT ERROR5"); return info;}


  ctx->t0=0;ctx->t1=0;ctx->t2=0;
  ctx->thessian=0;ctx->tfactor=0;ctx->tsolve=0;
  DSDPTime(&ctx->t0);
  DSDPFunctionReturn(0);
}

static void pvecscalem(double v1[], double v2[], double v3[], int n){
  int i;
  for (i=0;i<n;i++){
    v3[i] = (v1[i]*v2[i]);
  }
  return;
}

static void pvecscaled(double v1[], double v2[], double v3[], int n){
  int i;
  for (i=0;i<n;i++){
    v3[i] = (v1[i]/v2[i]);
  }
  return;
}

static int pmatzero(void*MM){
  scalapackM* ctx=(scalapackM*)MM;
  int i,nn1=ctx->maxmp*ctx->maxnp, nn2=ctx->maxmp2*ctx->maxnp2;
  DSDPFunctionBegin;
  DSDPTime(&ctx->t1);
  memset((void*)ctx->pB,0,nn1*sizeof(double));
  memset((void*)ctx->pB2,0,nn2*sizeof(double));
  for (i=0;i<ctx->m;i++){ctx->diag[i]=1.0;}
  DSDPFunctionReturn(0);
}

static int pmatonprocessor( void* MM,int row ,int*yesorno){
  scalapackM* ctx=(scalapackM*)MM;
  if (rowonp(ctx,row,0)){
    *yesorno=1;
  } else {
    *yesorno=0;
  }
  return 0;
}

static int pmatlocalvariables( void* MM,double vars[] ,int m){
  scalapackM* ctx=(scalapackM*)MM;
  int row;
  for (row=0;row<m;row++){
    if (rowonp(ctx,row,0)){
      vars[row]=1;
    } else {
      vars[row]=0;
    }
  }
  return 0;
}

static int pmatmult(void *MM, double x[], double y[], int m){
  scalapackM* ctx=(scalapackM*)MM;
  int IA=1,JA=1,IX=1,JX=1,IY=1,JY=1,INCX=1,INCY=1;
  int info,N=ctx->m;
  char UPLQ=ctx->UPLQ;
  double *A=ctx->pB;
  double ALPHA=1.0,BETA=0.0;
  double *X=ctx->pg,*Y=y;

  DSDPFunctionBegin;

  memset((void*)y,0,m*sizeof(double));
  pvecscaled(x,ctx->diag,ctx->rr,m);
  Cpdgemr2d(m,1,ctx->rr,1,1,ctx->descDy,X,1,1,ctx->descg,ctx->ictxt);

  pdsymv_(&UPLQ,&N,&ALPHA,A,&IA,&JA,ctx->descB,X,&IX,&JX,ctx->descg,&INCX,&BETA,Y,&IY,&JY,ctx->descDy,&INCY);
  
  info=pvecreduce(MM,y,m);DSDPCHKERR(info);
  pvecscaled(y,ctx->diag,y,m);
  
  DSDPFunctionReturn(0);
}

static int pmatrowcolumns(void* MM,int row ,double cols[],int*ncols, int nrows){
  scalapackM* ctx=(scalapackM*)MM;
  int i;
  DSDPFunctionBegin;
  memset((void*)cols,0,nrows*sizeof(int));
  if (rowonp(ctx,row,0)){
    *ncols=row+1;
    for (i=0;i<=row;i++){ cols[i]=1.0;}
  } else {
    *ncols=0;
  }
  DSDPFunctionReturn(0);
}

static int pmatshiftdiagonal(void* MM,double shift){
  scalapackM* ctx=(scalapackM*)MM;
  int mrow,row,m=ctx->m;
  DSDPFunctionBegin;
  if (shift==0) return 0;
  for (row=0;row<m;row++){
    if (rowonp(ctx,row,&mrow)){
      ctx->pB[mrow+row]+=10*shift;
    }
  }
  DSDPFunctionReturn(0);
}

static int pmataddelement(void* MM,int row, double shift){
  scalapackM* ctx=(scalapackM*)MM;
  int mrow;
  DSDPFunctionBegin;
  if (rowonp(ctx,row,&mrow)){
      ctx->pB[mrow+row]+=shift;
  }
  DSDPFunctionReturn(0);
}


static int pmatadddiagonal(void*MM, double val[],int n){
  scalapackM* ctx=(scalapackM*)MM;
  int mrow,row;

  DSDPFunctionBegin;
  for (row=0;row<n;row++){
    if (rowonp(ctx,row,&mrow)){
      ctx->pB[mrow+row]+=val[row];
    }
  }
  DSDPFunctionReturn(0);
}

static int pmatgetdiagonal(void*MM, double val[],int n){
  scalapackM* ctx=(scalapackM*)MM;
  int info,mrow,row;

  DSDPFunctionBegin;
  for (row=0;row<n;row++){
    if (rowonp(ctx,row,&mrow)){
      val[row]=fabs(ctx->pB[mrow+row]);
      if (val[row]==0) val[row]=1.0;
    } else {
      val[row]=0;
    }
  }
  info=pvecreduce(MM,val,n);DSDPCHKERR(info);
  for (row=0;row<n;row++){
    val[row]=sqrt(1.0/val[row]);
  }
  DSDPFunctionReturn(0);
}

static int pmatscale(void*MM){
  scalapackM* ctx=(scalapackM*)MM;
  int mrow,i,row,m=ctx->m;

  DSDPFunctionBegin;
  for (row=0;row<m;row++){
    if (rowonp(ctx,row,&mrow)){
      for (i=0;i<m;i++){
	ctx->pB[mrow+i]=ctx->pB[mrow+i]*ctx->diag[row]*ctx->diag[i];
	if (i!=row && fabs(ctx->pB[mrow+i])>1.0){printf("BAD ELEMENT: %d, %d, %4.10e\n",row,i,ctx->pB[mrow+i]);}
      }
    }
  }
  DSDPFunctionReturn(0);
}



/* Cache lines before inserting them into the matrix */
static int pmataddline(void*MM, int row, double dd, double vals[], int m){
  scalapackM* ctx=(scalapackM*)MM;
  int i,mrow;

  DSDPFunctionBegin;
  if (rowonp(ctx,row,&mrow)){
    for (i=0;i<=row;i++){
      ctx->pB[mrow+i]+=dd*vals[i];
    }
  }
  DSDPFunctionReturn(0);
}

static int pmatassemble(void*MM){
  scalapackM* ctx=(scalapackM*)MM;
  int info;
  DSDPFunctionBegin;

  DSDPTime(&ctx->t2);
  ctx->thessian+=ctx->t2-ctx->t1;
  info=pmatshiftdiagonal(MM,1.0e-12);DSDPCHKERR(info);

  DSDPFunctionReturn(0);
}

static int pmatdistributed(void*MM, int*flag){
  DSDPFunctionBegin;
  *flag=1;
  DSDPFunctionReturn(0);
}

static int pmatfactor(void*MM, int *flag){
  scalapackM* ctx=(scalapackM*)MM;
  int info=0,info2=0;
  int N=ctx->m,m=ctx->m;
  int IA=1,JA=1;
  double *A=ctx->pB2;
  char UPLQ=ctx->UPLQ;
  DSDPFunctionBegin;

  info=pmatgetdiagonal(MM,ctx->diag,ctx->m);DSDPCHKERR(info);
  info=pmatscale(MM);DSDPCHKERR(info);
  info=pmatshiftdiagonal(MM,1.0e-11);DSDPCHKERR(info);
  Cpdgemr2d(m,m,ctx->pB,1,1,ctx->descB,ctx->pB2,1,1,ctx->descB2,ctx->ictxt);

  PPDSDPPrintTime(ctx->rank,"Compute M",ctx->t2-ctx->t1,ctx->thessian);
  if (0){
    int IONE=1,IZERO=0;
    int ISTDOUT = 6;
    double tmpwork[100];
    pdlaprnt_(&m,&m,ctx->pB,&IONE,&IONE,ctx->descB,&IZERO,&IZERO,
	      "B0 = ", &ISTDOUT,tmpwork,5);
    pdlaprnt_(&m,&m,ctx->pB2,&IONE,&IONE,ctx->descB2,&IZERO,&IZERO,
	      "B2 = ", &ISTDOUT,tmpwork,5);
  }

  DSDPTime(&ctx->t1);
  *flag=0;
  info=pdpotrf_(&UPLQ,&N,A,&IA,&JA,ctx->descB2,&info2,1);
  if (info2!=0) {*flag=1;}
  if (info2!=0) {ctx->fflag=1;}
  if (info2!=0) {printf("BAD FACTOR : %d\n",info2); }
  DSDPTime(&ctx->t2);
  ctx->tfactor+=ctx->t2-ctx->t1;
  PPDSDPPrintTime(ctx->rank,"SCALAPACK: Factor M",ctx->t2-ctx->t1,ctx->tfactor);
  DSDPFunctionReturn(0);
}

static int pvecreduce(void*MM,double *v, int m){
  scalapackM* ctx=(scalapackM*)MM;
  int i;
  DSDPFunctionBegin;

  for (i=0;i<m;i++){ctx->rr[i]=v[i]; v[i]=0;}
  MPI_Allreduce((void*)ctx->rr,(void*)v,m,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  /*
    char top;
  void Cdgsum2d(int ConTxt, char *scope, char *top, int m, int n, double *A,
		int lda, int rdest, int cdest);
  top = *PB_Ctop(ctx->ictxt,"COMBINE","All", TOP_GET );
  Cdgsum2d(ctx->ictxt,"All",&top,m,1,v,m,-1,-1);
  */
  DSDPFunctionReturn(0);
}

static int pmatsolve(void* MM, double bb[], double xx[], int n){
  scalapackM* ctx=(scalapackM*)MM;
  int info=0,info2=0,m=ctx->m;
  int N=ctx->m, NRHS=1;
  int IA=1,JA=1,IB=1,JB=1;
  double *A=ctx->pB2;
  char UPLQ=ctx->UPLQ;
  int i;

  DSDPFunctionBegin;
  DSDPTime(&ctx->t1);

  /* 
     Copy RHS from DSDPVector to SCALAPACK vector. This assumes the entries in the local
     DSDP Vectors are not duplicated on multiple processors. It adds the first element of
     each vector together, second element, ... 
  */

  pvecscalem(bb,ctx->diag,ctx->rr,m);
  Cpdgemr2d(m,1,ctx->rr,1,1,ctx->descDy,ctx->pg2,1,1,ctx->descg2,ctx->ictxt2);
  
  if (0){
    int IONE=1,IZERO=0;
    int ISTDOUT = 6;
    double tmpwork[100];
    pdlaprnt_(&m,&IONE,ctx->pg2,&IONE,&IONE,ctx->descg2,&IZERO,&IZERO,
	      "rhs = ", &ISTDOUT,tmpwork,5);
  }

  info=pdpotrs_(&UPLQ,&N,&NRHS,A,&IA,&JA,ctx->descB2,ctx->pg2,&IB,&JB,ctx->descg2,&info2,1);
  if (info2!=0){ printf("CHOLESKY SOLVE INFO: %d\n",info2);return 10;}

  if (0){
    int IONE=1,IZERO=0;
    int ISTDOUT = 6;
    double tmpwork[100];
    pdlaprnt_(&m,&IONE,ctx->pg2,&IONE,&IONE,ctx->descg2,&IZERO,&IZERO,
	      "dy = ", &ISTDOUT,tmpwork,5);
  }

  Cpdgemr2d(m,1,ctx->pg2,1,1,ctx->descg2,xx,1,1,ctx->descDy,ctx->ictxt2);

  info=pvecreduce(MM,xx,m);DSDPCHKERR(info);
  if (xx[0]!=xx[0]){ printf("BAD CHOLESKY SOLVE.\n");return 10;}
  pvecscalem(xx,ctx->diag,xx,m);

  if (0){
    int IONE=1,IZERO=0;
    int ISTDOUT = 6;
    double tmpwork[100];
    pdlaprnt_(&m,&IONE,xx,&IONE,&IONE,ctx->descDy,&IZERO,&IZERO,
	      "dy1 = ", &ISTDOUT,tmpwork,5);
  }

  if (0){
    double derr=0,dd3,*rr3;
    rr3=(double*)malloc(n*sizeof(double));
    info=pmatmult(MM,xx,rr3,n);DSDPCHKERR(info);
    for (i=0;i<-m;i++){ printf("RANK:%d %d RHS1: %4.4e, RHS2: %4.4e, SCL: %4.4e\n",ctx->rank,i,bb[i],rr3[i],ctx->diag[i]);}
    for (i=0;i<m;i++){ dd3=rr3[i] - bb[i]; derr+=dd3*dd3; }
    free(rr3);
    printf("RANK: %d, DY ERROR: %4.4e\n",ctx->rank,sqrt(derr));
  }

  DSDPTime(&ctx->t2);
  ctx->tsolve+=ctx->t2-ctx->t1;
  /* PPDSDPPrintTime(ctx->rank,"Solve M",ctx->t2-ctx->t1,ctx->tsolve);*/
  DSDPFunctionReturn(0);
}

static int pmatdestroy(void*MM){
  scalapackM* ctx=(scalapackM*)MM;
  int info;
  DSDPFunctionBegin;

  Cblacs_barrier(ctx->ictxt,"All");
  Cblacs_gridexit(ctx->ictxt);
  Cblacs_barrier(ctx->ictxt2,"All");
  Cblacs_gridexit(ctx->ictxt2);
  /*  pbfreebuf(); */
  /*  Cblacs_exit(0); */
  PPDSDPPrint("PDSDP: Compute M: %5.5e \n",ctx->thessian);
  PPDSDPPrint("PDSDP: Factor M %4.5e seconds \n",ctx->tfactor);
  PPDSDPPrint("PDSDP: Solve M %4.5e seconds \n",ctx->tsolve);
  DSDPFREE(&ctx->pB,&info);DSDPCHKERR(info);
  DSDPFREE(&ctx->pB2,&info);DSDPCHKERR(info);
  DSDPFREE(&ctx->pg,&info);DSDPCHKERR(info);
  DSDPFREE(&ctx->pg2,&info);DSDPCHKERR(info);
  DSDPFREE(&ctx->rr,&info);DSDPCHKERR(info);
  DSDPFREE(&ctx->diag,&info);DSDPCHKERR(info);
  DSDPFREE(&ctx,&info);DSDPCHKERR(info);
  DSDPFunctionReturn(0);
}


static const char *scalapackmatname="SCALAPACK matrix";

static int DSDPScalapackOpsInit(struct DSDPSchurMat_Ops* sops){
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
  sops->pmatreduction=pvecreduce;
  sops->matsetup=pmatsetup;
  sops->matdestroy=pmatdestroy;
  sops->id=25;
  sops->matname=scalapackmatname;
  return 0;
}
