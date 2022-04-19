#include "../../src/solver/dsdpschurmat_impl.h"
#include "pdsdp5petsc.h"

/* The Schur Complement matrix has particular operations, and several implementations */
/* PDSDP  simply provides a different implementation of these operations.             */


typedef struct {
  MPI_Comm comm;
  int rank;
  int m, mlocal;
  int low,high;
  Mat M;
  Vec X,RHS,Vall;
  VecScatter scatter;
  KSP ksp;
  int *iptr;
  DSDP dsdp;
  PetscLogDouble t0,t1,thessian,tsolve;
  int hevent,dsdpevent;
  int hmethod;
} petscM;


static int pmatzero(void*MM){
  petscM* MC=(petscM*)MM;
  Mat M= MC->M;
  int i,m=MC->m,info;
  PetscFunctionBegin;
  info = PetscGetTime(&MC->t1);
  info=PetscLogEventBegin(MC->hevent,0,0,0,0);
  info=MatZeroEntries(M);CHKERRQ(info);
  for (i=0;i<m;i++) MC->iptr[i]=i;
  PetscFunctionReturn(0);  
}


static int onrow(petscM* MC,int row){
  if (MC->low<= row && row < MC->high){
    return 1;
  }
  return 0;
}

static int pmatlocalvariables( void* MM,double vars[] ,int m){
  petscM* MC=(petscM*)MM;
  int row;
  for (row=0;row<m;row++){
    if (onrow(MC,row)){
      vars[row]=1;
    } else {
      vars[row]=0;
    }
  }
  return 0;
}

static int pmatonprocessor( void* MM,int row ,int*yesorno){
  petscM* MC=(petscM*)MM;
  *yesorno=onrow(MC,row);
  return 0;
}
 
static int pmatidcolumns(void* MM, int row, double cols[], int *ncols, int m){

  petscM* MC=(petscM*)MM;
  int i,count;
  PetscFunctionBegin;

  if (MC->low<= row && row < MC->high){
    switch (MC->hmethod){
    case 2: 
      *ncols = m;
      for (i=0;i<m;i++){ cols[i]=1.0; }
      break;
    case 3:
      count=0;
      for (i=(row+1)%2 ;i<row; i+=2){
	cols[i]=1.0; count++;
      }
      for (i=row;i<m;i+=2){
	cols[i]=1.0; count++;
      }
      *ncols=count;
      break;
    default:
      *ncols = m-row;
      for (i=row;i<m;i++){ cols[i]=1.0; }
    } 
  } else {
    *ncols=0;
  }
  PetscFunctionReturn(0);  
}


static int pmataddline(void*MM, int row, double dd, double vals[],int m){
  petscM* MC=(petscM*)MM;
  int i,info;
  PetscFunctionBegin;
  for (i=0;i<m;i++)vals[i]*=dd;

  switch (MC->hmethod){
  case 2:
    info=MatSetValues(MC->M,1,&row,m,MC->iptr,vals,ADD_VALUES);CHKERRQ(info);
    break;
  case 3:
    info=MatSetValues(MC->M,1,&row,m,MC->iptr,vals,ADD_VALUES);CHKERRQ(info);
    vals[row]=0;
    info=MatSetValues(MC->M,m,MC->iptr,1,&row,vals,ADD_VALUES);CHKERRQ(info);
    break;
  default:
    info=MatSetValues(MC->M,1,&row,m-row,MC->iptr+row,vals+row,ADD_VALUES);CHKERRQ(info);
    info=MatSetValues(MC->M,m-row-1,MC->iptr+row+1,1,&row,vals+row+1,ADD_VALUES);CHKERRQ(info);
  }
  PetscFunctionReturn(0);  
}

static int pmatassemble(void*MM){
  petscM* MC=(petscM*)MM;
  Mat M= MC->M;
  PetscLogDouble t2;
  int info;
  PetscFunctionBegin;
  info=MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);CHKERRQ(info);
  info=MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);CHKERRQ(info);
  info=PetscLogEventEnd(MC->hevent,0,0,0,0);
  info = PetscGetTime(&t2);
  MC->thessian+=(t2-MC->t1);
  PetscLogInfo((M,"Compute Hessian: %10.6e, Total: %10.6e\n",t2-MC->t1,MC->thessian));
  PetscFunctionReturn(0);  
}

static int pmatfactor(void*MM,int*flag){
  petscM* MC=(petscM*)MM;
  int info;
  PetscFunctionBegin;
  info=KSPSetOperators(MC->ksp,MC->M,MC->M,SAME_NONZERO_PATTERN);CHKERRQ(info);
  *flag=0;
  PetscFunctionReturn(0);  
}

static int pmatreduce(void*MM,double *v, int m){
  petscM* MC=(petscM*)MM;
  Vec VP=MC->RHS,Vall=MC->Vall;
  int i,info;
  double *vv, zero=0.0;

  PetscFunctionBegin;
  info=VecGetArray(Vall,&vv);CHKERRQ(info);
  for (i=0;i<m;i++){ vv[i]=v[i]; }
  info=VecRestoreArray(Vall,&vv);CHKERRQ(info);

  info = VecSet(VP,zero);CHKERRQ(info);
  info = VecScatterBegin(Vall,VP,ADD_VALUES,SCATTER_REVERSE,MC->scatter);CHKERRQ(info);
  info = VecScatterEnd(Vall,VP,ADD_VALUES,SCATTER_REVERSE,MC->scatter);CHKERRQ(info);

  /* Copy PETSc Vector into DSDP vector */
  info = VecSet(Vall,zero);CHKERRQ(info);
  info = VecScatterBegin(VP,Vall,INSERT_VALUES,SCATTER_FORWARD,MC->scatter);CHKERRQ(info);
  info = VecScatterEnd(VP,Vall,INSERT_VALUES,SCATTER_FORWARD,MC->scatter);CHKERRQ(info);
  /*
  info=VecView(Vall,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(info);
  */
  info=VecGetArray(Vall,&vv);CHKERRQ(info);
  for (i=0;i<m;i++){ v[i]=vv[i]; }
  info=VecRestoreArray(Vall,&vv);CHKERRQ(info);
  PetscFunctionReturn(0);  
}

static int pmatsolve(void* MM, double bb[],double xx[], int m){
  petscM* MC=(petscM*)MM;
  Vec DX=MC->X,PRHS=MC->RHS, Vall=MC->Vall;
  int i,its,info,low,high;
  double *vv;
  double zero=0.0;
  KSP ksp=MC->ksp;
  KSPConvergedReason reason;
  PetscReal rnorm,rnorm2,rtol=1.0e-6,aatol=1.0e-12;  
  PetscLogDouble t0,t1;

  PetscFunctionBegin;
  info = PetscGetTime(&t0);

  info=VecSet(PRHS,zero);CHKERRQ(info);
  info=VecGetOwnershipRange(PRHS,&low,&high);CHKERRQ(info);
  info=VecGetArray(PRHS,&vv);CHKERRQ(info);
  for (i=low;i<high;i++){ vv[i-low]=bb[i]; }
  info=VecRestoreArray(PRHS,&vv);CHKERRQ(info);

  /* Solve the linear system */
  info=VecSet(DX,zero);CHKERRQ(info);

  info = KSPSetInitialGuessNonzero(ksp,PETSC_FALSE); CHKERRQ(info);
  info = KSPSetTolerances(ksp,rtol,aatol,1e+25,m/3+20);CHKERRQ(info);
  info = KSPSolve(ksp,PRHS,DX);CHKERRQ(info);
  info = KSPGetIterationNumber(ksp, &its); CHKERRQ(info);
  info = KSPGetConvergedReason(ksp,&reason); CHKERRQ(info);
  info = KSPGetResidualNorm(ksp,&rnorm); CHKERRQ(info);
  PetscLogInfo((MC->M,"DSDP : First Try Linear iterations: %d, ",its));
  PetscLogInfo((MC->M,"DSDP : Residual Norm: %4.2e, ",rnorm));
  PetscLogInfo((MC->M,"DSDP : Termination Reason: %d\n",reason));
  if (reason!=2 && reason!=3){
    info=KSPSetInitialGuessNonzero(ksp,PETSC_TRUE); CHKERRQ(info);
    info = KSPSetTolerances(ksp,rtol/rnorm,aatol,1e+25,m/3+20);CHKERRQ(info);
    info = KSPSolve(ksp,PRHS,DX);CHKERRQ(info);
    info = KSPGetIterationNumber(ksp, &its); CHKERRQ(info);
    info=KSPGetConvergedReason(ksp,&reason); CHKERRQ(info);
    info = KSPGetResidualNorm(ksp,&rnorm2); CHKERRQ(info);
    PetscLogInfo((MC->M,"DSDP : More Linear iterations: %d, ",its));
    PetscLogInfo((MC->M,"DSDP : Residual Norm: %4.2e, ",rnorm2));
    PetscLogInfo((MC->M,"DSDP : Termination Reason: %d\n",reason));
    if (reason!=2 && reason!=3 && rnorm2<rnorm/10){
      info=KSPSetInitialGuessNonzero(ksp,PETSC_TRUE); CHKERRQ(info);
      info = KSPSetTolerances(ksp,rtol/(rnorm*rnorm2),aatol,1e+25,m/3+20);CHKERRQ(info);
      info = KSPSolve(ksp,PRHS,DX);CHKERRQ(info);
      info = KSPGetIterationNumber(ksp, &its); CHKERRQ(info);
      info = KSPGetResidualNorm(ksp,&rnorm); CHKERRQ(info);
      info=KSPGetConvergedReason(ksp,&reason); CHKERRQ(info);
      PetscLogInfo((MC->M,"DSDP : STill More Linear iterations: %d, ",its));
      PetscLogInfo((MC->M,"DSDP : Residual Norm: %4.2e, ",rnorm));
      PetscLogInfo((MC->M,"DSDP : Termination Reason: %d\n",reason));
    }
  }

  /* Copy PETSc Vector into DSDP vector */
  info = VecSet(Vall,zero);CHKERRQ(info);
  info = VecScatterBegin(DX,Vall,INSERT_VALUES,SCATTER_FORWARD,MC->scatter);CHKERRQ(info);
  info = VecScatterEnd(DX,Vall,INSERT_VALUES,SCATTER_FORWARD,MC->scatter);CHKERRQ(info);
  /*
  info=VecView(Vall,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(info);
  */
  info=VecGetArray(Vall,&vv);CHKERRQ(info);
  for (i=0;i<m;i++){  xx[i]=vv[i]; }
  info=VecRestoreArray(Vall,&vv);CHKERRQ(info);
  info = PetscGetTime(&t1);
  MC->tsolve+=(t1-t0);

  PetscLogInfo((MC->M,"DSDP : Solve Linear system: %10.6e, Total: %10.6e\n",t1-t0,MC->tsolve));
  PetscFunctionReturn(0);  
}

static int pmatdestroy(void*MM){
  petscM* MC=(petscM*)MM;
  int info;
  PetscLogDouble tf;
  PetscFunctionBegin;
  info = PetscGetTime(&tf);
  PetscLogInfo((MC->M,"DSDP : DSDP Total: %10.6e\n",tf-MC->t0));
  info=MatDestroy(MC->M);CHKERRQ(info);
  info=VecScatterDestroy(MC->scatter);CHKERRQ(info);
  info=VecDestroy(MC->Vall);CHKERRQ(info);
  info=VecDestroy(MC->X);CHKERRQ(info);
  info=VecDestroy(MC->RHS);CHKERRQ(info);
  info=KSPDestroy(MC->ksp);CHKERRQ(info);
  info=PetscLogEventEnd(MC->dsdpevent,0,0,0,0);
  PetscFree(MC->iptr);
  /*  PetscFinalize(); */
  PetscFunctionReturn(0);CHKERRQ(info);
}

static int pmatview(void*MM){
  petscM* MC=(petscM*)MM;
  int info;
  PetscFunctionBegin;
  info=KSPView(MC->ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(info);
  info=MatView(MC->M,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(info);
  PetscFunctionReturn(0);  
}


#undef __FUNCT__
#define __FUNCT__ "DSDPLinearSolverCreate"
static int DSDPPetscLinearSolverCreate(void *MM, int m){
  int i,info,mlocal;
  PC pc;
  IS ISLocal,ISGlobal;
  petscM* MC=(petscM*)MM;

  PetscFunctionBegin;

  info = VecCreateSeq(PETSC_COMM_SELF,m,&MC->Vall);CHKERRQ(info);
  info = VecCreateMPI(MC->comm,MC->mlocal,m,&MC->X);CHKERRQ(info);
  info = VecDuplicate(MC->X,&MC->RHS);CHKERRQ(info);

  /* Create the parallel matrix */
  info = VecGetOwnershipRange(MC->X,&MC->low,&MC->high);CHKERRQ(info);
  info = VecGetLocalSize(MC->X,&mlocal);CHKERRQ(info);
  info = MatCreateMPIDense(MC->comm,mlocal,mlocal,m,m,0,&MC->M);CHKERRQ(info);

  info = MatSetOption(MC->M,MAT_COLUMNS_SORTED); CHKERRQ(info);
  info = MatSetOption(MC->M,MAT_ROWS_SORTED); CHKERRQ(info);

  MC->m=m; MC->mlocal=mlocal;

  /* Create the scatter between local and parallel vectors */
  info = ISCreateStride(PETSC_COMM_SELF,m,0,1,&ISLocal);CHKERRQ(info);
  info = ISCreateStride(MC->comm,m,0,1,&ISGlobal); CHKERRQ(info);
  info = VecScatterCreate(MC->X,ISGlobal,MC->Vall,ISLocal,&MC->scatter); CHKERRQ(info);
  info = ISDestroy(ISLocal);CHKERRQ(info);
  info = ISDestroy(ISGlobal);CHKERRQ(info);

  info = PetscMalloc(m*sizeof(int),&MC->iptr);
  for (i=0;i<m;i++) MC->iptr[i]=i;

  /* Create the Linear Solver */
  info = KSPCreate(MC->comm,&MC->ksp);CHKERRQ(info);
  info = KSPSetTolerances(MC->ksp,1.0e-8,1.0e-13,1e+25,1000);CHKERRQ(info);
  info = KSPGetPC(MC->ksp,&pc); CHKERRQ(info);
  info = KSPSetType(MC->ksp,KSPCG); CHKERRQ(info);
  info = PCSetType(pc,PCJACOBI); CHKERRQ(info);
  info = KSPSetFromOptions(MC->ksp); CHKERRQ(info);
  MC->hmethod=3;
  MC->t0=0; MC->t1=0; MC->thessian=0; MC->tsolve=0;
  info=PetscLogEventRegister(&MC->dsdpevent,"DSDPSolve",KSP_COOKIE);
  info=PetscLogEventRegister(&MC->hevent,"DSPPComputeM",KSP_COOKIE);

  info=PetscLogEventBegin(MC->dsdpevent,0,0,0,0);

  info = PetscGetTime(&MC->t0);
  
  PetscFunctionReturn(0);
}


static int TTTMatMult(void*MM,double xx[],double yy[],int m){
  petscM* MC=(petscM*)MM;
  Vec X=MC->X,Y=MC->RHS, Vall=MC->Vall;
  double *x,*y,*vv,zero=0.0;
  int i,info,low,high;
  PetscFunctionBegin;

  info=VecSet(X,zero);CHKERRQ(info);
  info=VecGetOwnershipRange(X,&low,&high);CHKERRQ(info);
  info=VecGetArray(X,&vv);CHKERRQ(info);
  for (i=low;i<high;i++){ vv[i-low]=xx[i]; }
  info=VecRestoreArray(X,&vv);CHKERRQ(info);

  info=MatMult(MC->M,X,Y);CHKERRQ(info);

  info = VecScatterBegin(Y,Vall,INSERT_VALUES,SCATTER_FORWARD,MC->scatter);CHKERRQ(info);
  info = VecScatterEnd(Y,Vall,INSERT_VALUES,SCATTER_FORWARD,MC->scatter);CHKERRQ(info);
  info=VecGetArray(Vall,&y);CHKERRQ(info);
  for (i=0;i<m;i++){  yy[i]=y[i]; }
  info=VecRestoreArray(Vall,&y);CHKERRQ(info);
  PetscFunctionReturn(0);
}

static int TTTMatShiftDiagonal(void *MM, double d){
  petscM* MC=(petscM*)MM;
  int info;
  PetscFunctionBegin;
  info=MatShift(MC->M,d);CHKERRQ(info);
  PetscFunctionReturn(0);
}

static int TTTMatAddDiagonal(void *MM, double diag[], int m){
  petscM* MC=(petscM*)MM;
  Vec D=MC->X;
  int i,info;
  double *d, zero=0.0;
  PetscFunctionBegin;
  info=VecSet(D,zero);CHKERRQ(info);
  info=VecGetArray(D,&d);CHKERRQ(info);
  for (i=MC->low;i<MC->high;i++) d[i-MC->low]=diag[i];
  info=VecRestoreArray(D,&d);CHKERRQ(info);
  info=MatDiagonalSet(MC->M,D,ADD_VALUES);CHKERRQ(info);
  PetscFunctionReturn(0);
}

static int pmataddelement(void* MM,int row, double shift){
  petscM* MC=(petscM*)MM;
  int info;
  PetscFunctionBegin;
  if (onrow(MC,row)){
    info=MatSetValues(MC->M,1,&row,1,&row,&shift,ADD_VALUES);CHKERRQ(info);
  }
  PetscFunctionReturn(0);
}

static int pmatdistributed(void*MM, int*flag){
  PetscFunctionBegin;
  *flag=1;
  PetscFunctionReturn(0);
}

static struct DSDPSchurMat_Ops petscdsdpops;
static const char* tmatname="PETSc";

static int DSDPPetscOpsInit(struct DSDPSchurMat_Ops* sops){
  int info;
  if (!sops) return 0;
  info=DSDPSchurMatOpsInitialize(sops); CHKERRQ(info);
  sops->matzero=pmatzero;
  sops->mataddrow=pmataddline;
  sops->matdestroy=pmatdestroy;
  sops->matfactor=pmatfactor;
  sops->matsolve=pmatsolve;
  sops->matrownonzeros=pmatidcolumns;
  sops->matassemble=pmatassemble;
  sops->matsetup=DSDPPetscLinearSolverCreate;
  sops->matadddiagonal=TTTMatAddDiagonal;
  sops->matshiftdiagonal=TTTMatShiftDiagonal;
  sops->matscaledmultiply=TTTMatMult;
  sops->matview=pmatview;
  sops->mataddelement=pmataddelement;
  sops->pmatwhichdiag=0;
  sops->pmatonprocessor=pmatonprocessor;
  sops->pmatlocalvariables=pmatlocalvariables;
  sops->pmatdistributed=pmatdistributed;
  sops->pmatreduction=pmatreduce;

  sops->id=22;
  sops->matname=tmatname;
  return 0;
}

int DSDPSetSchurMatOps(DSDP,struct DSDPSchurMat_Ops*,void*);

#undef __FUNCT__
#define __FUNCT__ "PDSDPUsePETScLinearSolver"
int PDSDPUsePETScLinearSolver(DSDP dsdp, MPI_Comm comm, int mlocal){
  int info; 
  petscM* MM;
  PetscFunctionBegin;
  /*  PetscInitialize(0,0,0,0); */
  info=PetscMalloc(sizeof(petscM),&MM);
  MM->dsdp=dsdp;
  MM->comm=comm;
  MM->mlocal=mlocal;
  MPI_Comm_rank(comm,&MM->rank);
  info=DSDPPetscOpsInit(&petscdsdpops);CHKERRQ(info);
  info=DSDPSetSchurMatOps(dsdp,&petscdsdpops,(void*)MM);CHKERRQ(info);

  PetscFunctionReturn(0);
}
