
#define PRINT_LOAD_BALANCE 0
  // for debug use

#include <pthread.h>
#include <sched.h>

#include "sdpa_newton.h"
#include "sdpa_parts.h"
#include "sdpa_jordan.h"
#include "sdpa_linear.h"  
#include "sdpa_algebra.h"
#include "sdpa_mpicopy.h"
#include "sdpa_scalapack.h"

  // from sdpa_rpdpotrf.cpp
extern "C" {
  int rpdpotrf_(char* uplo, int* n, double* A,
		int* ia, int* ja, int* desca, int* info);
};

namespace sdpa {

pthread_mutex_t Newton::job_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  Newton::job_cond  = PTHREAD_COND_INITIALIZER;
int  Newton::Column_Number = 0;
bool Newton::mutex_flag    = false;
int  Newton::Calc_F1       = 0;

Newton::Newton()
{
  useFormula = NULL;
  bMat_type  = DENSE;
  m = 0;
  // Caution: if SDPA doesn't use sparse bMat, 
  //          following variables are indefinite.
  this->SDP_nBlock = -1;
  SDP_number = NULL;  SDP_location_sparse_bMat = NULL;
  SDP_constraint1 = NULL;  SDP_constraint2 = NULL;
  SDP_blockIndex1 = NULL;  SDP_blockIndex2 = NULL;
  this->SOCP_nBlock = -1;
  SOCP_number = NULL;  SOCP_location_sparse_bMat = NULL;
  SOCP_constraint1 = NULL;  SOCP_constraint2 = NULL;
  SOCP_blockIndex1 = NULL;  SOCP_blockIndex2 = NULL;
  this->LP_nBlock = -1;
  LP_number = NULL;  LP_location_sparse_bMat = NULL;
  LP_constraint1 = NULL;  LP_constraint2 = NULL;
  LP_blockIndex1 = NULL;  LP_blockIndex2 = NULL;

  // For Parallel Distributed Sparse Schur Matrix
  mySchurStart    = -1;
  mySchurEnd      = -1;
  mySchurLength   = 0;

  diagonalIndex = NULL;

  // parallel start
  pB  = NULL;
  pg  = NULL;
  pB2 = NULL;
  pg2 = NULL;
  maxnp  = 0;
  maxmp2 = 0;
  maxnp2 = 0;
  symmetric_b = false;
  
  NUM_THREADS  = 1;
  NUM_GOTOBLAS = 1;
  threadSchurStart = NULL;
  threadSchurEnd   = NULL;
  threadSchurLoad  = NULL;
  threadSDPLocStart = NULL;
  threadSDPLocEnd   = NULL;
  threadLPLocStart = NULL;
  threadLPLocEnd   = NULL;
}

Newton::Newton(int m, BlockStruct& bs)
{
  initialize(m, bs);
}

Newton::~Newton()
{
  terminate();
}

void Newton::initialize(int m, BlockStruct& bs)
{
  this->m = m;
  
  // gVec.initialize(m);
  // for parallel only sparse gVec is used
  
  SDP_nBlock  = bs.SDP_nBlock;
  SOCP_nBlock = bs.SOCP_nBlock;
  LP_nBlock   = bs.LP_nBlock;
  
  DxMat.initialize(bs);
  DyVec.initialize(m);
  DzMat.initialize(bs);
  r_zinvMat.initialize(bs);
  x_rd_zinvMat.initialize(bs);

  // Memory allocation of useFormula is moved to computeFormula_SDP
  // NewArray(useFormula,FormulaType,m*SDP_nBlock);


  bMat_type = DENSE;
  // Caution: if SDPA doesn't use sparse bMat, 
  //          following variables are indefinite.
  this->SDP_nBlock = -1;
  SDP_number = NULL;  SDP_location_sparse_bMat = NULL;
  SDP_constraint1 = NULL;  SDP_constraint2 = NULL;
  SDP_blockIndex1 = NULL;  SDP_blockIndex2 = NULL;
  this->SOCP_nBlock = -1;
  SOCP_number = NULL;  SOCP_location_sparse_bMat = NULL;
  SOCP_constraint1 = NULL;  SOCP_constraint2 = NULL;
  SOCP_blockIndex1 = NULL;  SOCP_blockIndex2 = NULL;
  this->LP_nBlock = -1;
  LP_number = NULL;  LP_location_sparse_bMat = NULL;
  LP_constraint1 = NULL;  LP_constraint2 = NULL;
  LP_blockIndex1 = NULL;  LP_blockIndex2 = NULL;

  // For Parallel Distributed Sparse Schur Matrix
  mySchurStart    = -1;
  mySchurEnd      = -1;
  mySchurLength   = 0;
  threadSchurStart  = NULL;
  threadSchurEnd    = NULL;
  threadSchurLoad   = NULL;
  threadSDPLocStart = NULL;
  threadSDPLocEnd   = NULL;
  threadLPLocStart = NULL;
  threadLPLocEnd   = NULL;
  
  diagonalIndex = NULL;

  // parallel start
  pB  = NULL;
  pg  = NULL;
  pB2 = NULL;
  pg2 = NULL;
  maxnp  = 0;
  maxmp2 = 0;
  maxnp2 = 0;
  symmetric_b = false;
  
}

void Newton::terminate()
{

  if (bMat_type == SPARSE){

    if (SDP_location_sparse_bMat && SDP_constraint1 && SDP_constraint2
	&& SDP_blockIndex1 && SDP_blockIndex2) {
      for (int l=0; l<SDP_nBlock; ++l) {
	DeleteArray(SDP_location_sparse_bMat[l]);
	DeleteArray(SDP_constraint1[l]);
	DeleteArray(SDP_constraint2[l]);
	DeleteArray(SDP_blockIndex1[l]);
	DeleteArray(SDP_blockIndex2[l]);
      }
      DeleteArray(SDP_number);
      DeleteArray(SDP_location_sparse_bMat);
      DeleteArray(SDP_constraint1);
      DeleteArray(SDP_constraint2);
      DeleteArray(SDP_blockIndex1);
      DeleteArray(SDP_blockIndex2);
      for (int nThreads = 0; nThreads < NUM_THREADS; ++nThreads) {
	DeleteArray(threadSDPLocStart[nThreads]);
	DeleteArray(threadSDPLocEnd[nThreads]);
      }
      DeleteArray(threadSDPLocStart);
      DeleteArray(threadSDPLocEnd);
    }
#if 0
    if (SOCP_location_sparse_bMat && SOCP_constraint1 && SOCP_constraint2
	&& SOCP_blockIndex1 && SOCP_blockIndex2) {
      for (int l=0; l<SOCP_nBlock; ++l) {
	DeleteArray(SOCP_location_sparse_bMat[l]);
	DeleteArray(SOCP_constraint1[l]);
	DeleteArray(SOCP_constraint2[l]);
	DeleteArray(SOCP_blockIndex1[l]);
	DeleteArray(SOCP_blockIndex2[l]);
      }
      DeleteArray(SOCP_number);
      DeleteArray(SOCP_location_sparse_bMat);
      DeleteArray(SOCP_constraint1);
      DeleteArray(SOCP_constraint2);
      DeleteArray(SOCP_blockIndex1);
      DeleteArray(SOCP_blockIndex2);
    }
#endif
    if (LP_location_sparse_bMat && LP_constraint1 && LP_constraint2
	&& LP_blockIndex1 && LP_blockIndex2) {
      for (int l=0; l<LP_nBlock; ++l) {
	DeleteArray(LP_location_sparse_bMat[l]);
	DeleteArray(LP_constraint1[l]);
	DeleteArray(LP_constraint2[l]);
	DeleteArray(LP_blockIndex1[l]);
	DeleteArray(LP_blockIndex2[l]);
      }
      DeleteArray(LP_number);
      DeleteArray(LP_location_sparse_bMat);
      DeleteArray(LP_constraint1);
      DeleteArray(LP_constraint2);
      DeleteArray(LP_blockIndex1);
      DeleteArray(LP_blockIndex2);
      for (int nThreads = 0; nThreads < NUM_THREADS; ++nThreads) {
	DeleteArray(threadLPLocStart[nThreads]);
	DeleteArray(threadLPLocEnd[nThreads]);
      }
      DeleteArray(threadLPLocStart);
      DeleteArray(threadLPLocEnd);
    }

    DeleteArray(diagonalIndex);
    sparse_bMat.terminate();

    DeleteArray(threadSchurStart);
    DeleteArray(threadSchurEnd);
    DeleteArray(threadSchurLoad);

  } else { // bMat_type == DENSE
    // bMat.terminate();
    DeleteArray(pB);
    DeleteArray(pg);
    DeleteArray(pB2);
  }

  DeleteArray(pg2);
  gVec.terminate();
  DxMat.terminate();
  DyVec.terminate();
  DzMat.terminate();
  r_zinvMat.terminate();
  x_rd_zinvMat.terminate();

  if (useFormula) {
    for (int j=0; j<m; ++j) {
      DeleteArray(useFormula[j]);
    }
    DeleteArray(useFormula);
  }
}

void Newton::initialize_bMat_parallel(int m)
{
  maxnp = numroc_(&m,&IONE,&MpiSt::mycol,
		  &IZERO,&MpiSt::npcol);
  if (maxnp < 1) {
    maxnp = 1;
  }
  maxmp2 = numroc_(&m,&MpiSt::mb2,&MpiSt::myrow2,
		   &IZERO,&MpiSt::nprow2);
  if (maxmp2 < 1) {
    maxmp2 = 1;
  }
  maxnp2 = numroc_(&m,&MpiSt::mb2,&MpiSt::mycol2,
		   &IZERO,&MpiSt::npcol2);
  if (maxnp2 < 1) {
    maxnp2 = 1;
  }

  int info = 0;
  descinit_(descB,&m,&m,&m,&IONE,&IZERO,&IZERO,
            &MpiSt::ictxt,&m,&info);
  if (info!=0 && MpiSt::iam == 0) {
    rMessage("info " << info);
  }
  descinit_(descg,&m,&IONE,&m,&IONE,&IZERO,&IZERO,
            &MpiSt::ictxt,&m,&info);
  if (info!=0 && MpiSt::iam == 0) {
    rMessage("info " << info);
  }
  descinit_(descB2,&m,&m,&MpiSt::mb2,&MpiSt::mb2,&IZERO,&IZERO,
            &MpiSt::ictxt2,&maxmp2,&info);
  if (info!=0 && MpiSt::iam == 0) {
    rMessage("info " << info);
  }
  descinit_(descg2,&m,&IONE,&MpiSt::mb2,&IONE,&IZERO,&IZERO,
            &MpiSt::ictxt2,&maxmp2,&info);
  if (info!=0 && MpiSt::iam == 0) {
    rMessage("info " << info);
  }
  descinit_(descDy,&m,&IONE,&m,&IONE,&IZERO,&IZERO,
            &MpiSt::ictxt,&m,&info);
  if (info!=0 && MpiSt::iam == 0) {
    rMessage("info " << info);
  }
  // parallel end
  // No use of pg in this version;
  // NewArray(pg , double, maxmp *1);
  if (MpiSt::mycol2 == 0) {
    NewArray(pg2, double, maxmp2*1);
  }
}

void Newton::initialize_dense_bMat(int m)
{
  //  bMat_type = DENSE;
  //  printf("DENSE computations\n");
  // bMat.initialize(m,m,DenseMatrix::DENSE);

  NewArray(pB , double, m     *maxnp);
  NewArray(pB2, double, maxmp2*maxnp2);
}

  // 2008/03/12 kazuhide nakata
void Newton::initialize_sparse_bMat(int m)
{

  //  bMat_type = SPARSE;
  //  printf("SPARSE computation\n");

  // for parallel, gVec is used only for sparse case
  // gVec.initialize(m);
  // copy pg to Dy directly
  
  // initialize sparse_bMat by Chordal::makeGraph
  //  sparse_bMat.display();

  // make index of diagonalIndex
  NewArray(diagonalIndex,int,m+1);
  int k=0;
  for (int index=0; index<sparse_bMat.NonZeroCount; index++){
    if (sparse_bMat.row_index[index] == sparse_bMat.column_index[index]) {
      diagonalIndex[k] = index;
      k++;
    }
  }
  diagonalIndex[m] = sparse_bMat.NonZeroCount;
  #if 0
  rMessage("diagonalIndex = ");
  for (int index=0; index <m; ++index) {
    printf(" [%d:%d]",index,diagonalIndex[index]);
  }
  printf("\n");
  #endif
}

  // 2008/03/12 kazuhide nakata
void Newton::initialize_bMat(int m, Chordal& chordal,
			     InputData& inputData,
			     FILE* Display,
                             FILE* fpOut)
{
  /* Create clique tree */
  initialize_bMat_parallel(m);
  switch (chordal.best) {
  case SELECT_DENSE: 
    bMat_type = DENSE;
    if (Display) {
      fprintf(Display,"Schur computation : DENSE \n");
    }
    if (fpOut) {
      fprintf(fpOut,"Schur computation : DENSE \n");
    }
    initialize_dense_bMat(m);
    // Here, we release MUMPS and sparse_bMat
    chordal.terminate();
    break;
  case SELECT_MUMPS_BEST: 
    bMat_type = SPARSE;
    if (Display) {
      fprintf(Display,"Schur computation : SPARSE \n");
    }
    if (fpOut) {
      fprintf(fpOut,"Schur computation : SPARSE \n");
    }
    initialize_sparse_bMat(m);
    make_aggrigateIndex(inputData);
    break;
  default: 
    rError("Wrong Ordering Obtained");
    break;
  }

}


int Newton::binarySearchIndex(int i, int j)
{
  // binary search for index of sparse_bMat 
  int t = -1;
  // We store only the lower triangular
  int ii = i, jj = j;
  if (i<j) {
    jj = i;
    ii = j;
  }
  int begin  = diagonalIndex[jj]; 
  int end    = diagonalIndex[jj+1]-1;
  int target = (begin + end) / 2;
  while (end - begin > 1){
    if (sparse_bMat.row_index[target] < ii+1){
      begin = target;
      target = (begin + end) / 2;
    } else if (sparse_bMat.row_index[target] > ii+1){
      end = target;
      target = (begin + end) / 2;
    } else if (sparse_bMat.row_index[target] == ii+1) {
      t = target;
      break;
    }
  }
  if (t == -1){
    if (sparse_bMat.row_index[begin] == ii+1){
      t = begin;
    } else if (sparse_bMat.row_index[end] == ii+1){
      t = end;
    } else {
      rError("Newton::make_aggrigateIndex program bug");
    }
  } 
  return t;
}
  

void Newton::make_aggrigateIndex_SDP(InputData& inputData)
{
  SDP_nBlock = inputData.SDP_nBlock;
  NewArray(SDP_number,int,SDP_nBlock);

  // memory allocate for aggrigateIndex
  NewArray(SDP_constraint1,int*,SDP_nBlock);
  NewArray(SDP_constraint2,int*,SDP_nBlock);
  NewArray(SDP_blockIndex1,int*,SDP_nBlock);
  NewArray(SDP_blockIndex2,int*,SDP_nBlock);
  NewArray(SDP_location_sparse_bMat,int*,SDP_nBlock);

  for (int l=0; l<SDP_nBlock; l++){
    const int size = (inputData.SDP_nConstraint[l] + 1) 
      * inputData.SDP_nConstraint[l] / 2;
    SDP_number[l] = size;
    NewArray(SDP_constraint1[l],int,size);
    NewArray(SDP_constraint2[l],int,size);
    NewArray(SDP_blockIndex1[l],int,size);
    NewArray(SDP_blockIndex2[l],int,size);
    NewArray(SDP_location_sparse_bMat[l],int,size);
  }

  for (int l = 0; l<SDP_nBlock; l++){
    int NonZeroCount = 0;

    for (int k1=0; k1<inputData.SDP_nConstraint[l]; k1++){
      int j = inputData.SDP_constraint[l][k1];
      int jb = inputData.SDP_blockIndex[l][k1];
      int jnz = inputData.A[j].SDP_sp_block[jb].NonZeroEffect;

      for (int k2=0; k2<inputData.SDP_nConstraint[l]; k2++){
	int i = inputData.SDP_constraint[l][k2];
	int ib = inputData.SDP_blockIndex[l][k2];
	int inz = inputData.A[i].SDP_sp_block[ib].NonZeroEffect;

	if ((jnz < inz) || ((inz == jnz) && (i < j))){
	  continue;
	}
	// set index which A_i and A_j are not zero matrix
	SDP_constraint1[l][NonZeroCount] = i;
	SDP_constraint2[l][NonZeroCount] = j;
	SDP_blockIndex1[l][NonZeroCount] = ib;
	SDP_blockIndex2[l][NonZeroCount] = jb;
	int target = binarySearchIndex(i,j);
	SDP_location_sparse_bMat[l][NonZeroCount] = target;
	NonZeroCount++;
      }
    } // for k1
  } //for l  lth block
}


void Newton::make_aggrigateIndex_SOCP(InputData& inputData)
{
  SOCP_nBlock = inputData.SOCP_nBlock;
  NewArray(SOCP_number,int,SOCP_nBlock);
  if (SOCP_number == NULL) {
    rError("Newton::make_aggrigateIndex_SOCP memory exhausted ");
  }

  NewArray(SOCP_constraint1,int*,SOCP_nBlock);
  NewArray(SOCP_constraint2,int*,SOCP_nBlock);
  NewArray(SOCP_blockIndex1,int*,SOCP_nBlock);
  NewArray(SOCP_blockIndex2,int*,SOCP_nBlock);
  NewArray(SOCP_location_sparse_bMat,int*,SOCP_nBlock);

  for (int l=0; l<SOCP_nBlock; l++){
    int size = (inputData.SOCP_nConstraint[l]+1)
      * inputData.SOCP_nConstraint[l]/2;
    SOCP_number[l] = size;
    NewArray(SOCP_constraint1[l],int,size);
    NewArray(SOCP_constraint2[l],int,size);
    NewArray(SOCP_blockIndex1[l],int,size);
    NewArray(SOCP_blockIndex2[l],int,size);
    NewArray(SOCP_location_sparse_bMat[l],int,size);
  }

  for (int l = 0; l<SOCP_nBlock; l++){
    int NonZeroCount = 0;

    for (int k1=0; k1<inputData.SOCP_nConstraint[l]; k1++){
      int j = inputData.SOCP_constraint[l][k1];
      int jb = inputData.SOCP_blockIndex[l][k1];
      int jnz = inputData.A[j].SOCP_sp_block[jb].NonZeroEffect;

      for (int k2=0; k2<inputData.SOCP_nConstraint[l]; k2++){
	int i = inputData.SOCP_constraint[l][k2];
	int ib = inputData.SOCP_blockIndex[l][k2];
	int inz = inputData.A[i].SOCP_sp_block[ib].NonZeroEffect;

	if ((jnz < inz) || ((inz == jnz) && (i < j))){
	  continue;
	}

	// set index which A_i and A_j are not zero matrix
	SOCP_constraint1[l][NonZeroCount] = i;
	SOCP_constraint2[l][NonZeroCount] = j;
	SOCP_blockIndex1[l][NonZeroCount] = ib;
	SOCP_blockIndex2[l][NonZeroCount] = jb;

	int target = binarySearchIndex(i,j);
	SOCP_location_sparse_bMat[l][NonZeroCount] = target;
	NonZeroCount++;
      }
    } // for k1
  } //for l  lth block
}

void Newton::make_aggrigateIndex_LP(InputData& inputData)
{
  LP_nBlock = inputData.LP_nBlock;
  NewArray(LP_number,int,LP_nBlock);

  // memory allocate for aggrigateIndex
  NewArray(LP_constraint1,int*,LP_nBlock);
  NewArray(LP_constraint2,int*,LP_nBlock);
  NewArray(LP_blockIndex1,int*,LP_nBlock);
  NewArray(LP_blockIndex2,int*,LP_nBlock);
  NewArray(LP_location_sparse_bMat,int*,LP_nBlock);

  for (int l=0; l<LP_nBlock; l++){
    int size = (inputData.LP_nConstraint[l]+1)
      * inputData.LP_nConstraint[l]/2;
    LP_number[l] = size;
    NewArray(LP_constraint1[l],int,size);
    NewArray(LP_constraint2[l],int,size);
    NewArray(LP_blockIndex1[l],int,size);
    NewArray(LP_blockIndex2[l],int,size);
    NewArray(LP_location_sparse_bMat[l],int,size);
  }

  for (int l = 0; l<LP_nBlock; l++){
    int NonZeroCount = 0;

    for (int k1=0; k1<inputData.LP_nConstraint[l]; k1++){
      int j = inputData.LP_constraint[l][k1];
      int jb = inputData.LP_blockIndex[l][k1];

      for (int k2=0; k2<inputData.LP_nConstraint[l]; k2++){
	int i = inputData.LP_constraint[l][k2];
	int ib = inputData.LP_blockIndex[l][k2];

	if (i < j){
	  continue;
	}

	// set index which A_i and A_j are not zero matrix
	LP_constraint1[l][NonZeroCount] = i;
	LP_constraint2[l][NonZeroCount] = j;
	LP_blockIndex1[l][NonZeroCount] = ib;
	LP_blockIndex2[l][NonZeroCount] = jb;
	
	int target = binarySearchIndex(i,j);
	LP_location_sparse_bMat[l][NonZeroCount] = target;
	NonZeroCount++;
      }
    } // for k1
  } //for l  lth block
}

void Newton::make_aggrigateIndex(InputData& inputData)
{
  make_aggrigateIndex_SDP(inputData);
  //  make_aggrigateIndex_SOCP(inputData);
  make_aggrigateIndex_LP(inputData);
}

void Newton::computeFormula_SDP(InputData& inputData,
				double DenseRatio, double Kappa)
{
  int m = inputData.b.nDim;
  int SDP_nBlock = inputData.SDP_nBlock;

  NewArray(useFormula, FormulaType*, m);
  for (int j=0; j<m; ++j) {
    NewArray(useFormula[j], FormulaType, inputData.A[j].SDP_sp_nBlock);
  }
  // need not to initialize useFormula
  
  int** upNonZeroCount;
  NewArray(upNonZeroCount, int*, m);
  for (int j=0; j<m; ++j) {
    NewArray(upNonZeroCount[j], int, inputData.A[j].SDP_sp_nBlock);
  }
  for (int j=0; j<m; ++j) {
    for (int jb=0; jb < inputData.A[j].SDP_sp_nBlock; ++jb) {
      upNonZeroCount[j][jb] = 0;
    }
  }

  SparseLinearSpace* A = inputData.A;

  #if 0
  for (int k=0; k<m; ++k) {
    for (int l=0; l<inputData.A[0].nBlock; ++l) {
      rMessage("A[" << k << "].ele[" << l << "] ="
	       << inputData.A[k].ele[l].NonZeroEffect);
    }
  }
  #endif

  // Count sum of number of elements
  // that each number of elements are less than own.

  for (int l=0; l<SDP_nBlock; ++l) {
    for (int k1=0; k1 < inputData.SDP_nConstraint[l];k1++){
      int j = inputData.SDP_constraint[l][k1];
      int jb = inputData.SDP_blockIndex[l][k1];
      int jnz = A[j].SDP_sp_block[jb].NonZeroEffect;
      int up = jnz;
      // rMessage("up = " << up);

      for (int k2=0; k2 < inputData.SDP_nConstraint[l];k2++){
	int i = inputData.SDP_constraint[l][k2];
	int ib = inputData.SDP_blockIndex[l][k2];
	int inz = A[i].SDP_sp_block[ib].NonZeroEffect;
	//	printf("%d %d %d %d %d %d\n",i,ib,inz, j, jb,jnz);
	if (inz < jnz) {
	  up += inz;
	}
#if 1
	else if ((jnz == inz) && i>j ) {
	  up += inz;
	}
#endif
      }
      upNonZeroCount[j][jb] = up;
      // rMessage("up = " << up);
    }
  }

  // Determine which formula
  double ff=0, ff1=0, ff2=0, ff3=0;
  Calc_F1 = 0;
  for (int l=0; l<SDP_nBlock; ++l) {
    int countf1,countf2,countf3;
    countf1 = countf2 = countf3 = 0;
    for (int k=0; k < inputData.SDP_nConstraint[l]; k++){
      int j =  inputData.SDP_constraint[l][k];
      int jb =  inputData.SDP_blockIndex[l][k];
      double jnz = inputData.A[j].SDP_sp_block[jb].NonZeroEffect;

      double f1,f2,f3;
      double n       = inputData.A[j].SDP_sp_block[jb].nRow;
      double up      = upNonZeroCount[j][jb];

#if 1
      f1 = Kappa*n*jnz + n*n*n + Kappa*up;
      f2 = Kappa*n*jnz + Kappa*(n+1)*up;
      #if 0
      f3 = Kappa*(2*Kappa*jnz+1)*up/Kappa;
      #else
      f3 = Kappa*(2*Kappa*jnz+1)*up;
      #endif

#endif
      ff1 += f1;
      ff2 += f2;
      ff3 += f3;
      //rMessage("up = " << up << " nonzero = " << nonzero);
      //rMessage("f1=" << f1 << " f2=" << f2 << " f3=" << f3);
      //printf("%d %d %lf %lf %lf\n",k,l,f1,f2,f3);
      //printf("%d %d %le %le %le\n",k,l,f1,f2,f3);
      if (inputData.A[j].SDP_sp_block[jb].type == SparseMatrix::DENSE) {
	// if DENSE, we use only F1 or F2,
	// that is we don't use F3
	if (f1<f2) {
	  useFormula[j][jb] = F1;
	  countf1++;
	  ff += f1;
	} else {
	  useFormula[j][jb] = F2;
	  countf2++;
	  ff += f2;
	}
      } else {
	//printf("n = %d, inz = %d, n * n = %d, inz / (n * n) = %lf\n", (int)n, (int)inz, (int)(n*n), (double)inz/(n*n));
	// this case is SPARSE
	if (f1<f2 && f1<f3) {
	  if ((n <= 200) && (2.0 * n >= jnz)) {
	    //rMessage("line " << k << " is F3:1");
	    useFormula[j][jb] = F3;
	    countf3++;
	    ff += f3;
	  }
	  else {
	    //rMessage("line " << k << " is F1");
	    useFormula[j][jb] = F1;
	    countf1++;
	    ff += f1;
	  }
	} else if (f2<f3) {
	  //rMessage("line " << k << " is F2");
	  useFormula[j][jb] = F2;
	  countf2++;
	  ff += f2;
	} else {
	  //rMessage("line " << k << " is F3:2");
	  useFormula[j][jb] = F3;
	  countf3++;
	  ff += f3;
	}
      }
    }

    Calc_F1 += countf1;
    // rMessage("Kappa = " << Kappa);
    #if 0
    rMessage("count f1 = " << countf1
	     << ":: count f2 = " << countf2
	     << ":: count f3 = " << countf3);
    #endif

    // if we need to symmetrize Schur complement matrix,
    // that is, we compute both of uppper and lower part
    // and adjust to symmetric matrix,
    // we change the variable 'symmetric_b' to true.
    // This effect for only DENSE case

    if (bMat_type == DENSE && countf1 !=m && countf2 !=m
	&& (countf1 >= 1 || countf2 >= 2) ) {
      // rMessage("need symmetric");
      symmetric_b = true;
    }

    
  } // end of 'for (int l)'

  for (int j=0; j<m; ++j) {
    DeleteArray(upNonZeroCount[j]);
  }
  DeleteArray(upNonZeroCount);

  return;
}

void Newton::accumulateFormulaCost(InputData& inputData,
				   double Kappa)
{
  // Acculated cost is temporaraly stored into 
  // the elements of Schur complement
  // set sparse_bMat zero
  sdpa_dset(sparse_bMat.NonZeroCount, 0.0, sparse_bMat.sp_ele, 1);
  accumulateFormulaCost_SDP(inputData,Kappa);
  // accumulateFormulaCost_SOCP(inputData,Kappa);
  accumulateFormulaCost_LP(inputData,Kappa);

  MpiSt::barrier();
}

void Newton::accumulateFormulaCost_SDP(InputData& inputData,
				       double Kappa)
{
  int SDP_nBlock = inputData.SDP_nBlock;

  SparseLinearSpace* A = inputData.A;
  for (int l=0; l<SDP_nBlock; ++l) {
    int previous_j = -1;
    for (int iter = 0; iter < SDP_number[l]; ++iter) {
      double currVal = 0.0;
      int j = SDP_constraint2[l][iter];
      int jb = SDP_blockIndex2[l][iter];
      double jnz = inputData.A[j].SDP_sp_block[jb].NonZeroEffect;
      double n = inputData.A[j].SDP_sp_block[jb].nRow;
      FormulaType formula = useFormula[j][jb];
      if (j != previous_j) {
	if (formula==F1) {
	  currVal += n*n*n + Kappa*n*jnz;
	}
	if (formula==F2) {
	  currVal += Kappa*n*jnz;
	}
	previous_j = j;
      }
      
      int i = SDP_constraint1[l][iter];
      int ib = SDP_blockIndex1[l][iter];
      double inz = A[i].SDP_sp_block[ib].NonZeroEffect;
      
      if (formula==F3) {
	#if 0
	currVal += (2*Kappa*jnz+1)*inz/Kappa;
	#else
	currVal += (2*Kappa*jnz+1)*inz*1.0;
	#endif
      }
      else if (formula==F2) {
	currVal += (n+1)*Kappa*inz;
      }
      else if (formula==F1) {
	currVal += Kappa*inz;
      }
      sparse_bMat.sp_ele[SDP_location_sparse_bMat[l][iter]] += currVal;
      #if 0
      rMessage("Location = " << SDP_location_sparse_bMat[l][iter]
	       << " => (" << i << " , " << j << ")"
	       << " at iter " << iter);
      #endif
    }
  } // end of 'for (int l)'
  return;
}

void Newton::accumulateFormulaCost_SOCP(InputData& inputData,
					double Kappa)
{
  // What Shall I do ?
}

void Newton::accumulateFormulaCost_LP(InputData& inputData,
				      double Kappa)
{
  double tau = 1.0;
  // tau decides the ratio between SDP cost and LP cost.
  // Large tau means LP cost is maginified.
  // So, if the bottleneck thread incurres heavier LP cost,
  //   change tau to smaller values.
  for (int l=0; l<LP_nBlock; ++l) {
    for (int iter = 0; iter < LP_number[l]; iter++){
      double value = tau * Kappa;
      sparse_bMat.sp_ele[LP_location_sparse_bMat[l][iter]] += value;
    } // end of 'for (int iter)
  } // end of 'for (int l)'
  return;
}

void Newton::computeSchurIndices()
{
  int start, end;
  double load = 0;
  int CPU     = MpiSt::nprocs;
  int id      = MpiSt::iam;
  int N       = sparse_bMat.NonZeroCount;
  double* array = sparse_bMat.sp_ele; 

  // #define array2(x,i,CPUindex) (x * pow((N-i+1)/(double)N,2.0) )
  #define array2(x,i,CPUindex) (x)
  double sum = 0.0;
  for (int i=0; i<N; ++i) {
    sum += array2(array[i],i,CPUindex);
  }
  double ave = sum / CPU;
  #if 0
  if (MpiSt::iam == 0) {
    rMessage(" sum = " << sum << " : ave = " << ave);
  }
  #endif
  
  // if too small, all elements are computed in id==0
  if (2*CPU > N) {
    if (id == 0) {
      start = 0;
      end = N;
      load = sum;
    }
    else {
      start = end = N;
      load = 0.0;
    }
  }
  else { // We distribute the Schur to multiple matrices
    int index = 0;
    start = end = 0;
    double currentCPULoad = 0.0;
    for (int currentCPU = 0; currentCPU <= id; ++currentCPU) {
      currentCPULoad = 0.0;
      start = end;
      while (index < N) {
	if (currentCPULoad
	    + array2(array[index],index,currentCPU)/2 > ave) {
	  end = index;
	  break;
	}
	currentCPULoad += array2(array[index],index,
				 currentCPU);
	index++;
      }
      if (index == N) {
	end = N;
      }
    }
    if (id == CPU-1 && end < N) {
      for (int index2=end; index2<N; ++index2) {
	currentCPULoad += array2(array[index2],index,CPU-1);
      }
      end = N;
    }
    load = currentCPULoad;
  }

  mySchurStart  = start;
  mySchurEnd    = end;
  mySchurLength = end - start;

  // The same scheme for threads
  
  NewArray(threadSchurStart, int, NUM_THREADS); 
  NewArray(threadSchurEnd,   int, NUM_THREADS);
  NewArray(threadSchurLoad,  double, NUM_THREADS);

  NewArray(threadSDPLocStart, int*, NUM_THREADS);
  NewArray(threadSDPLocEnd,   int*, NUM_THREADS);
  NewArray(threadLPLocStart, int*, NUM_THREADS);
  NewArray(threadLPLocEnd,   int*, NUM_THREADS);
  for (int nThreads=0; nThreads<NUM_THREADS; ++nThreads) {
    NewArray(threadSDPLocStart[nThreads], int, SDP_nBlock);
    NewArray(threadSDPLocEnd  [nThreads], int, SDP_nBlock);
    NewArray(threadLPLocStart[nThreads], int, LP_nBlock);
    NewArray(threadLPLocEnd  [nThreads], int, LP_nBlock);
  }

  ave = load / NUM_THREADS;

  // if too small, all elements are computed in id==0
  if (2*NUM_THREADS > mySchurLength) {
    threadSchurStart[0] = mySchurStart;
    threadSchurEnd  [0] = mySchurEnd;
    threadSchurLoad [0] = load;
    for (int nThreads=1; nThreads < NUM_THREADS; ++nThreads) {
      threadSchurStart[nThreads] = mySchurEnd;
      threadSchurEnd  [nThreads] = mySchurEnd;
      threadSchurLoad [nThreads] = 0.0;
    }
  }
  else { // We distribute the Schur to multiple matrices
    int index = mySchurStart;
    start = end = mySchurStart;
    double currentThreadLoad = 0.0;
    for (int currentThread = 0; currentThread < NUM_THREADS;
	 ++currentThread) {
      currentThreadLoad = 0.0;
      start = end;
      while (index < mySchurEnd) {
	if (currentThreadLoad
	    + array2(array[index],index,currentThead)/2 > ave) {
	  end = index;
	  break;
	}
	currentThreadLoad += array2(array[index],index,
				    currentThread);
	index++;
      }
      if (index == mySchurEnd) {
	end = mySchurEnd;
      }
      threadSchurStart[currentThread] = start;
      threadSchurEnd  [currentThread] = end;
      threadSchurLoad [currentThread] = currentThreadLoad;
    }
    if (end < mySchurEnd) {
      for (int index2=end; index2<mySchurEnd; ++index2) {
	currentThreadLoad += array2(array[index2],index,CPU-1);
      }
      threadSchurEnd [NUM_THREADS-1] = mySchurEnd;
      threadSchurLoad[NUM_THREADS-1] = currentThreadLoad;
    }
  }

  // clean up the Schur complement matrix
  sdpa_dset(sparse_bMat.NonZeroCount, 0.0, sparse_bMat.sp_ele, 1);

  #if PRINT_LOAD_BALANCE
  MpiSt::barrier();
  if (MpiSt::iam == 0) {
    rMessage("Schur NonZero = " << N);
  }
  for (int idtmp = 0; idtmp < CPU; ++idtmp) {
    MpiSt::barrier();
    if (MpiSt::iam == idtmp) {
      rMessage("mySchur is [" << mySchurStart
	       << "," << mySchurEnd << "] : Length = "
	       << mySchurLength << " : Load = " << load);
      for (int nThreads=0; nThreads < NUM_THREADS; ++nThreads) {
	rMessage("Thread(" << nThreads << " on "
		 << MpiSt::iam << ") :: " 
		 << " [" << threadSchurStart[nThreads]
		 << "," << threadSchurEnd[nThreads]
		 << "] : Length = "
		 << (threadSchurEnd[nThreads]
		     - threadSchurStart[nThreads])
		 << " : Load = "
		 << threadSchurLoad[nThreads]);
      }
    }
    MpiSt::barrier();
  }
  MpiSt::barrier();
  #endif
}

void Newton::slimUpAggrigateIndex()
{
  // We pickup only indices whose target is between
  // mySchurStart and mySchurEnd
  
  slimUpAggrigateIndex_SDP();
  // slimUpAggrigateIndex_SOCP();
  slimUpAggrigateIndex_LP();


  // Display slimed AggrigateIndex
  #if 0
  for (int id=0; id<MpiSt::nprocs; ++id) {
    if (id == MpiSt::iam) {
      rMessage("Indices on CPU " << MpiSt::nprocs);
      display_index();
    }
    MpiSt::barrier();
  }
  #endif
}

void Newton::slimUpAggrigateIndex_SDP()
{
  if (mySchurStart == mySchurEnd) {
    // This processor do not compute Schur elements
    for (int l=0; l<SDP_nBlock; ++l) {
      SDP_number[l] = 0;
    }
    return;
  }
    
  for (int l=0; l<SDP_nBlock; ++l) {
    int searchIndex  = 0;
    // This step can be slightly faster
    // because 'iter' between start and end are successive
    for (int iter=0; iter < SDP_number[l]; ++iter) {
      if (mySchurStart <= SDP_location_sparse_bMat[l][iter]
	  && SDP_location_sparse_bMat[l][iter] < mySchurEnd) {
	SDP_constraint1[l][searchIndex] = SDP_constraint1[l][iter];
	SDP_constraint2[l][searchIndex] = SDP_constraint2[l][iter];
	SDP_blockIndex1[l][searchIndex] = SDP_blockIndex1[l][iter];
	SDP_blockIndex2[l][searchIndex] = SDP_blockIndex2[l][iter];
	SDP_location_sparse_bMat[l][searchIndex]
	  = SDP_location_sparse_bMat[l][iter];
	searchIndex++;
      }
    }
    SDP_number[l] = searchIndex;
  }
  for (int l=0; l<SDP_nBlock; ++l) {
    int searchIndex = 0;
    int end = SDP_number[l];
    // rMessage("SDP_number["<<l<<"] = " << SDP_number[l]);
    for(int nThreads = 0; nThreads < NUM_THREADS; ++nThreads) {
      threadSDPLocStart[nThreads][l] = searchIndex;
      while (searchIndex < end
	     && SDP_location_sparse_bMat[l][searchIndex]
	     < threadSchurEnd[nThreads]) {
	searchIndex++;
      }
      threadSDPLocEnd[nThreads][l] = searchIndex;
    }
  }
  MpiSt::barrier();
  #if 0
  for (int l=0; l<SDP_nBlock; ++l) {
    for (int iamtmp=0; iamtmp<MpiSt::nprocs; ++iamtmp) {
      MpiSt::barrier();
      if (MpiSt::iam == iamtmp) {
	for(int nThreads = 0; nThreads < NUM_THREADS; ++nThreads) {
	  rMessage("SDPLoc bl("<<l<<") Thread ["
		   << nThreads <<"/" << MpiSt::iam
		   << "] = [" << threadSDPLocStart[nThreads][l]
		   << " , " << threadSDPLocEnd[nThreads][l]
		   << "]");
	}
	fflush(stdout);
      }
      MpiSt::barrier();
    }
  }
  MpiSt::barrier();
  #endif
}

void Newton::slimUpAggrigateIndex_SOCP()
{
  if (mySchurStart == mySchurEnd) {
    // This processor do not compute Schur elements
    for (int l=0; l<SOCP_nBlock; ++l) {
      SOCP_number[l] = 0;
    }
    return;
  }
    
  for (int l=0; l<SOCP_nBlock; ++l) {
    int searchIndex  = 0;
    for (int iter=0; iter < SOCP_number[l]; ++iter) {
      if (mySchurStart <= SOCP_location_sparse_bMat[l][iter]
	  && SOCP_location_sparse_bMat[l][iter] < mySchurEnd) {
	SOCP_constraint1[l][searchIndex] = SOCP_constraint1[l][iter];
	SOCP_constraint2[l][searchIndex] = SOCP_constraint2[l][iter];
	SOCP_blockIndex1[l][searchIndex] = SOCP_blockIndex1[l][iter];
	SOCP_blockIndex2[l][searchIndex] = SOCP_blockIndex2[l][iter];
	SOCP_location_sparse_bMat[l][searchIndex]
	  = SOCP_location_sparse_bMat[l][iter];
	searchIndex++;
      }
    }
    SOCP_number[l] = searchIndex;
  }
}

void Newton::slimUpAggrigateIndex_LP()
{
  if (mySchurStart == mySchurEnd) {
    // This processor do not compute Schur elements
    for (int l=0; l<LP_nBlock; ++l) {
      LP_number[l] = 0;
    }
    return;
  }
    
  for (int l=0; l<LP_nBlock; ++l) {
    int searchIndex  = 0;
    for (int iter=0; iter < LP_number[l]; ++iter) {
      if (mySchurStart <= LP_location_sparse_bMat[l][iter]
	  && LP_location_sparse_bMat[l][iter] < mySchurEnd) {
	LP_constraint1[l][searchIndex] = LP_constraint1[l][iter];
	LP_constraint2[l][searchIndex] = LP_constraint2[l][iter];
	LP_blockIndex1[l][searchIndex] = LP_blockIndex1[l][iter];
	LP_blockIndex2[l][searchIndex] = LP_blockIndex2[l][iter];
	LP_location_sparse_bMat[l][searchIndex]
	  = LP_location_sparse_bMat[l][iter];
	searchIndex++;
      }
    }
    LP_number[l] = searchIndex;
  }
  for (int l=0; l<LP_nBlock; ++l) {
    int searchIndex = 0;
    int end = LP_number[l];
    // rMessage("LP_number["<<l<<"] = " << LP_number[l]);
    for(int nThreads = 0; nThreads < NUM_THREADS; ++nThreads) {
      threadLPLocStart[nThreads][l] = searchIndex;
      while (searchIndex < end
	     && LP_location_sparse_bMat[l][searchIndex]
	     < threadSchurEnd[nThreads]) {
	searchIndex++;
      }
      threadLPLocEnd[nThreads][l] = searchIndex;
    }
  }
}


void Newton::compute_rMat(Newton::WHICH_DIRECTION direction,
			  AverageComplementarity& mu,
			  DirectionParameter& beta,
			  Solutions& currentPt,
			  WorkVariables& work)
{

  //     CORRECTOR ::  r_zinv = (-XZ -dXdZ + mu I)Z^{-1}
  // not CORRECTOR ::  r_zinv = (-XZ + mu I)Z^{-1}
  double target = beta.value*mu.current;
  Lal::let(r_zinvMat,'=',currentPt.invzMat,'*',&target);
  Lal::let(r_zinvMat,'=',r_zinvMat,'+',currentPt.xMat,&DMONE);

  if (direction == CORRECTOR) {
    // work.DLS1 = Dx Dz Z^{-1}
    Jal::ns_jordan_triple_product(work.DLS1,DxMat,DzMat,
				  currentPt.invzMat,work.DLS2);
    Lal::let(r_zinvMat,'=',r_zinvMat,'+',work.DLS1,&DMONE);
  }

  //  rMessage("r_zinvMat = ");
  //  r_zinvMat.display();
}

void Newton::Make_gVec(Newton::WHICH_DIRECTION direction,
		       InputData& inputData,
		       Solutions& currentPt,
		       Residuals& currentRes,
		       AverageComplementarity& mu,
		       DirectionParameter& beta,
		       Phase& phase,
		       WorkVariables& work,
		       ComputeTime& com)
{
  TimeStart(START1);
  // rMessage("mu = " << mu.current);
  // rMessage("beta = " << beta.value);
  compute_rMat(direction,mu,beta,currentPt,work);

  TimeEnd(END1);
  com.makerMat += TimeCal(START1,END1);

  TimeStart(START2);
  TimeStart(START_GVEC_MUL);

  // work.DLS1 = R Z^{-1} - X D Z^{-1} = r_zinv - X D Z^{-1}
  if (phase.value == SolveInfo:: pFEAS
      || phase.value == SolveInfo::noINFO) {

    if (direction == CORRECTOR) {
      // x_rd_zinvMat is computed in PREDICTOR step
      Lal::let(work.DLS1,'=',r_zinvMat,'+',x_rd_zinvMat,&DMONE);
    } else {
      // currentPt is infeasilbe, that is the residual
      // dualMat is not 0.
      //      x_rd_zinvMat = X D Z^{-1}
      Jal::ns_jordan_triple_product(x_rd_zinvMat,currentPt.xMat,
				    currentRes.dualMat,currentPt.invzMat,
				    work.DLS2);
      Lal::let(work.DLS1,'=',r_zinvMat,'+',x_rd_zinvMat,&DMONE);
    } // if (direction == CORRECTOR)

  } else {
    // dualMat == 0
    work.DLS1.copyFrom(r_zinvMat);
  }
  
  //  rMessage("work.DLS1");
  //  work.DLS1.display();

  TimeEnd(END_GVEC_MUL);
  com.makegVecMul += TimeCal(START_GVEC_MUL,END_GVEC_MUL);

  if (MpiSt::mycol2 == 0) {
    // -----------------------------------------------
    int PP = m/MpiSt::mb2;
    int QQ = PP / MpiSt::nprow2;
    if (MpiSt::myrow2 < PP % MpiSt::nprow2) {
      QQ ++;
    }
    for (int pp = 0; pp < QQ ; ++pp) {
      for (int qq = 0; qq < MpiSt::mb2; ++qq) {
	int k_inter = pp*MpiSt::mb2 + qq;
	int k_outer = (pp*MpiSt::nprow2 + MpiSt::myrow2)
	  * MpiSt::mb2 + qq;
	// rMessage(" k_inter = " << k_inter << " k_outer = " << k_outer);
	double ip;
	Lal::let(ip,'=',inputData.A[k_outer],'.',work.DLS1);
	ip = -ip + currentRes.primalVec.ele[k_outer];
	pg2[k_inter] = ip;
      }
    }
    if (MpiSt::myrow2 == PP % MpiSt::nprow2) {
      int amari = m % MpiSt::mb2;
      if (amari != 0) {
	for (int qq = 0; qq < amari; ++qq) {
	  int k_inter = QQ*MpiSt::mb2 + qq;
	  int k_outer = (QQ*MpiSt::nprow2 + MpiSt::myrow2)
	    * MpiSt::mb2 + qq;
	  // rMessage(" k_inter = " << k_inter << " k_outer = " << k_outer);
	  double ip;
	  Lal::let(ip,'=',inputData.A[k_outer],'.',work.DLS1);
	  ip = -ip + currentRes.primalVec.ele[k_outer];
	  pg2[k_inter] = ip;
	}
      }
    }
  }
  // -----------------------------------------------
  
  TimeEnd(END2);
  com.makegVec += TimeCal(START2,END2);

  TimeStart(START_COPYGVEC);
  if (bMat_type == DENSE) {
    // Nothing Paticular
    // display_pg2();
  }
  else {
    // bMat_type == SPARSE
    // copy pg2 to dy of zero-th processor
    Cpdgemr2d(m,1,pg2,1,1,descg2,DyVec.ele,1,1,descDy,MpiSt::ictxt);
  }
  TimeStart(END_COPYGVEC);
  com.copygVec += TimeCal(START_COPYGVEC,END_COPYGVEC);
    
}

void Newton::calF1(double& ret, DenseMatrix& G,
		    SparseMatrix& Ai)
{
  Lal::let(ret,'=',Ai,'.',G);
}

void Newton::calF2(double& ret,
		   DenseMatrix& F, DenseMatrix& G,
		   DenseMatrix& invZ, SparseMatrix& Ai,
		    bool& hasF2Gcal)
{
  int alpha,beta;
  double value1,value2;

  int n    = Ai.nRow;
  // rMessage(" using F2 ");
  switch (Ai.type) {
  case SparseMatrix::SPARSE:
    // rMessage("F2::SPARSE  " << Ai.NonZeroCount);
    ret = 0.0;
    for (int index = 0; index < Ai.NonZeroCount; ++index) {
#if DATA_CAPSULE
      alpha  = Ai.DataS[index].vRow;
      beta   = Ai.DataS[index].vCol;
      value1 = Ai.DataS[index].vEle;
#else
      alpha  = Ai.row_index[index];
      beta   = Ai.column_index[index];
      value1 = Ai.sp_ele[index];
#endif      
      value2 = ddot_f77(&n, &F.de_ele[alpha+n*0], &n,
			&invZ.de_ele[0+n*beta], &IONE);
      ret += value1*value2;
      if (alpha!=beta) {
	value2 = ddot_f77(&n, &F.de_ele[beta+n*0], &n,
			  &invZ.de_ele[0+n*alpha], &IONE);
	ret += value1*value2;
      }
    }
    break;
  case SparseMatrix::DENSE:
    // G is temporary matrix
    // rMessage("F2::DENSE");
    if (hasF2Gcal == false) {
      // rMessage(" using F2 changing to F1");
      Lal::let(G,'=',F,'*',invZ);
      hasF2Gcal = true;
    }
    Lal::let(ret,'=',Ai,'.',G);
    break;
  } // end of switch
}

inline void Newton::calF3(double& ret,
		    DenseMatrix& X, DenseMatrix& invZ,
		    SparseMatrix& Ai, SparseMatrix& Aj)
{
  // Ai and Aj are SPARSE
  ret = 0.0;
  double sum;
  const int nCol = X.nCol;
  // rMessage("Aj.NonZeroCount = " << Aj.NonZeroCount);
  for (int index1=0; index1<Aj.NonZeroCount; ++index1) {
#if DATA_CAPSULE
    int alpha = Aj.DataS[index1].vRow;
    int beta  = Aj.DataS[index1].vCol;
    double value1 = Aj.DataS[index1].vEle;
#else
    int alpha = Aj.row_index[index1];
    int beta  = Aj.column_index[index1];
    double value1 = Aj.sp_ele[index1];
#endif
    sum = 0.0;
    for (int index2=0; index2<Ai.NonZeroCount; ++index2) {
#if DATA_CAPSULE
      int gamma = Ai.DataS[index2].vRow;
      int delta  = Ai.DataS[index2].vCol;
      double value2 = Ai.DataS[index2].vEle;
#else
      int gamma = Ai.row_index[index2];
      int delta  = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
#endif
      double plu = value2*invZ.de_ele[delta+nCol*beta]
        * X.de_ele[gamma+nCol*alpha];
      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2*invZ.de_ele[gamma+nCol*beta]
          * X.de_ele[delta+nCol*alpha];
        sum += plu2;
      }
    }
    ret += value1*sum;
    if (alpha==beta) {
      continue;
    }
    sum = 0.0;
    for (int index2=0; index2<Ai.NonZeroCount; ++index2) {
#if DATA_CAPSULE
      int gamma = Ai.DataS[index2].vRow;
      int delta  = Ai.DataS[index2].vCol;
      double value2 = Ai.DataS[index2].vEle;
#else
      int gamma = Ai.row_index[index2];
      int delta  = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
#endif
      double plu = value2*invZ.de_ele[delta+nCol*alpha]
        * X.de_ele[gamma+nCol*beta];
      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2*invZ.de_ele[gamma+nCol*alpha]
          * X.de_ele[delta+nCol*beta];
        sum += plu2;
      }
    }
    ret += value1*sum;
  } // end of 'for (index1)'
  return;
}

inline void Newton::calF1_thread(double& ret, DenseMatrix& G,
				 SparseMatrix& Ai)
{
  Lal::let(ret,'=',Ai,'.',G);
}

void Newton::calF2_thread(double& ret,
			  DenseMatrix& F, DenseMatrix& G,
			  DenseMatrix& invZ, SparseMatrix& Ai,
			  bool& hasF2Gcal)
{
  int alpha,beta;
  double value1,value2;

  int n    = Ai.nRow;
  // rMessage(" using F2 ");
  switch (Ai.type) {
  case SparseMatrix::SPARSE:
    // rMessage("F2::SPARSE  " << Aj.NonZeroCount);
    ret = 0.0;
    for (int index = 0; index < Ai.NonZeroCount; ++index) {
#if DATA_CAPSULE
      alpha  = Ai.DataS[index].vRow;
      beta   = Ai.DataS[index].vCol;
      value1 = Ai.DataS[index].vEle;
#else
      alpha  = Ai.row_index[index];
      beta   = Ai.column_index[index];
      value1 = Ai.sp_ele[index];
#endif
      value2 = ddot_f77(&n, &F.de_ele[alpha+n*0], &n,
			&invZ.de_ele[0+n*beta], &IONE);
      ret += value1*value2;
      if (alpha!=beta) {
	value2 = ddot_f77(&n, &F.de_ele[beta+n*0], &n,
			  &invZ.de_ele[0+n*alpha], &IONE);
	ret += value1*value2;
      }
    }
    break;
  case SparseMatrix::DENSE:
    // G is temporary matrix
    // rMessage("F2::DENSE");
    if (hasF2Gcal == false) {
      // rMessage(" using F2 changing to F1");
      Lal::let(G,'=',F,'*',invZ);
      hasF2Gcal = true;
    }
    Lal::let(ret,'=',Ai,'.',G);
    break;
  } // end of switch
}

#if 0
inline void Newton::calF3_thread(double& ret,
				 DenseMatrix& X, DenseMatrix& invZ,
				 SparseMatrix& Ai, SparseMatrix& Aj)
{
  // Ai and Aj are SPARSE
  ret = 0.0;
  double sum;
  const int nCol = X.nCol;
  // rMessage("Aj.NonZeroCount = " << Aj.NonZeroCount);
  for (int index1=0; index1<Aj.NonZeroCount; ++index1) {
#if DATA_CAPSULE
    int alpha = Aj.DataS[index1].vRow;
    int beta  = Aj.DataS[index1].vCol;
    double value1 = Aj.DataS[index1].vEle;
#else
    int alpha = Aj.row_index[index1];
    int beta  = Aj.column_index[index1];
    double value1 = Aj.sp_ele[index1];
#endif
    sum = 0.0;
    
    const int nCol = X.nCol;
    double *Xalpha, *Xbeta, *invZalpha, *invZbeta;
    Xalpha = &X.de_ele[nCol*alpha];
    Xbeta  = &X.de_ele[nCol*beta];
    invZalpha = &invZ.de_ele[nCol*alpha];
    invZbeta  = &invZ.de_ele[nCol*beta];
    
    for (int index2=0; index2<Ai.NonZeroCount; ++index2) {
#if DATA_CAPSULE
      int gamma = Ai.DataS[index2].vRow;
      int delta  = Ai.DataS[index2].vCol;
      double value2 = Ai.DataS[index2].vEle;
#else
      int gamma = Ai.row_index[index2];
      int delta  = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
#endif
      double plu = value2*invZbeta[delta] * Xalpha[gamma];

      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2*invZbeta[gamma] * Xalpha[delta];
        sum += plu2;
      }
    }
    ret += value1*sum;
    if (alpha==beta) {
      continue;
    }
    sum = 0.0;
    for (int index2=0; index2<Ai.NonZeroCount; ++index2) {
#if DATA_CAPSULE
      int gamma = Ai.DataS[index2].vRow;
      int delta  = Ai.DataS[index2].vCol;
      double value2 = Ai.DataS[index2].vEle;
#else
      int gamma = Ai.row_index[index2];
      int delta  = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
#endif      
      double plu = value2*invZalpha[delta] * Xbeta[gamma];
      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2*invZalpha[gamma] * Xbeta[delta];
        sum += plu2;
      }
    }
    ret += value1*sum;
  } // end of 'for (index1)'
  return;
}
#endif

inline void Newton::calF3_thread_1x1(double& ret,
				     DenseMatrix& X, DenseMatrix& invZ,
				     SparseMatrix& Ai, SparseMatrix& Aj)
{
  // Ai and Aj are SPARSE
  ret = 0.0;
  double sum = 0.0;
  int index1 = 0, index2 = 0;

  // rMessage("Aj.NonZeroCount = " << Aj.NonZeroCount);
#if DATA_CAPSULE
  int alpha = Aj.DataS[index1].vRow;
  int beta  = Aj.DataS[index1].vCol;
  double value1 = Aj.DataS[index1].vEle;
  int gamma = Ai.DataS[index2].vRow;
  int delta  = Ai.DataS[index2].vCol;
  double value2 = Ai.DataS[index2].vEle;
#else
  int alpha = Aj.row_index[index1];
  int beta  = Aj.column_index[index1];
  double value1 = Aj.sp_ele[index1];
  int gamma = Ai.row_index[index2];
  int delta  = Ai.column_index[index2];
  double value2 = Ai.sp_ele[index2];
#endif

  const int nCol = X.nCol;
  double *Xalpha, *Xbeta, *invZalpha, *invZbeta;
  Xalpha = &X.de_ele[nCol*alpha];
  Xbeta  = &X.de_ele[nCol*beta];
  invZalpha = &invZ.de_ele[nCol*alpha];
  invZbeta  = &invZ.de_ele[nCol*beta];
  
  double plu = value2*invZbeta[delta] * Xalpha[gamma];
  sum += plu;
  if (gamma!=delta) {
    double plu2 = value2*invZbeta[gamma] * Xalpha[delta];
    sum += plu2;
  }
  ret += value1*sum;
  if (alpha==beta) {
    return;
  }
  sum = 0.0;
  plu = value2*invZalpha[delta] * Xbeta[gamma];
  sum += plu;
  if (gamma!=delta) {
    double plu2 = value2*invZalpha[gamma] * Xbeta[delta];
    sum += plu2;
  }
  ret += value1*sum;

  return;
}

#if 0
inline void Newton::calF3_thread(double& ret,
				 DenseMatrix& X, DenseMatrix& invZ,
				 SparseMatrix& Ai, SparseMatrix& Aj)
{
  // Ai and Aj are SPARSE
  ret = 0.0;
  int index1, index2;
  int alpha, beta, gamma, delta;
  double value1, value2, plu, plu2;
  double sum;
  const int nCol = X.nCol;
  double *Xalpha, *Xbeta, *invZalpha, *invZbeta;
  for (index1=0; index1<Aj.NonZeroCount; ++index1) {
    alpha = Aj.DataS[index1].vRow;
    beta  = Aj.DataS[index1].vCol;
    value1 = Aj.DataS[index1].vEle;
    sum = 0.0;
    
    Xalpha = &X.de_ele[nCol*alpha];
    Xbeta  = &X.de_ele[nCol*beta];
    invZalpha = &invZ.de_ele[nCol*alpha];
    invZbeta  = &invZ.de_ele[nCol*beta];
    
    for (index2=0; index2<Ai.NonZeroCount; ++index2) {
      #if DATA_CAPSULE
      gamma = Ai.DataS[index2].vRow;
      delta  = Ai.DataS[index2].vCol;
      value2 = Ai.DataS[index2].vEle;
      #else
      gamma = Ai.row_index[index2];
      delta  = Ai.column_index[index2];
      value2 = Ai.sp_ele[index2];
      #endif
      plu = value2*invZbeta[delta] * Xalpha[gamma];
      sum += plu;
      if (gamma!=delta) {
        plu2 = value2*invZbeta[gamma] * Xalpha[delta];
        sum += plu2;
      }

      if (alpha==beta) {
	continue;
      }

      plu = value2*invZalpha[delta] * Xbeta[gamma];
      sum += plu;
      if (gamma!=delta) {
        plu2 = value2*invZalpha[gamma] * Xbeta[delta];
        sum += plu2;
      }
    }
    ret += value1*sum;
  } // end of 'for (index1)'
  return;
}
#endif

#if 1
inline void Newton::calF3_thread(double& ret,
				 DenseMatrix& X, DenseMatrix& invZ,
				 SparseMatrix& Ai, SparseMatrix& Aj)
{
  // Ai and Aj are SPARSE
  ret = 0.0;
  double sum;
  double *Xalpha, *Xbeta, *invZalpha, *invZbeta;
  const int nCol = X.nCol;
  //  rMessage("Aj.NonZeroCount = " << Aj.NonZeroCount);
  for (int index1=0; index1<Aj.NonZeroCount; ++index1) {
    #if DATA_CAPSULE
    int alpha = Aj.DataS[index1].vRow; 
    int beta  = Aj.DataS[index1].vCol;
    double value1 = Aj.DataS[index1].vEle;
    #else
    int alpha = Aj.row_index[index1]; 
    int beta  = Aj.column_index[index1];
    double value1 = Aj.sp_ele[index1];
    #endif
    Xalpha =    &X.de_ele[nCol*alpha];
    Xbeta  =    &X.de_ele[nCol*beta];
    invZalpha = &invZ.de_ele[nCol*alpha];
    invZbeta  = &invZ.de_ele[nCol*beta];
    sum = 0.0;
    //    rMessage("Ai.NonZeroCount = " << Ai.NonZeroCount);
    for (int index2=0; index2<Ai.NonZeroCount; ++index2) {
      #if DATA_CAPSULE
      int gamma = Ai.DataS[index2].vRow;
      int delta  = Ai.DataS[index2].vCol;
      double value2 = Ai.DataS[index2].vEle;
      #else
      int gamma = Ai.row_index[index2]; 
      int delta = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
      #endif
      double plu = value2 * invZbeta[delta] * Xalpha[gamma];
      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2 * invZbeta[gamma] * Xalpha[delta];
        sum += plu2;
      }
    }
    ret += value1*sum;
    if (alpha==beta) {
      continue;
    }
    sum = 0.0;
    for (int index2=0; index2<Ai.NonZeroCount; ++index2) {
      #if DATA_CAPSULE
      int gamma = Ai.DataS[index2].vRow;
      int delta  = Ai.DataS[index2].vCol;
      double value2 = Ai.DataS[index2].vEle;
      #else
      int gamma = Ai.row_index[index2]; 
      int delta  = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
      #endif
      double plu = value2 * invZalpha[delta] * Xbeta[gamma];
      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2 * invZalpha[gamma] * Xbeta[delta];
        sum += plu2;
      }
    }
    ret += value1*sum;
  } // end of 'for (index1)'
  return;
}
#endif

inline void Newton::calF3_thread_2(double& ret,
				   DenseMatrix& X, DenseMatrix& invZ,
				   SparseMatrix& Ai, SparseMatrix& Aj)
{
  // Ai and Aj are SPARSE
  ret = 0.0;
  double sum = 0.0;
  int index1, index2, alpha, beta, gamma, delta;
  double value1, plu;

  double *Xalpha, *Xbeta, *invZalpha, *invZbeta;
  const int nCol = X.nCol;

  // rMessage("Aj.NonZeroCount = " << Aj.NonZeroCount);
  for (index1=0; index1<Aj.NonZeroCount; ++index1) {
    #if DATA_CAPSULE
    alpha  = Aj.DataS[index1].vRow;
    beta   = Aj.DataS[index1].vCol;
    value1 = Aj.DataS[index1].vEle;
    #else
    alpha = Aj.row_index[index1]; 
    beta  = Aj.column_index[index1];
    value1 = Aj.sp_ele[index1];
    #endif

    Xalpha    = &X.de_ele[nCol*alpha];
    Xbeta     = &X.de_ele[nCol*beta];
    invZalpha = &invZ.de_ele[nCol*alpha];
    invZbeta  = &invZ.de_ele[nCol*beta];

    if (alpha != beta) {
      sum = 0.0;
      for (index2=0; index2 < Ai.NonZeroCount; index2++) {
	#if DATA_CAPSULE
	gamma = Ai.DataS[index2].vRow;
	delta  = Ai.DataS[index2].vCol;
	#else
	gamma = Ai.row_index[index2];
	delta  = Ai.column_index[index2];
	#endif
	
	if (gamma != delta) {
	  plu = invZbeta [delta] * Xalpha[gamma]
	      + invZbeta [gamma] * Xalpha[delta]
	      + invZalpha[delta] * Xbeta [gamma]
	      + invZalpha[gamma] * Xbeta [delta];
	} else {
	  plu = invZbeta [delta] * Xalpha[gamma]
	      + invZalpha[delta] * Xbeta [gamma];
	}
	#if DATA_CAPSULE
	sum += Ai.DataS[index2].vEle * plu;
	#else
	sum += Ai.sp_ele[index2] * plu;
	#endif
      }
      ret += value1 * sum;
    } else {
      sum = 0.0;
      for (index2=0; index2 < Ai.NonZeroCount; index2++) {
	#if DATA_CAPSULE
	gamma = Ai.DataS[index2].vRow;
	delta  = Ai.DataS[index2].vCol;
	#else
	gamma = Ai.row_index[index2];
	delta  = Ai.column_index[index2];
	#endif
	
	if (gamma != delta) {
	  plu = invZbeta[delta] * Xalpha[gamma]
	      + invZbeta[gamma] * Xalpha[delta];
	} else {
	  plu = invZbeta[delta] * Xalpha[gamma];
	}
	#if DATA_CAPSULE
	sum += Ai.DataS[index2].vEle * plu;
	#else
	sum += Ai.sp_ele[index2] * plu;
	#endif	
      }
      ret += value1 * sum;
    }
  }
  return;
}

void Newton::compute_bMat_dense_SDP_thread(InputData& inputData,
                                    Solutions& currentPt,
                                    WorkVariables& work,
                                    ComputeTime& com)
{
  // rMessage("NUM_THREADS = " << NUM_THREADS);
  pthread_t*  handle;
  NewArray(handle,pthread_t,NUM_THREADS);
  thread_arg_t* targ;
  NewArray(targ,thread_arg_t,NUM_THREADS);
  #ifdef PRINT_LOAD_BALANCE
  double* ThreadTime;
  NewArray(ThreadTime,double,NUM_THREADS);
  for (int k=0; k<NUM_THREADS; k++) {
    ThreadTime[k] = 0.0;
  }
  #endif
  
#ifdef GOTO_BLAS
  if (Calc_F1 > 1) {
    goto_set_num_threads(NUM_GOTOBLAS);
  }
#endif

  int ret;
  ret = pthread_mutex_init(&job_mutex, NULL);
  if (ret  != 0) {
    rError("pthread_mutex_init error");
  }
  ret = pthread_cond_init(&job_cond, NULL);
  if (ret != 0) {
    rError("pthread_cond_init error");
  }
  int m = currentPt.mDim;
  int SDP_nBlock = inputData.SDP_nBlock;

  for (int k=0; k<NUM_THREADS; k++) {
    targ[k].mDIM        = m;
    targ[k].SDP_nBlock  = SDP_nBlock;
    // targ[k].bMat        = &bMat;
    targ[k].useFormula  = useFormula;
    targ[k].inputData   = &inputData;
    targ[k].currentPt   = &currentPt;
    targ[k].work        = &work;
    targ[k].com         = &com;
    targ[k].symmetric_b = symmetric_b;
    targ[k].pB          = pB;
    targ[k].maxnp       = maxnp;
  }

  for (int l=0; l<SDP_nBlock; l++) {
    Column_Number = 0;
    for (int k=0; k<NUM_THREADS; k++) {
      targ[k].Block_Number = l;
      targ[k].thread_num   = k;
      targ[k].pBtime      = 0.0;
      pthread_create(&handle[k], NULL,
		     compute_bMat_dense_SDP_thread_func,
		     (void *)&targ[k]);
    }
    
    for (int k=0; k<NUM_THREADS; k++) {
      pthread_join(handle[k], NULL);
#if PRINT_LOAD_BALANCE
      ThreadTime[k] += targ[k].pBtime;
#endif
    }
  }
  DeleteArray(handle);
  DeleteArray(targ);
#if PRINT_LOAD_BALANCE
  for (int k=0; k<NUM_THREADS; k++) {
    rMessage("Make bMat on [" << k
	     << "/" << MpiSt::iam << "] (SDP) = "
	     << ThreadTime[k]);
  }
  DeleteArray(ThreadTime);
#endif

  ret = pthread_mutex_destroy(&job_mutex);
  if (ret != 0) {
    rError("pthread_mutex_destroy error in sdpa_newton.cpp");
  }
  ret = pthread_cond_destroy(&job_cond);
  if (ret != 0) {
    rError("pthread_cond_destroy error in sdpa_newton.cpp");
  }

#ifdef GOTO_BLAS
  if (Calc_F1 > 1) {
    goto_set_num_threads(-1);
  }
#endif

}

void* Newton::compute_bMat_dense_SDP_thread_func(void *arg)
{
#if PRINT_LOAD_BALANCE
  TimeStart(START1);
#endif
  int l, m;
  int k1;
  int SDP_nBlock;
  DenseMatrix work1, work2;

  thread_arg_t *targ = (thread_arg_t *)arg;
 
  l = targ->Block_Number;
  m = targ->mDIM;
  SDP_nBlock = targ->SDP_nBlock;

  //  printf("targ-> Block_Number = %d\n", targ-> Block_Number); 

  // DenseMatrix& xMat = targ->currentPt->xMat.SDP_block[l];
  // DenseMatrix& invzMat = targ->currentPt->invzMat.SDP_block[l];
  work1.initialize(targ->work->DLS1.SDP_block[l].nRow,
		   targ->work->DLS1.SDP_block[l].nCol,
		   DenseMatrix::DENSE);
  work2.initialize(targ->work->DLS2.SDP_block[l].nRow,
		   targ->work->DLS2.SDP_block[l].nCol,
		   DenseMatrix::DENSE);

  while(1) {
    pthread_mutex_lock(&job_mutex);
    k1 = Column_Number++;
    pthread_mutex_unlock(&job_mutex);

    if (k1 >= targ->inputData->SDP_nConstraint[l]) {
      break;
    }

    int j = targ->inputData->SDP_constraint[l][k1];
    if (j % MpiSt::nprocs != MpiSt::iam) {
      continue;
    }
    int jb = targ->inputData->SDP_blockIndex[l][k1];
    int jnz = targ->inputData->A[j].SDP_sp_block[jb].NonZeroEffect;
    SparseMatrix& Aj = targ->inputData->A[j].SDP_sp_block[jb];
    
    FormulaType formula = targ->useFormula[j][jb];

    TimeStart(B_NDIAG_START1);
    TimeStart(B_NDIAG_START2);

    bool hasF2Gcal = false;
    #if 0
    printf("thread_num = %d, mutex_flag = %d\n",
	   targ->thread_num, mutex_flag);
    #endif

    if (formula==F1) {
      pthread_mutex_lock(&job_mutex);
      Lal::let(work1,'=',targ->currentPt->xMat.SDP_block[l],'*',Aj);
      Lal::let(work2,'=',work1,'*',targ->currentPt->invzMat.SDP_block[l]);
      pthread_mutex_unlock(&job_mutex);
    } else if (formula==F2) {
      pthread_mutex_lock(&job_mutex);
      // Lal::let(work1,'=',Ai,'*',targ->currentPt->invzMat.SDP_block[l]);
      Lal::let(work1,'=',targ->currentPt->xMat.SDP_block[l],'*',Aj);
      pthread_mutex_unlock(&job_mutex);
      hasF2Gcal = false;
    }
    TimeEnd(B_NDIAG_END2);
    targ->com->B_PRE += TimeCal(B_NDIAG_START2,B_NDIAG_END2);

    for (int k2=targ->inputData->SDP_nConstraint[l]-1; k2 >= 0; k2--) {
      int i = targ->inputData->SDP_constraint[l][k2];
      int ib = targ->inputData->SDP_blockIndex[l][k2];
      int inz = targ->inputData->A[i].SDP_sp_block[ib].NonZeroEffect;
      SparseMatrix& Ai = targ->inputData->A[i].SDP_sp_block[ib];

      if (targ->symmetric_b) {
	if ((jnz < inz) || ( (inz == jnz) && (i<j))) {
	  continue;
	}
      }
      else {
	if (i<j) {
	  continue;
	}
      }
      double value;
      switch (formula) {
      case F3:
	// rMessage("calF3");
#if DEBUG
	printf("F3 in %d\n", targ->thread_num);
#endif
	if ((Ai.NonZeroCount == 1) && (Aj.NonZeroCount == 1))
	  calF3_thread_1x1(value,
			   targ->currentPt->xMat.SDP_block[l],
			   targ->currentPt->invzMat.SDP_block[l],
			   Ai,Aj);
	else
	  calF3_thread_2(value,
			 targ->currentPt->xMat.SDP_block[l],
			 targ->currentPt->invzMat.SDP_block[l],
			 Ai,Aj);
	break;
      case F1:
	// rMessage("calF1");
#if DEBUG
	printf("F1 in %d\n", targ->thread_num);
#endif
	calF1_thread(value,work2,Ai);
	break;
      case F2:
	// rMessage("calF2 ");
#if DEBUG
	printf("F2 in %d\n", targ->thread_num);
#endif
	calF2_thread(value,work1,work2,
		     targ->currentPt->invzMat.SDP_block[l],Ai,hasF2Gcal);
	break;
      } // end of switch
      targ->pB[i + m*(j/MpiSt::nprocs)] += value;
    } // end of 'for (int j)'

    TimeEnd(B_NDIAG_END1);
    double t = TimeCal(B_NDIAG_START1,B_NDIAG_END1);
    switch (formula) {
    case F1: targ->com->B_F1 += t; break;
    case F2: targ->com->B_F2 += t; break;
    case F3: targ->com->B_F3 += t; break;
    }
  }

  work1.terminate();
  work2.terminate();

#if PRINT_LOAD_BALANCE
  TimeEnd(END1);
  targ->pBtime = TimeCal(START1,END1);
#endif
  return NULL;
}

void Newton::compute_bMat_dense_SDP(InputData& inputData,
                                    Solutions& currentPt,
                                    WorkVariables& work,
                                    ComputeTime& com)
{
  TimeStart(B_NDIAG_START1);
  // int m = currentPt.mDim;
  int SDP_nBlock = inputData.SDP_nBlock;

  for (int l=0; l<SDP_nBlock; ++l) {
    DenseMatrix& xMat = currentPt.xMat.SDP_block[l];
    DenseMatrix& invzMat = currentPt.invzMat.SDP_block[l];
    DenseMatrix& work1 = work.DLS1.SDP_block[l];
    DenseMatrix& work2 = work.DLS2.SDP_block[l];

    for (int k1=0; k1<inputData.SDP_nConstraint[l]; k1++) {
      int j = inputData.SDP_constraint[l][k1];
      if (j % MpiSt::nprocs != MpiSt::iam) {
	continue;
      }
      int jb = inputData.SDP_blockIndex[l][k1];
      int jnz = inputData.A[j].SDP_sp_block[jb].NonZeroEffect;
      SparseMatrix& Aj = inputData.A[j].SDP_sp_block[jb];

      FormulaType formula = useFormula[j][jb];
      // ---------------------------------------------------
      // formula = F3; // this is force change
      // ---------------------------------------------------
      TimeStart(B_NDIAG_START2);

      bool hasF2Gcal = false;
      if (formula==F1) {
	// Lal::let(work1,'=',Ai,'*',invzMat);
	// Lal::let(work2,'=',xMat,'*',work1);
	Lal::let(work1,'=',xMat,'*',Aj);
	Lal::let(work2,'=',work1,'*',invzMat);
      } else if (formula==F2) {
	// Lal::let(work1,'=',Ai,'*',invzMat);
	Lal::let(work1,'=',xMat,'*',Aj);
	hasF2Gcal = false;
	// Lal::let(gMat.ele[l],'=',xMat.ele[l],'*',fMat.ele[l]);
      }
      TimeEnd(B_NDIAG_END2);
      com.B_PRE += TimeCal(B_NDIAG_START2,B_NDIAG_END2);

      for (int k2=0; k2<inputData.SDP_nConstraint[l]; k2++) {
	int i = inputData.SDP_constraint[l][k2];
	int ib = inputData.SDP_blockIndex[l][k2];
	int inz = inputData.A[i].SDP_sp_block[ib].NonZeroEffect;
	SparseMatrix& Ai = inputData.A[i].SDP_sp_block[ib];
          
	// Select the formula A[i] or the formula A[j].
	// Use formula that has more NonZeroEffects than others.
	// We must calculate i==j.
	if (symmetric_b) {
	  if ((jnz < inz) || ( (inz == jnz) && (i<j))) {
	    continue;
	  }
	}
	else {
	  if (i<j) {
	    continue;
	  }
	}

	double value;
	switch (formula) {
	case F3:
	  // rMessage("calF3");
	  calF3(value,xMat,invzMat,Ai,Aj);
	  break;
	case F1:
	  // rMessage("calF1");
	  calF1(value,work2,Ai);
	  break;
	case F2:
	  // rMessage("calF2 ");
	  calF2(value,work1,work2,invzMat,Ai,hasF2Gcal);
	  // calF1(value2,gMat.ele[l],A[j].ele[l]);
	  // rMessage("calF2:  " << (value-value2));
	  break;
	} // end of switch

	#if 0
	rMessage("i = "<< i << " j = " << j 
		 << " v = " << value
		 << " loc = " << ( i/MpiSt::nprocs + maxmp*j));
	#endif

	pB[i + m*(j/MpiSt::nprocs)] += value;
      } // end of 'for (int j)'

    } // end of 'for (int i)'
  } // end of 'for (int l)'
  TimeEnd(B_NDIAG_END1);
  // double t = TimeCal(B_NDIAG_START1,B_NDIAG_END1);
  // com.makebMat += t; // This addition will be done in Make_bMat
}

void Newton::computeStartEndIndices(int& start, int& end,
				    int iam, int nprocs, int length)
{
  int shou  = length / nprocs;
  int amari = length % nprocs;
  start = iam * shou;
  if (iam < amari) {
    start += iam;
  }
  else {
    start += amari;
  }
  end = start + shou;
  if (iam < amari) {
    end++;
  }
}

void Newton::setNumThreads(FILE* Display, FILE* fpOut, int NumThreads)
{
  if (NumThreads == 0) { // Automatic from OMP_NUM_THREADS
    char* env1      = NULL;
    env1 = getenv("OMP_NUM_THREADS");
    if (env1 != NULL) {
      NUM_THREADS = atoi(env1);
    }
    else {
      NUM_THREADS = 1;
    }
  }
  else {
    NUM_THREADS = NumThreads;
  }
  if (Display) {
    fprintf(Display,"NumNodes    is set as %d\n", MpiSt::nprocs);
    fprintf(Display,"NumThreads  is set as %d\n", NUM_THREADS);
  }
  if (fpOut) {
    fprintf(fpOut,  "NumNodes    is set as %d\n", MpiSt::nprocs);
    fprintf(fpOut,  "NumThreads  is set as %d\n", NUM_THREADS);
  }
  
  // rMessage("MyNumThreads = " << NUM_THREADS);
  
  #ifdef GOTO_BLAS
  NUM_GOTOBLAS = NUM_THREADS / 2;
  if (NUM_GOTOBLAS <= 0) {
    NUM_GOTOBLAS = 1;
  }
  goto_set_num_threads(NUM_GOTOBLAS);
  if (Display) {
    fprintf(Display,"GotoThreads is set as %d\n", NUM_GOTOBLAS);
  }
  if (fpOut) {
    fprintf(fpOut,  "GotoThreads is set as %d\n", NUM_GOTOBLAS);
  }
  #endif

}

void Newton::compute_bMat_sparse_thread_invoke(InputData& inputData,
					       Solutions& currentPt,
					       WorkVariables& work,
					       ComputeTime& com)
{
  pthread_t*  handle;
  NewArray(handle,pthread_t,NUM_THREADS);
  thread_arg_t* targ;
  NewArray(targ,thread_arg_t,NUM_THREADS);

#ifdef GOTO_BLAS
  if (Calc_F1 > 1) {
    goto_set_num_threads(NUM_GOTOBLAS);
  }
#endif
  int m = currentPt.mDim;
  int SDP_nBlock = inputData.SDP_nBlock;

  for (int k=0; k<NUM_THREADS; k++) {
    targ[k].thread_num      = k;
    targ[k].mDIM            = m;

    targ[k].SDP_nBlock      = SDP_nBlock;
    targ[k].SDP_number      = SDP_number;
    targ[k].SDP_constraint1 = SDP_constraint1;
    targ[k].SDP_constraint2 = SDP_constraint2;
    targ[k].SDP_blockIndex1 = SDP_blockIndex1;
    targ[k].SDP_blockIndex2 = SDP_blockIndex2;
    targ[k].SDP_location_sparse_bMat = SDP_location_sparse_bMat;    

    targ[k].LP_nBlock       = LP_nBlock;
    targ[k].LP_number       = LP_number;
    targ[k].LP_constraint1  = LP_constraint1;
    targ[k].LP_constraint2  = LP_constraint2;
    targ[k].LP_blockIndex1  = LP_blockIndex1;
    targ[k].LP_blockIndex2  = LP_blockIndex2;
    targ[k].LP_location_sparse_bMat = LP_location_sparse_bMat;

    targ[k].sparse_bMat     = &sparse_bMat;
    targ[k].useFormula      = useFormula;
    targ[k].inputData       = &inputData;
    targ[k].currentPt       = &currentPt;
    targ[k].work            = &work;
    targ[k].com             = &com;

    targ[k].threadSchurStart = threadSchurStart[k];
    targ[k].threadSchurEnd   = threadSchurEnd[k];
    targ[k].threadSDPLocStart = threadSDPLocStart[k];
    targ[k].threadSDPLocEnd   = threadSDPLocEnd[k];
    targ[k].threadLPLocStart  = threadLPLocStart[k];
    targ[k].threadLPLocEnd    = threadLPLocEnd[k];
   }

  for (int k=0; k<NUM_THREADS; k++) {
    // targ[k].Block_Number = l;
    if (threadSchurStart[k] < threadSchurEnd[k]) {
      // rMessage("invoke No" << k << " thread");
      pthread_create(&handle[k], NULL,
		     compute_bMat_sparse_thread_func,
		     (void *)&targ[k]);
    }
  }
    
  for (int k=0; k<NUM_THREADS; k++) {
    if (threadSchurStart[k] < threadSchurEnd[k]) {
      pthread_join(handle[k], NULL);
    }
  }
  DeleteArray(handle);
  DeleteArray(targ);
  #ifdef GOTO_BLAS
  goto_set_num_threads(-1);
  #endif
}

void* Newton::compute_bMat_sparse_thread_func(void* arg)

{
#if PRINT_LOAD_BALANCE
  TimeStart(START1);
#endif
  compute_bMat_sparse_SDP_thread(arg);
#if PRINT_LOAD_BALANCE
  TimeEnd(END1);
#endif
  //   compute_bMat_sparse_SOCP(inputData,currentPt,work,com);
#if PRINT_LOAD_BALANCE
  TimeStart(START3);
#endif
  compute_bMat_sparse_LP_thread(arg);
#if PRINT_LOAD_BALANCE
  TimeEnd(END3);
  rMessage("Make bMat on [" << (((thread_arg_t *)arg)->thread_num)
	   << "/" << MpiSt::iam << "] (SDP) = " << TimeCal(START1,END1)
	   << " (LP) = " << TimeCal(START3,END3));
#endif
  return NULL;
}

void Newton::compute_bMat_sparse_SDP_thread(void *arg)
{
  #if PRINT_LOAD_BALANCE
  TimeStart(START1);
  #endif
  thread_arg_t *targ = (thread_arg_t *)arg;
  const int SDP_nBlock = targ->SDP_nBlock;
  DenseMatrix work1, work2;
  double* bMat_sp_ele = targ->sparse_bMat->sp_ele;
  for (int l=0; l<SDP_nBlock; ++l) {
    //  printf("targ-> Block_Number = %d\n", targ-> Block_Number);
    int previous_j = -1;

    // DenseMatrix& xMat = targ->currentPt->xMat.SDP_block[l];
    // DenseMatrix& invzMat = targ->currentPt->invzMat.SDP_block[l];
    work1.initialize(targ->work->DLS1.SDP_block[l].nRow,
		     targ->work->DLS1.SDP_block[l].nCol,
		     DenseMatrix::DENSE);
    work2.initialize(targ->work->DLS2.SDP_block[l].nRow,
		     targ->work->DLS2.SDP_block[l].nCol,
		     DenseMatrix::DENSE);

    // TimeStart(B_NDIAG_START1);
    const int startLoc = targ->threadSDPLocStart[l];
    const int endLoc   = targ->threadSDPLocEnd[l];
    #if 0
    rMessage("Loc [" << startLoc << "," <<  endLoc << "]"
	     << " block("<<l<<") by ("
	     << targ->thread_num << "/" << MpiSt::iam << ")");
    #endif
    for (int iter = startLoc; iter < endLoc; ++iter) {
      const int j = targ->SDP_constraint2[l][iter];
      const int jb = targ->SDP_blockIndex2[l][iter];
      SparseMatrix& Aj = targ->inputData->A[j].SDP_sp_block[jb];
    
      const FormulaType formula = targ->useFormula[j][jb];
      // rMessage("FormulaType = " << formula);

      if (j != previous_j){
	// TimeStart(B_NDIAG_START2);

	if (formula==F1) {
	  pthread_mutex_lock(&job_mutex);
	  Lal::let(work1,'=',targ->currentPt->xMat.SDP_block[l],'*',Aj);
	  Lal::let(work2,'=',work1,'*',targ->currentPt->invzMat.SDP_block[l]);
	  pthread_mutex_unlock(&job_mutex);
	} else if (formula==F2) {
	  // Lal::let(work1,'=',Ai,'*',targ->currentPt->invzMat.SDP_block[l]);
	  Lal::let(work1,'=',targ->currentPt->xMat.SDP_block[l],'*',Aj);
	}
	// TimeEnd(B_NDIAG_END2);
	// targ->com->B_PRE += TimeCal(B_NDIAG_START2,B_NDIAG_END2);
      }
    
      const int i = targ->SDP_constraint1[l][iter];
      const int ib = targ->SDP_blockIndex1[l][iter];
      SparseMatrix& Ai = targ->inputData->A[i].SDP_sp_block[ib];

      double value;
      bool dummyHasF2Gcal = true;
      switch (formula) {
      case F3:
	//      rMessage("calF3");
	if ((Ai.NonZeroCount == 1) && (Aj.NonZeroCount == 1))
	  calF3_thread_1x1(value,
			   targ->currentPt->xMat.SDP_block[l],
			   targ->currentPt->invzMat.SDP_block[l],
			   Ai,Aj);
	else
	  calF3_thread_2(value,
			 targ->currentPt->xMat.SDP_block[l],
			 targ->currentPt->invzMat.SDP_block[l],
			 Ai,Aj);
	break;
      case F1:
	//      rMessage("calF1");
	calF1_thread(value,work2,Ai);
	break;
      case F2:
	//      rMessage("calF2 ");
	calF2_thread(value,work1,work2,
		     targ->currentPt->invzMat.SDP_block[l],Ai,dummyHasF2Gcal);
	break;
      } // end of switch
      bMat_sp_ele[targ->SDP_location_sparse_bMat[l][iter]] += value;
      previous_j = j;
    }
#if 0
    // TimeEnd(B_NDIAG_END1);
    const double t = TimeCal(B_NDIAG_START1,B_NDIAG_END1);
    switch (formula) {
    case F1: targ->com->B_F1 += t; break;
    case F2: targ->com->B_F2 += t; break;
    case F3: targ->com->B_F3 += t; break;
    }
#endif

    work1.terminate();
    work2.terminate();
  }

  #if 0
  TimeEnd(END1);
  rMessage("SDP cone = " << TimeCal(START1,END1)
	   << " by (" << targ->thread_num << "/" << MpiSt::iam << ")");
  #endif
}

void Newton::compute_bMat_sparse_SDP(InputData& inputData,
				     Solutions& currentPt,
				     WorkVariables& work,
				     ComputeTime& com)
{
  TimeStart(B_NDIAG_START1);
  // TimeStart(B_NDIAG_START2);
  for (int l=0; l<SDP_nBlock; ++l) {
    DenseMatrix& xMat = currentPt.xMat.SDP_block[l];
    DenseMatrix& invzMat = currentPt.invzMat.SDP_block[l];
    DenseMatrix& work1 = work.DLS1.SDP_block[l];
    DenseMatrix& work2 = work.DLS2.SDP_block[l];
    int previous_j = -1;

    for (int iter = 0; iter < SDP_number[l]; iter++){
      //      TimeStart(B_NDIAG_START1);
      int j = SDP_constraint2[l][iter];
      int jb = SDP_blockIndex2[l][iter];
      SparseMatrix& Aj = inputData.A[j].SDP_sp_block[jb];
      FormulaType formula = useFormula[j][jb];
      
      if (j != previous_j){
	// ---------------------------------------------------
	// formula = F3; // this is force change
	// ---------------------------------------------------
	
	if (formula==F1) {
	  // Lal::let(work1,'=',Ai,'*',invzMat);
	  // Lal::let(work2,'=',xMat,'*',work1);
	  Lal::let(work1,'=',xMat,'*',Aj);
	  Lal::let(work2,'=',work1,'*',invzMat);
	} else if (formula==F2) {
	  // Lal::let(work1,'=',Ai,'*',invzMat);
	  Lal::let(work1,'=',xMat,'*',Aj);
	  // Lal::let(gMat.ele[l],'=',xMat.ele[l],'*',fMat.ele[l]);
	}
	// TimeEnd(B_NDIAG_END2);
	// com.B_PRE += TimeCal(B_NDIAG_START2,B_NDIAG_END2);
      }
      
      int i = SDP_constraint1[l][iter];
      int ib = SDP_blockIndex1[l][iter];
      SparseMatrix& Ai = inputData.A[i].SDP_sp_block[ib];
      // rMessage("B(" << i <<","<< j <<") is computed");

      #if 0
      printf("iter = %d, loc = %d, i = %d, ib = %d, j = %d, jb = %d on CPU %d\n",
	     iter, SDP_location_sparse_bMat[l][iter], i, ib, j, jb,
	     MpiSt::iam);
      #endif
      
      double value;
      bool dummyHasF2Gcal = true;
      switch (formula) {
      case F3:
	// rMessage("calF3");
	if ((Ai.NonZeroCount == 1) && (Aj.NonZeroCount == 1))
	  calF3_thread_1x1(value,xMat,invzMat,Ai,Aj);
	else
	  calF3_thread(value,xMat,invzMat,Ai,Aj);
        // calF3(value,work1,work2,xMat,invzMat,Ai,Aj);
	break;
      case F1:
	// rMessage("calF1");
	calF1(value,work2,Ai);
	break;
      case F2:
	// rMessage("calF2 ");
	calF2(value,work1,work2,invzMat,Ai,dummyHasF2Gcal);
	// calF1(value2,gMat.ele[l],A[j].ele[l]);
	// rMessage("calF2:  " << (value-value2));
	break;
      } // end of switch
      sparse_bMat.sp_ele[SDP_location_sparse_bMat[l][iter]] += value;
      previous_j = j;
    } // end of 'for (int index)'
  } // end of 'for (int l)'
#if 0
  TimeEnd(B_NDIAG_END1);
  double t = TimeCal(B_NDIAG_START1,B_NDIAG_END1);
  switch (formula) {
  case F1: com.B_F1 += t; break;
  case F2: com.B_F2 += t; break;
  case F3: com.B_F3 += t; break;
  }
#else
  TimeEnd(B_NDIAG_END1);
  // double t = TimeCal(B_NDIAG_START1,B_NDIAG_END1);
  // com.makebMat += t; // This addition will be done in Make_bMat
#endif
}

#if 0
void Newton::compute_bMat_dense_SCOP(InputData& inputData,
				     Solutions& currentPt,
				     WorkVariables& work,
				     ComputeTime& com)
{
    rError("current version does not support SOCP");
}

void Newton::compute_bMat_sparse_SOCP(InputData& inputData,
				      Solutions& currentPt,
				      WorkVariables& work,
				      ComputeTime& com)
{
    rError("current version does not support SOCP");
}
#endif

void Newton::compute_bMat_dense_LP(InputData& inputData,
				   Solutions& currentPt,
				   WorkVariables& work,
				   ComputeTime& com)
{
  // int m = currentPt.mDim;
  int LP_nBlock = inputData.LP_nBlock;

  TimeEnd(B_DIAG_START1);
  for (int l=0; l<LP_nBlock; ++l) {
    double xMat = currentPt.xMat.LP_block[l];
    double invzMat = currentPt.invzMat.LP_block[l];

      for (int k1=0; k1<inputData.LP_nConstraint[l]; k1++) {
	int j = inputData.LP_constraint[l][k1];
	if (j % MpiSt::nprocs != MpiSt::iam) {
	  continue;
	}
	int jb = inputData.LP_blockIndex[l][k1];
	//	int inz = inputData.A[i].LP_sp_block[ib].NonZeroEffect;
	double Aj = inputData.A[j].LP_sp_block[jb];

	// NOTE: for parallel  k2 start from 0
	// NOTE: for serial    k2 start from k1
	for (int k2=0; k2<inputData.LP_nConstraint[l]; k2++) {
	  int i = inputData.LP_constraint[l][k2];
	  int ib = inputData.LP_blockIndex[l][k2];
	  if (i < j) {
	    continue; // symmetric_b
	  }
	  //	  int jnz = inputData.A[j].LP_sp_block[jb].NonZeroEffect;
	  double Ai = inputData.A[i].LP_sp_block[ib];

	  double value;
	  value = xMat * invzMat * Ai * Aj;
	  pB[i+m*(j/MpiSt::nprocs)] += value;

	} // end of 'for (int j)'
      } // end of 'for (int i)'
  } // end of 'for (int l)'
  TimeEnd(B_DIAG_END1);
  com.B_DIAG += TimeCal(B_DIAG_START1,B_DIAG_END1);

}

void Newton::compute_bMat_sparse_LP_thread(void* arg)
{
  thread_arg_t *targ = (thread_arg_t *)arg;
  const int LP_nBlock = targ->LP_nBlock;
  double* bMat_sp_ele = targ->sparse_bMat->sp_ele;
  int** LP_location   = targ->LP_location_sparse_bMat;
  
  // TimeEnd(B_DIAG_START1);
  for (int l=0; l<LP_nBlock; ++l) {
    
    const double xMat = targ->currentPt->xMat.LP_block[l];
    const double invzMat = targ->currentPt->invzMat.LP_block[l];
    const int startLoc = targ->threadLPLocStart[l];
    const int endLoc   = targ->threadLPLocEnd[l];
    
    for (int iter = startLoc; iter < endLoc; iter++){
      const int j = targ->LP_constraint2[l][iter];
      const int jb = targ->LP_blockIndex2[l][iter];
      const double Aj = targ->inputData->A[j].LP_sp_block[jb];

      const int i = targ->LP_constraint1[l][iter];
      const int ib = targ->LP_blockIndex1[l][iter];
      const double Ai = targ->inputData->A[i].LP_sp_block[ib];
      
      const double value = xMat * invzMat * Ai * Aj;
      bMat_sp_ele[LP_location[l][iter]] += value;
    } // end of 'for (int iter)
  } // end of 'for (int l)'
  // TimeEnd(B_DIAG_END1);
  // targ->com->B_DIAG += TimeCal(B_DIAG_START1,B_DIAG_END1);
}



void Newton::Make_bMat(InputData& inputData,
		       Solutions& currentPt,
		       WorkVariables& work,
		       ComputeTime& com)
{
  MpiSt::barrier();
  TimeStart(START3);
  if (bMat_type == SPARSE){
    // set sparse_bMat zero
    // sdpa_dset(sparse_bMat.NonZeroCount, 0.0, sparse_bMat.sp_ele, 1);
    if (mySchurLength != 0) {
      double* targetAddress = &sparse_bMat.sp_ele[mySchurStart];
      sdpa_dset(mySchurLength, 0.0, targetAddress,1);
      compute_bMat_sparse_thread_invoke(inputData, currentPt,
					work, com);     
    }
   } else {
    // pB = 0
    sdpa_dset(m*maxnp, 0.0, pB, 1);
    compute_bMat_dense_SDP_thread(inputData,currentPt,work,com);
    //    compute_bMat_dense_SOCP(inputData,currentPt,work,com);
    compute_bMat_dense_LP(inputData,currentPt,work,com);
  }
  TimeEnd(END3);
  MpiSt::barrier();
  #if PRINT_LOAD_BALANCE
  rMessage("Make bMat iteration = " << TimeCal(START3,END3)
	   << " seconds" );
  #endif
  TimeEnd(END4);
  // The elapsed time to add to makebMat is measured after barrier();
  com.makebMat += TimeCal(START3,END4);
  
  #if 0
  rMessage("bMat =  ");
  if (bMat_type == DENSE) {
    display_pB();
  }
  else {
    display_sparse_bMat();
  }
  #endif

  if (bMat_type == DENSE) {
    TimeStart(START_COPYBMAT);
    MpiSt::barrier();
    Cpdgemr2d(m,m,pB,1,1,descB,pB2,1,1,descB2,MpiSt::ictxt);
    MpiSt::barrier();
    TimeStart(END_COPYBMAT);
    com.copybMat += TimeCal(START_COPYBMAT,END_COPYBMAT);
    
    TimeStart(START_SYMMBMAT);
    if (symmetric_b) {
      int mm1 = m-1;
      int two = 2;
      pdtradd_((char*)"Lower",(char*)"Transpose",&mm1,&mm1,&DONE,
	       pB2,&IONE,&two,descB2,&DONE,pB2,&two,&IONE,descB2);
    }
    // display_pB2();
    MpiSt::barrier();
    TimeStart(END_SYMMBMAT);
    com.symmetrisebMat += TimeCal(START_SYMMBMAT,END_SYMMBMAT);
    // display_pB();
    // display_pB2();
  }
  else {
    #if 0
    // Nothing is Needed for SPARSE
    TimeStart(START_COPYBMAT);
    TimeStart(END_COPYBMAT);
    com.copybMat += TimeCal(START_COPYBMAT,END_COPYBMAT);
    #endif
  }
}

bool Newton::compute_DyVec(Newton::WHICH_DIRECTION direction,
			   int m,
			   InputData& inputData,
			   Chordal& chordal,
			   Solutions& currentPt,
			   WorkVariables& work,
			   ComputeTime& com,
			   FILE* Display, FILE* fpOut)
{
  if (direction == PREDICTOR) {
    TimeStart(START3_2);
    
    if (bMat_type == SPARSE){
      bool ret = chordal.factorizeSchur(m,diagonalIndex, Display, fpOut);
      if (ret == SDPA_FAILURE) {
	return SDPA_FAILURE;
      }
    } else {

      #if 0
      bool ret = Lal::choleskyFactorWithAdjust(bMat);
      if (ret == SDPA_FAILURE) {
        return SDPA_FAILURE;
      }
      #endif

      int info = 0;
      rpdpotrf_((char*)"Lower",&m,pB2,&IONE,&IONE,descB2,&info);
      if (info != 0) {
	if (MpiSt::iam == 0) {
	  rMessage("Cholesky miss :: info = " << info);
	}
	return SDPA_FAILURE;
      }
    }
    // rMessage("Cholesky of bMat =  ");
    // bMat.display();
    // sparse_bMat.display();
    TimeEnd(END3_2);
    com.choleskybMat += TimeCal(START3_2,END3_2);
  }
  // bMat is already cholesky factorized.


  TimeStart(START4);
  if (bMat_type == SPARSE){
    // already copied to DyVec;
    // DyVec.copyFrom(gVec);
    chordal.solveSchur(DyVec);
  } else {
    #if 0
    Lal::let(DyVec,'=',bMat,'/',gVec);
    #endif

    // parallel start
    pdtrsv_((char*)"Lower",(char*)"NoTrans",(char*)"NonDiag",
	    &m,pB2,&IONE,&IONE,descB2,
	    pg2,&IONE,&IONE,descg2,&IONE);
    pdtrsv_((char*)"Lower",(char*)"Trans",(char*)"NonDiag",
	    &m,pB2,&IONE,&IONE,descB2,
	    pg2,&IONE,&IONE,descg2,&IONE); 
    // parallel end
  }
  TimeEnd(END4);
  com.solve += TimeCal(START4,END4);

  // display_pg2();
  // parallel start
  TimeStart(START_COPYDYVEC);
  if (bMat_type == DENSE){
    Cpdgemr2d(m,1,pg2,1,1,descg2,DyVec.ele,1,1,descDy,MpiSt::ictxt);
  }
  MpiCopy::allSendRecieveD(m,DyVec.ele);
  TimeEnd(END_COPYDYVEC);
  com.copyDyVec += TimeCal(START_COPYDYVEC,END_COPYDYVEC);
  // parallel end
  // rMessage("DyVec =  ");
  // DyVec.display();

  return SDPA_SUCCESS;
}

void Newton::compute_DzMat(InputData& inputData,
			   Residuals& currentRes,
			   Phase& phase,
			   ComputeTime& com)
{
  TimeStart(START_SUMDZ);
  inputData.multi_plusToA(DyVec, DzMat);
  Lal::let(DzMat,'=',DzMat,'*',&DMONE);
  if (phase.value == SolveInfo:: pFEAS
      || phase.value == SolveInfo::noINFO) {
    Lal::let(DzMat,'=',DzMat,'+',currentRes.dualMat);
  }
  TimeEnd(END_SUMDZ);
  com.sumDz += TimeCal(START_SUMDZ,END_SUMDZ);
  // rMessage("DzMat = ");
  // DzMat.display();
}

void Newton::compute_DxMat(Solutions& currentPt,
			   WorkVariables& work,
			   ComputeTime& com)
{
  TimeStart(START_DX);
  // work.DLS1 = dX dZ Z^{-1}
  Jal::ns_jordan_triple_product(work.DLS1,currentPt.xMat,DzMat,
				currentPt.invzMat,work.DLS2);
  // dX = R Z^{-1} - dX dZ Z^{-1}
  Lal::let(DxMat,'=',r_zinvMat,'+',work.DLS1,&DMONE);
  TimeEnd(END_DX);
  TimeStart(START_SYMM);
  Lal::getSymmetrize(DxMat);
  TimeEnd(END_SYMM);
  // rMessage("DxMat =  ");
  // DxMat.display();
  com.makedX += TimeCal(START_DX,END_DX);
  com.symmetriseDx += TimeCal(START_SYMM,END_SYMM);
  // rMessage("DxMat = ");
  // DxMat.display();
}


bool Newton::Mehrotra(Newton::WHICH_DIRECTION direction,
		      int m,
		      InputData& inputData,
		      Chordal& chordal,
		      Solutions& currentPt,
		      Residuals& currentRes,
		      AverageComplementarity& mu,
		      DirectionParameter& beta,
		      Switch& reduction,
		      Phase& phase,
		      WorkVariables& work,
		      ComputeTime& com,
		      FILE* Display, FILE* fpOut)
{
  //   rMessage("xMat, yVec, zMat =  ");
  //   currentPt.xMat.display();
  //   currentPt.yVec.display();
  //   currentPt.zMat.display();
  Make_gVec(direction, inputData, currentPt, currentRes,
	    mu, beta, phase, work, com);
  //  gVec.display();

  if (direction == PREDICTOR) {
    Make_bMat(inputData, currentPt, work, com);

    #if 0
    static int iter = 0;
    char filename[100];
    sprintf(filename,"B%02d.m",iter);
    iter++;
    chordal.writeSchur(filename);
    #endif
    // bMat.display();
    // display_sparse_bMat();
    // display_index();
  }

  bool ret = compute_DyVec(direction,
			   m, inputData, chordal,
			   currentPt, work, com, Display, fpOut);
  if (ret == SDPA_FAILURE) {
    return SDPA_FAILURE;
  }
  //  rMessage("cholesky factorization =  ");
  //  sparse_bMat.display();

  TimeStart(START5);

  compute_DzMat(inputData, currentRes, phase, com);
  compute_DxMat(currentPt, work, com);

  TimeEnd(END5);
  com.makedXdZ += TimeCal(START5,END5);

  // rMessage("DxMat, DyVec, DzMat =  ");
  //   DxMat.display();
  //   DyVec.display();
  //   DzMat.display();

  return true;
}

void Newton::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  fprintf(fpout,"rNewton.DxMat = \n");
  DxMat.display(fpout);
  fprintf(fpout,"rNewton.DyVec = \n");
  DyVec.display(fpout);
  fprintf(fpout,"rNewton.DzMat = \n");
  DzMat.display(fpout);
}

void Newton::display_index(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  printf("display_index: %d %d %d\n",SDP_nBlock,SOCP_nBlock,LP_nBlock);

  for (int l=0; l<SDP_nBlock; l++){
    printf("SDP:%dth block\n",l);
    for (int k=0; k<SDP_number[l]; k++){
      printf("SDP(i=%d,ib=%d; j=%d,jb=%d) for target = %d\n",
	     SDP_constraint1[l][k],SDP_blockIndex1[l][k],
	     SDP_constraint2[l][k],SDP_blockIndex2[l][k], 
	     SDP_location_sparse_bMat[l][k]);
    }
  }

  for (int l=0; l<SOCP_nBlock; l++){
    printf("SOCP:%dth block\n",l);
    for (int k=0; k<SOCP_number[l]; k++){
      printf("SOCP(i=%d,ib=%d; j=%d,jb=%d) for target = %d\n",
	     SOCP_constraint1[l][k],SOCP_blockIndex1[l][k],
	     SOCP_constraint2[l][k],SOCP_blockIndex2[l][k], 
	     SOCP_location_sparse_bMat[l][k]);
    }
  }

  for (int l=0; l<LP_nBlock; l++){
    printf("LP:%dth block\n",l);
    for (int k=0; k<LP_number[l]; k++){
      printf("LP(i=%d,ib=%d; j=%d,jb=%d) for target = %d\n",
	     LP_constraint1[l][k],LP_blockIndex1[l][k],
	     LP_constraint2[l][k],LP_blockIndex2[l][k], 
	     LP_location_sparse_bMat[l][k]);
    }

  }

}

void Newton::display_sparse_bMat(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"{\n");
  for (int id=0; id<MpiSt::nprocs; ++id) {
    if (id == MpiSt::iam) {
      for (int index=mySchurStart; index<mySchurEnd; ++index) {
	int i        = sparse_bMat.row_index[index];
	int j        = sparse_bMat.column_index[index];
	double value = sparse_bMat.sp_ele[index];
	fprintf(fpout,"val[%d,%d] = %e on CPU %d\n", i,j,value, id);
      }
    }
  }
  fprintf(fpout,"}\n");
}

void Newton::display_pB()
{
  MpiSt::barrier();
  int ISTDOUT = 6;
  double* tmpwork;
  NewArray(tmpwork,double,m);
  printf("pB = zeros(%d,%d)",m,m);
  pdlaprnt_(&m,&m,pB,&IONE,&IONE,descB,&IZERO,&IZERO,
	    (char*)"pB", &ISTDOUT,tmpwork,2);
  #if 0
  if (MpiSt::iam == 0) {
    printf("pB = tril(pB)");
    printf("pB = pB + pB' - diag(diag(pB))");
  }
  #endif
  DeleteArray(tmpwork);
  MpiSt::barrier();
}

void Newton::display_pB2()
{
  MpiSt::barrier();
  int ISTDOUT = 6;
  double* tmpwork;
  NewArray(tmpwork,double,maxmp2);
  printf("pB = zeros(%d,%d)",m,m);
  pdlaprnt_(&m,&m,pB2,&IONE,&IONE,descB2,&IZERO,&IZERO,
	    (char*)"pB2", &ISTDOUT,tmpwork,3);
  #if 0
  if (MpiSt::iam == 0) {
    printf("pB = tril(pB)");
    printf("pB = pB + pB' - diag(diag(pB))");
  }
  #endif
  DeleteArray(tmpwork);
  MpiSt::barrier();
}

void Newton::display_pg()
{
  MpiSt::barrier();
  int ISTDOUT = 6;
  double* tmpwork;
  NewArray(tmpwork,double,m);
  pdlaprnt_(&m,&IONE,pg,&IONE,&IONE,descg,&IZERO,&IZERO,
	    (char*)"pg", &ISTDOUT,tmpwork,2);
  DeleteArray(tmpwork);
  MpiSt::barrier();
}

void Newton::display_pg2()
{
  int ISTDOUT = 6;
  double* tmpwork;
  MpiSt::barrier();
  NewArray(tmpwork,double,maxmp2);
  MpiSt::barrier();
  pdlaprnt_(&m,&IONE,pg2,&IONE,&IONE,descg2,&IZERO,&IZERO,
	    (char*)"pg2", &ISTDOUT,tmpwork,3);
  MpiSt::barrier();
  DeleteArray(tmpwork);
  MpiSt::barrier();
}

} // end of namespace 'sdpa'

