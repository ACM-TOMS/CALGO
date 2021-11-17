/*--------------------------------------------------
  sdpa_mpicopy.cpp
--------------------------------------------------*/

#include "sdpa_mpist.h"
#include "sdpa_mpicopy.h"
#include "sdpa_scalapack.h"
#include "sdpa_linear.h"

namespace sdpa {

void MpiCopy::allSendRecieveC(int length, char* source)
{
  if (MpiSt::iam == 0) {
    Cigebs2d(MpiSt::ictxt,(char*)"All",(char*)" ",length,1,
	     (char*)source,length*sizeof(char)/sizeof(int));
  }
  else {
    Cigebr2d(MpiSt::ictxt,(char*)"All",(char*)" ",length,1,
	     (char*)source,length*sizeof(char)/sizeof(int),0,0);
  }
}

void MpiCopy::allSendRecieveI(int length, int* source)
{
  if (MpiSt::iam == 0) {
    Cigebs2d(MpiSt::ictxt,(char*)"All",(char*)" ",length,1,
	     (char*)source,length);
  }
  else {
    Cigebr2d(MpiSt::ictxt,(char*)"All",(char*)" ",length,1,
	     (char*)source,length,0,0);
  }
}

void MpiCopy::allSendRecieveD(int length, double* source)
{
  if (MpiSt::iam == 0) {
    Cdgebs2d(MpiSt::ictxt,(char*)"All",(char*)" ",length,1,
	     (char*)source,length);
  }
  else {
    Cdgebr2d(MpiSt::ictxt,(char*)"All",(char*)" ",length,1,
	     (char*)source,length,0,0);
  }
}

void MpiCopy::allAccumulateD(int length, double* source)
{
  Cdgsum2d(MpiSt::ictxt, (char*)"All", (char*)" ",
     length, 1, (char*)source , length, 0, 0);
}

void MpiCopy::allAccumulateDirstibuteD(int length, double* source)
{
  Cdgsum2d(MpiSt::ictxt, (char*)"All", (char*)" ",
     length, 1, (char*)source , length, 0, 0);
  allSendRecieveD(length,source);
}

void MpiCopy::sendToHostI(int length, int* source)
{
  Cigesd2d(MpiSt::ictxt, length, 1, (char*)source, length, 0, 0);
}

void MpiCopy::sendToHostD(int length, double* source)
{
  Cdgesd2d(MpiSt::ictxt, length, 1, (char*)source, length, 0, 0);
}

void MpiCopy::receiveAtHostI(int from, int length, int* source)
{
  // Note : if ictxt is row-wise, receive from 'from,0'
  //        but column-wisee, receive from '0,from'
  // Cigerv2d(MpiSt::ictxt, length, 1, (char*)source, length, from, 0);
  Cigerv2d(MpiSt::ictxt, length, 1, (char*)source, length, 0, from);
}

void MpiCopy::receiveAtHostD(int from, int length, double* source)
{
  // Cdgerv2d(MpiSt::ictxt, length, 1, (char*)source, length, from, 0);
  Cdgerv2d(MpiSt::ictxt, length, 1, (char*)source, length, 0, from);
}

void MpiCopy::copy(Parameter& param)
{
  allSendRecieveI(1,&param.maxIteration);
  double array[9];
  if (MpiSt::iam == 0) {
    array[0] = param.epsilonStar;
    array[1] = param.lambdaStar;
    array[2] = param.omegaStar;
    array[3] = param.lowerBound;
    array[4] = param.upperBound;
    array[5] = param.betaStar;
    array[6] = param.betaBar;
    array[7] = param.gammaStar;
    array[8] = param.epsilonDash;
  }
  allSendRecieveD(9,array);
  if (MpiSt::iam != 0) {
    param.epsilonStar  = array[0];
    param.lambdaStar   = array[1];
    param.omegaStar    = array[2];
    param.lowerBound   = array[3];
    param.upperBound   = array[4];
    param.betaStar     = array[5];
    param.betaBar      = array[6];
    param.gammaStar    = array[7];
    param.epsilonDash  = array[8];
  }
  allSendRecieveC(PRINT_DEFAULT_LENGTH,param.xPrint);
  allSendRecieveC(PRINT_DEFAULT_LENGTH,param.XPrint);
  allSendRecieveC(PRINT_DEFAULT_LENGTH,param.YPrint);
  allSendRecieveC(PRINT_DEFAULT_LENGTH,param.infPrint);
  // param.display(stdout);
}

void MpiCopy::copy(int& m, int& nBlock, BlockStruct& bs)
{
  int array[3];
  if (MpiSt::iam == 0) {
    array[0] = m;
    array[1] = nBlock;
  }
  allSendRecieveI(2,array);
  if (MpiSt::iam != 0) {
    m = array[0];
    nBlock = array[1];
    bs.initialize(nBlock);
  }
  allSendRecieveI(nBlock,bs.blockStruct);
  allSendRecieveI(nBlock,(int*)bs.blockType);
  if (MpiSt::iam != 0) {
    // trick part for MPI so that bs.makeInternalStructure works correctly
    for (int l=0; l<nBlock; ++l) {
      if (bs.blockType[l] == BlockStruct::btLP) {
	bs.blockStruct[l] = - bs.blockStruct[l];
      }
    }

    bs.makeInternalStructure();
  }
}

void MpiCopy::copy(InputData& inputData, int m, BlockStruct& bs)
{
  if (MpiSt::iam != 0) {
    inputData.initialize(bs);
  }
  #if 0
  for (int np=0; np<MpiSt::nprocs; ++np) {
    MpiSt::barrier();
    if (MpiSt::iam == np) {
      rMessage("blockStruct");
      bs.display();
    }
    MpiSt::barrier();
  }
  #endif
    
  MpiCopy::copy(inputData.b, m);
  MpiCopy::copy(inputData.C);
  // inputData.C.display();
  if (MpiSt::iam != 0) {
    NewArray(inputData.A, SparseLinearSpace, m);
  }
  for (int k=0; k<m; ++k) {
    MpiCopy::copy(inputData.A[k]);
    // inputData.A[k].display();
  }
  // inputData.display();
}

void MpiCopy::copy(Vector& yVec, int m)
{
  if (MpiSt::iam != 0) {
    yVec.initialize(m);
  }
  allSendRecieveD(m,yVec.ele);
}

void MpiCopy::copy(SparseLinearSpace& C)
{
  int array[3];
  if (MpiSt::iam == 0) {
    array[0] = C.SDP_sp_nBlock;
    array[1] = C.SOCP_sp_nBlock;
    array[2] = C.LP_sp_nBlock;
  }
  allSendRecieveI(3,array);
  if (MpiSt::iam != 0) {
    C.SDP_sp_nBlock  = array[0];
    C.SOCP_sp_nBlock = array[1];
    C.LP_sp_nBlock   = array[2];
  }
  if (C.SDP_sp_nBlock > 0){
    if (MpiSt::iam != 0) {
      NewArray(C.SDP_sp_index,int,C.SDP_sp_nBlock);
      NewArray(C.SDP_sp_block,SparseMatrix,C.SDP_sp_nBlock);
    }
    allSendRecieveI(C.SDP_sp_nBlock,C.SDP_sp_index);
    for (int l=0; l<C.SDP_sp_nBlock; ++l) {
      MpiCopy::copy(C.SDP_sp_block[l]);
    }
  }

  // for SOCP
#if 0
  if (C.SOCP_sp_nBlock > 0){
    if (MpiSt::iam != 0) {
      NewArray(C.SOCP_sp_index,int,C.SOCP_sp_nBlock);
      NewArray(C.SOCP_sp_block,SparseMatrix,C.SOCP_sp_nBlock);
    }
    allSendRecieveI(C.SOCP_sp_nBlock,C.SOCP_sp_index);
    for (int l=0; l<C.SOCP_sp_nBlock, ++l) {
      MpiCopy::copy(C.SOCP_sp_block[l]);
    }
  }
#endif
  // for LP
  if (C.LP_sp_nBlock > 0){
    if (MpiSt::iam != 0) {
      NewArray(C.LP_sp_index,int,C.LP_sp_nBlock);
      NewArray(C.LP_sp_block,double,C.LP_sp_nBlock);
    }
    allSendRecieveI(C.LP_sp_nBlock,C.LP_sp_index);
    allSendRecieveD(C.LP_sp_nBlock,C.LP_sp_block);
  }
}

void MpiCopy::copy(SparseMatrix& C)
{
  // At this step, C should be sparse
  // because this step is before changeToDense in initializeSolve
  int array[5];
  if (MpiSt::iam == 0) {
    if (C.type == SparseMatrix::DENSE) {
      rError("We should consider DENSE case");
    }
    array[0] = C.nRow;
    array[1] = C.nCol;
    array[2] = C.NonZeroNumber;
    array[3] = C.NonZeroCount;
    array[4] = C.NonZeroEffect;
  }
  allSendRecieveI(5,array);
  if (MpiSt::iam != 0) {
    C.initialize(array[0],array[1], SparseMatrix::SPARSE, array[2]);
    C.NonZeroCount  = array[3];
    C.NonZeroEffect = array[4];
  }
  #if DATA_CAPSULE
  for (int i=0; i<C.NonZeroCount; ++i) {
    allSendRecieveI(1,&C.DataS[i].vRow);
    allSendRecieveI(1,&C.DataS[i].vCol);
    allSendRecieveD(1,&C.DataS[i].vEle);
  }    
  #else
  allSendRecieveI(C.NonZeroCount,C.row_index);
  allSendRecieveI(C.NonZeroCount,C.column_index);
  allSendRecieveD(C.NonZeroCount,C.sp_ele);
  #endif
}

void MpiCopy::copy(DenseLinearSpace& xMat)
{
  int array[3];
  if (MpiSt::iam == 0) {
    array[0] = xMat.SDP_nBlock;
    array[1] = xMat.SOCP_nBlock;
    array[2] = xMat.LP_nBlock;
  }
  allSendRecieveI(3,array);
  if (MpiSt::iam != 0) {
    xMat.SDP_nBlock  = array[0];
    xMat.SOCP_nBlock = array[1];
    xMat.LP_nBlock   = array[2];
  }
  
  if (xMat.SDP_nBlock > 0){
    if (MpiSt::iam != 0) {
      NewArray(xMat.SDP_block,DenseMatrix,xMat.SDP_nBlock);
    }
    for (int l=0; l<xMat.SDP_nBlock; ++l) {
      MpiCopy::copy(xMat.SDP_block[l]);
    }
  }

  // for SOCP
#if 0
  if (xMat.SOCP_nBlock > 0){
    if (MpiSt::iam != 0) {
      NewArray(xMat.SOCP_block,DenseMatrix,xMat.SOCP_nBlock);
    }
    for (int l=0; l<xMat.SOCP_nBlock; ++l) {
      MpiCopy::copy(xMat.SOCP_block[l]);
    }
  }

#endif

  // for LP
  if (xMat.LP_nBlock > 0){
    if (MpiSt::iam != 0) {
      NewArray(xMat.LP_block,double,xMat.LP_nBlock);
    }
    allSendRecieveD(xMat.LP_nBlock,xMat.LP_block);
  }
}

void MpiCopy::copy(DenseMatrix& xMat)
{
  int array[2];
  if (MpiSt::iam == 0) {
    if (xMat.type == DenseMatrix::COMPLETION) {
      rError("We should consider DENSE case");
    }
    array[0] = xMat.nRow;
    array[1] = xMat.nCol;
  }
  allSendRecieveI(2,array);
  if (MpiSt::iam != 0) {
    xMat.initialize(array[0],array[1], DenseMatrix::DENSE);
  }
  allSendRecieveD(xMat.nRow * xMat.nCol, xMat.de_ele);
}

} // end of namespace 'sdpa'
