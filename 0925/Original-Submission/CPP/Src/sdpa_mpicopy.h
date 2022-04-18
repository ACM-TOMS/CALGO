/*--------------------------------------------------
  sdpa_mpicopy.h
--------------------------------------------------*/

#ifndef __sdpa_mpicopy_h__
#define __sdpa_mpicopyh__

#include "sdpa_struct.h"
#include "sdpa_io.h"

namespace sdpa {

class MpiCopy 
{
public:
  static void allSendRecieveC(int length, char*   source);
  static void allSendRecieveI(int length, int*    source);
  static void allSendRecieveD(int length, double* source);
  static void allAccumulateD(int length, double* source);
  static void allAccumulateDirstibuteD(int length, double* source);
  static void sendToHostI(int length, int* source);
  static void sendToHostD(int length, double* source);
  static void receiveAtHostI(int from, int length, int* source);
  static void receiveAtHostD(int from, int length, double* source);

  static void copy(Parameter& param);
  static void copy(int& m, int& nBlock, BlockStruct& bs);
  static void copy(InputData& inputData, int m, BlockStruct& bs);
  static void copy(Vector& yVec, int m);
  static void copy(SparseLinearSpace& C);
  static void copy(SparseMatrix& C);
  static void copy(DenseLinearSpace& xMat);
  static void copy(DenseMatrix& xMat);

};

} // end of namespace 'sdpa'

#endif // __sdpa_mpicopy_h__
