/*-----------------------------------------
  sdpa_chordal.h
-----------------------------------------*/

#ifndef __sdpa_chordal_h__
#define __sdpa_chordal_h__

#include "sdpa_dataset.h"
#include <dmumps_c.h>
  
#define SELECT_MUMPS_BEST  7 // MUMPS selects automatically when 7
#define SELECT_DENSE      -1 // This value must be minus

namespace sdpa {

class Chordal {
public:

  // condition of sparse computation
  // m_threshold < mDim, 
  // b_threshold < nBlock, 
  // aggregate_threshold >= aggrigated sparsity ratio
  // extend_threshold    >= extended sparsity ratio
  int m_threshold;
  int b_threshold;
  double aggregate_threshold;
  double extend_threshold;

  int   best;
/* indicates the used ordering method */
  /* -1: dense  computation */
  /*  7: sparse computation by MUMPS */

  SparseMatrix* sparse_bMat_ptr;
  DMUMPS_STRUC_C mumps_id;
  bool mumps_usage;

  int mySchurStart, mySchurEnd, mySchurLength;

  Chordal(void);
  ~Chordal();
  void initialize(SparseMatrix* sparse_bMat_ptr);
  void terminate();

  // merge array1 to array2
  void mergeArray(int na1, int* array1, int na2, int* array2);
  void catArray(int na1, int* array1, int na2, int* array2);
  void slimArray(int i, int length, int* array, int& slimedLength);
 
  void makeGraph(InputData& inputData, int m);
  
  void ordering_bMat(int m, int nBlock,
                     InputData& inputData, FILE* Display,
                     FILE* fpOut);
  double analysisAndcountLowerNonZero(int m);
  void setSchurIndices(int mySchurStart, int mySchurEnd,
		       int mySchurLength);
  bool factorizeSchur(int m, int* diagonalIndex,
		      FILE* Display, FILE* fpOut);
  bool solveSchur(Vector& rhs);

  void writeSchur(char* filename);
  void writeSchurPart(char* filename);
};

} // end of namespace 'sdpa'

#endif // __sdpa_chordal_h__
