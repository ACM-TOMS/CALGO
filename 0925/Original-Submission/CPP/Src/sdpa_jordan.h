
#ifndef __sdpa_jordan_h__
#define __sdpa_jordan_h__

#include "sdpa_struct.h"

namespace sdpa {

class WorkVariables;
  
class Jal
{
public:

  static double trace(DenseLinearSpace& aMat);

  // calculate the minimum eigen value of lMat*xMat*(lMat^T)
  // by Lanczos methods.
  // lMat is lower triangular, xMat is symmetric
  // block size > 20   : Lanczos method
  // block size <= 20  : QR method
  static double getMinEigen(DenseLinearSpace& lMat,
			    DenseLinearSpace& xMat,
			    WorkVariables& work);

  // calculate the minimum eigen value of xMat by QR method.
  static double getMinEigen(DenseLinearSpace& xMat,
			    WorkVariables& work);

  static bool getInvChol(DenseLinearSpace& invCholMat,
			 DenseLinearSpace& aMat,
			 DenseLinearSpace& workMat);

  static bool getInvCholAndInv(DenseLinearSpace& invCholMat,
			       DenseLinearSpace& inverseMat,
			       DenseLinearSpace& aMat,
			       DenseLinearSpace& workMat);

  static bool multiply(DenseLinearSpace& retMat,
		       DenseLinearSpace& aMat,
		       DenseLinearSpace& bMat,
		       double* scalar = NULL);
#if 0
// CAUTION!!! We don't initialize retMat to zero matrix for efficiently.
  static bool multiply(DenseLinearSpace& retMat,
		       SparseLinearSpace& aMat,
		       DenseLinearSpace& bMat,
		       double* scalar = NULL);
// CAUTION!!! We don't initialize retMat to zero matrix for efficiently.
  static bool multiply(DenseLinearSpace& retMat,
		       DenseLinearSpace& aMat,
		       SparseLinearSpace& bMat,
		       double* scalar = NULL);
#endif

  //  retMat = L_{A} B = (A * B + B * A)/2
  static bool jordan_product(DenseLinearSpace& retMat,
			     DenseLinearSpace& aMat,
			     DenseLinearSpace& bMat);
  //  retMat = A * B
  static bool ns_jordan_product(DenseLinearSpace& retMat,
				DenseLinearSpace& aMat,
				DenseLinearSpace& bMat);
  
  //  retMat = P_{A} B = A * B * A 
  static bool jordan_quadratic_product(DenseLinearSpace& retMat,
				       DenseLinearSpace& aMat,
				       DenseLinearSpace& bMat,
				       DenseLinearSpace& work);
  
  //  retMat = Q_{A,C} B = (A * B * C + C * B * A)/2
  static bool jordan_triple_product(DenseLinearSpace& retMat,
				    DenseLinearSpace& aMat,
				    DenseLinearSpace& bMat,
				    DenseLinearSpace& cMat,
				    DenseLinearSpace& work);

  //  retMat = A * B * C
  static bool ns_jordan_triple_product(DenseLinearSpace& retMat,
				       DenseLinearSpace& aMat,
				       DenseLinearSpace& bMat,
				       DenseLinearSpace& cMat,
				       DenseLinearSpace& work);

};

} // end of namespace 'sdpa'

#endif // __sdpa_jordan_h__
