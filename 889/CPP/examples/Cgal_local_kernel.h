#ifndef CGAL_LOCAL_KERNEL_H
#define CGAL_LOCAL_KERNEL_H

#include <CGAL/Linear_algebraCd.h>
#include <CGAL/eigen.h>
#include <CGAL/Cartesian_converter.h>

////////////////////////class Cgal_local_kernel/////////////////////
template< class FT >
class Cgal_local_kernel:public CGAL::Cartesian< FT > {
  
 public:
  //2*2 and 3*3 matrices, *, inverse, constructors
  typedef typename CGAL::Linear_algebraCd<FT>::Matrix LKMatrix;
  
  //used for 3*3 matrices
  int sign_of_determinant (LKMatrix M)
  { return CGAL::Linear_algebraCd<FT>::sign_of_determinant (M); }
  
  //used for 3*3 matrices
  LKMatrix inverse (LKMatrix M, FT& D)
  { return CGAL::Linear_algebraCd<FT>::inverse (M, D); }
  
  //must provide constructor, inverse, *
  typedef typename CGAL::Aff_transformation_3< CGAL::Cartesian< FT > > Aff_transformation;
  
  FT Lsqrt( FT x) { return CGAL::sqrt<FT>(x); }
  
  void eigen_symmetric (const FT *mat, 
			const int n, 
			FT *eigen_vectors, 
			FT *eigen_values)
  { CGAL::CGALi::eigen_symmetric<FT>(mat,n,eigen_vectors,eigen_values); }
};

//////////////////////end class Cgal_local_kernel/////////////////////

////////////////////////class Cgal_Kernel_converter/////////////////////
//the two member converters must provide an operator() for the number
//types, points and vectors

template< class K1, class K2 >
  class Cgal_kernel_converters 
{
 public:
 CGAL::Cartesian_converter< K1, K2 > D2L_converter;
 CGAL::Cartesian_converter< K2, K1 > L2D_converter;
};
////////////////////////end class Cgal_Kernel_converter/////////////////////
 
#endif // CGAL_LOCAL_KERNEL_H
