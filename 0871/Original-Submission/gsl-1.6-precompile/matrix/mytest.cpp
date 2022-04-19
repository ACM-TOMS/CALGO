#include "gsl_matrix.h"

int main(){
  MpIeee::fpEnv.setRadix(2);

  //using gsl matrix alloc functions
  gsl_matrix* M=gsl_matrix_alloc( 3,3 ); //allocate 3 by 3
  cout << "allocated M(3,3)"<<endl;


  gsl_matrix_set_identity( M ); //make it the identity matrix
  cout << "set M to identity" << endl;

//  gsl_matrix_set_zero( M ); //make zero matrix, also Abort :(
//  cout << "set M to zeros" << endl;
  
  gsl_matrix_set( M, 1, 1, MpIeee("7") );
  
  cout << "showing M="<<endl;  
  for( int row=0;row<3;row++){
    for( int col = 0; col<3; col++){
      cout << gsl_matrix_get( M, row, col ) << "\t" << flush;
    }
    cout << endl;
  }  
  return 0;
}


