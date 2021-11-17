#include <vector>
#include <iostream>
using namespace std;
typedef vector< vector <float> > matrix;

void init(matrix& A){
  A[1][1]=2.0;  A[1][2]=3.0;  A[1][3]=-4.0; 
  A[2][1]=1.0;  A[2][2]=-1.0; A[2][3]=0;    
  A[3][1]=3.0;  A[3][2]=3.0;  A[3][3]=4.0;  
}

int main(){
  int i,matrix_size;
  matrix_size=3;
  matrix A;
  float b=0;
  A.resize(matrix_size+1);
  for(int i=0;i<=matrix_size;i++) A[i].resize(matrix_size+1);
  init(A); 
  if( A[3][3] > 3 ) b=A[3][3];
  if( b == 4.0  ) cout<<"ok, b="<<b<<endl;
  else cout<<"failed, b="<<b<<endl;
  return 0;
}

