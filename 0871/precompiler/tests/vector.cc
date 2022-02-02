#include <vector>
#include <stdio.h>

int main(){
  int matrix_size=4;
  
  vector< vector<float> > A;
  vector< float > B;
  
  A.resize(matrix_size+1);
  for(int i=0;i<=matrix_size;i++) A[i].resize(matrix_size+1);
  
  B.resize(matrix_size+1);



  B[1]=-2.0;    B[2]=-7.0;    B[3]=20.0;    B[4]=13.0; 
  printf("array B stored\n");

  A[1][1]=2.0;  A[1][2]=3.0;  A[1][3]=-4.0; A[1][4]=1.0;
  A[2][1]=1.0;  A[2][2]=-1.0; A[2][3]=0;    A[2][4]=-2.0; 
  A[3][1]=3.0;  A[3][2]=3.0;  A[3][3]=4.0;  A[3][4]=3.0;
  A[4][1]=4.0;  A[4][2]=1.0;  A[4][3]=0;    A[4][4]=4.0;
  printf("array A stored\n");


  
  return 0;
}