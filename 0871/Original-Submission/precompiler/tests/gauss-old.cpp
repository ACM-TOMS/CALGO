#include <vector>
#include <iostream>

using namespace std;

typedef vector< vector <float> > matrix;
typedef vector< float > col;

typedef vector< int > icol;





void show(matrix& A){
  for(int i=1;i<A.size();i++){
    for(int j=1;j<A.size();j++){
      cout<<"\t"<<A[i][j];
    }
    cout<<endl;
  }
}

void show(col& L){
  for(int i=1;i<L.size();i++){
      cout<<"\t"<<L[i]<<endl;
  }
}

void show(icol& L){
  for(int i=1;i<L.size();i++){
      cout<<"\t"<<L[i]<<endl;
  }
}


void show(col& X ,icol& L){
  for( int i=1; i<X.size(); i++ )
    cout<<"\t"<<X[L[i]]<<endl;
  cout<<endl;
  
}




//no partial pivoting...
void forward_elim (matrix& A, col& B, col& X ) {
  int i, j, k;
  float xmult;
  
  int matrix_size = A.size()-1;

  for (k = 1; k < matrix_size; k++) {
    for (i = k+1; i <= matrix_size; i++) {
      if(A[k][k] == 0 ) {
        cerr<<"Matrix is singular!"<<endl;
        exit(0);
      }
      xmult = A[i][k] / A[k][k];
      A[i][k]=0; //spil to 0
      for (j = k+1 ; j <= matrix_size; j++){
        A[i][j] -= xmult * A[k][j];
      }
      B[i] -= xmult * B[k];      
    }
  }
  
  X[matrix_size] = B[matrix_size] / A[matrix_size][matrix_size];
}

void back_subst (matrix& A, col& B, col& X) {
  int i, j;

  float sum;

  for (i = A.size()-2; i >= 1; i--) {
    sum = B[i];
    for (j = i+1; j < A.size(); j++)
      sum -= A[i][j] * X[j];
    X[i] = sum / A[i][i];
  }

}



float abs(float a){
  if(a<0) a=-a;
  return a;
}

void swap(icol& L,int i, int j){
  int temp=L[i];
  L[i]=L[j];
  L[j]=temp;
}

//with partial pivoting
void forward_elim (matrix& A, col& B, col& X , icol& L ) {
  for (int i=1;i<L.size();i++) L[i]=i;
  
  int i, j, k,l,mpos, row;
  int tmpl,tmpk;
  float xmult,max=0,tmp;
  

  int matrix_size = A.size()-1;
  
  for (k = 1; k < matrix_size; k++) {

    cout<<"A="<<endl;show(A);
    max=A[1][k];
    mpos=1;
    for(int l=1;l<= matrix_size;l++){
      tmpl=L[l];
      tmp=A[tmpl][k];
      tmp=abs(tmp);
      if( tmp > max ){
        max=tmp;
        mpos=tmpl;
      }
    }
    cout<<" -> max="<<max<<", mpos="<<mpos<<endl;
    if(k!=mpos) swap(L,k,mpos);
    
    if( max == 0.0 ){
      cerr << "Matrix is singular!" << endl;
      exit(0);
    }
        
    for (i = k+1; i <= matrix_size; i++) {
      tmpk=L[k];
      
      row=L[i];
      xmult = A[row][k] / A[tmpk][k];
      A[row][k]=0; //spil to 0
      for (j = k+1 ; j <= matrix_size; j++){
        A[row][j] -= xmult * A[tmpk][j];
      }
      B[row] -= xmult * B[tmpk];      
    }
  }

  X[matrix_size] = B[matrix_size] / A[matrix_size][matrix_size];
}


// with partial pivoting
void back_subst (matrix& A, col& B, col& X, icol& L) {
  int i, j,tmpi;
  float sum;

  for (i = A.size()-2; i >= 1; i--) {
    tmpi= L[i];
    sum = B[tmpi];
    for (j = i+1; j < A.size(); j++)
      sum -= A[tmpi][j] * X[L[j]];
    X[tmpi] = sum / A[tmpi][i];
  }

}



int main(){

  int i,j,matrix_size;
  matrix_size=4;
  
  matrix A;
  col B,X;
  icol L;
  
  //allocate matrix and vectors
  A.resize(matrix_size+1);
  for(int i=0;i<=matrix_size;i++) A[i].resize(matrix_size+1);
  
  B.resize(matrix_size+1);
  X.resize(matrix_size+1);
  L.resize(matrix_size+1);


/*
  A[1][1]=2.0;  A[1][2]=3.0;  A[1][3]=-4.0; A[1][4]=1.0;
  A[2][1]=1.0;  A[2][2]=-1.0; A[2][3]=0;    A[2][4]=-2.0; 
  A[3][1]=3.0;  A[3][2]=3.0;  A[3][3]=4.0;  A[3][4]=3.0;
  A[4][1]=4.0;  A[4][2]=1.0;  A[4][3]=0;    A[4][4]=4.0;

  B[1]=-2.0;    B[2]=-7.0;    B[3]=20.0;    B[4]=13.0; 
*/

  A[1][1]=2.0;  A[1][2]=4.0;  A[1][3]=1.0;  A[1][4]=4.0;
  A[2][1]=7.0;  A[2][2]=3.0;  A[2][3]=-5.0; A[2][4]=8.0; 
  A[3][1]=9.0;  A[3][2]=6.0;  A[3][3]=2.0;  A[3][4]=1.0;
  A[4][1]=1.0;  A[4][2]=-3.0; A[4][3]=4.0;  A[4][4]=6.0;
    
  B[1]=15.0;    B[2]=13.0;    B[3]=6.0;     B[4]=2.0; 
 

/*
  A[1][1]=1.0;  A[1][2]=1.0;  A[1][3]=1.0;  A[1][4]=1.0;
  A[2][1]=2.0;  A[2][2]=2.0; A[2][3]=2.0;    A[2][4]=2.0; 
  
  A[3][1]=3.0;  A[3][2]=-1.0;  A[3][3]=8.0;  A[3][4]=3.0;
  A[4][1]=1.0;  A[4][2]=4.0;  A[4][3]=5.0;    A[4][4]=4.0;
    
  B[1]=1.0;    B[2]=2.0;    B[3]=3.0;    B[4]=4.0; 
*/
  
  forward_elim( A, B, X ,L );   
  cout<<"Matrix A after forward elimination"<<endl;
  show(A);
  
  
  back_subst ( A, B, X ,L );
  cout<<"Matrix A after back subst"<<endl;
  show(A);
  
  cout<<" Vector X = "<<endl;  
  show(X,L);  
  
  
  matrix_size=3;
  matrix A2;
  col B2,X2;
  icol L2;
  A2.resize(matrix_size+1);
  for(int i=0;i<=matrix_size;i++) A2[i].resize(matrix_size+1);
  
  B2.resize(matrix_size+1);
  X2.resize(matrix_size+1);
  L2.resize(matrix_size+1);

  A2[1][1]=1.0;  A2[1][2]=1.0;  A2[1][3]=2.0;  
  A2[2][1]=2.0;  A2[2][2]=1.0;  A2[2][3]=4.0; 
  A2[3][1]=1.0;  A2[3][2]=1.0;  A2[3][3]=1.0;  
    
  B2[1]=6.0;    B2[2]=2.0;    B2[3]=4.0;  
  
  forward_elim(A2,B2,X2,L2);
  back_subst(A2,B2,X2,L2);
  
  cout<<"A2="<<endl;
  show(A2);
  cout<<"X2="<<endl;
  show(X2,L2);
  cout<<"L2="<<endl;
  show(L2);
  
  return 0;
}

