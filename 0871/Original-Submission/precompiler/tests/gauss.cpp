#include <vector>
#include <iostream>

using namespace std;

typedef vector< vector <double> > matrix;
typedef vector< double > col;
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
  cout<<"Matrix A forward elimination"<<endl;
  
  col pivots;
  pivots.clear();
  
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
      pivots.push_back(xmult);
      A[i][k]=0; //spil to 0
      for (j = k+1 ; j <= matrix_size; j++){
        A[i][j] -= xmult * A[k][j];
      }
      //B[i] -= xmult * B[k];      
    }
  }
  
  //X[matrix_size] = B[matrix_size] / A[matrix_size][matrix_size];

  cout<<"===================="<<endl;  
  cout<<"pivots ="<<endl;
  show(pivots);
  cout<<"===================="<<endl;
  
}






float mabs(float a){
  if(a<0) a=-a;
  return a;
}


void swap(icol& L,int i, int j){
  int temp=L[i];
  L[i]=L[j];
  L[j]=temp;
}


void swap(matrix& A, col& B, int i, int j){
  float temp;
  for(int col=1;col<A.size();col++){
    temp=A[i][col];
    A[i][col]=A[j][col];
    A[j][col]=temp;
    
  }
  
  //swap B also!
  temp=B[i];
  B[i]=B[j];
  B[j]=temp;
}




//partial pivoting (no lookup vector used)
//just swap complete rows in A
void forward_elim_pp (matrix& A, col& B, col& X ) {
  cout<<"Matrix A forward elimination pp"<<endl;
  int i, j, k, mpos;
  float xmult,tmp,max;
  
  int matrix_size = A.size()-1;

  for (k = 1; k < matrix_size; k++) {
    max=A[k][k];
    mpos=k;
    for(int l=k+1;l<= matrix_size;l++){
      tmp=A[l][k];
      tmp=mabs(tmp);
      if( tmp > max ){
        max=tmp;
        mpos=l;
      }
    }
    
    if(k!=mpos) swap(A,B,k,mpos);
    
    if(max == 0 ) {
        cerr<<"Matrix is singular!"<<endl;
        exit(0);
    }
    
    for (i = k+1; i <= matrix_size; i++) {
     
      xmult = A[i][k] / A[k][k];
      A[i][k]=0; //spil to 0
      for (j = k+1 ; j <= matrix_size; j++){
        A[i][j] -= xmult * A[k][j];
      }
      B[i] -= xmult * B[k];      
    }
  }
  
  
}


void back_subst (matrix& A, col& B, col& X) {
  int i, j;
  float sum;

  int matrix_size=A.size()-1;
  X[matrix_size] = B[matrix_size] / A[matrix_size][matrix_size];
  
  for (i = A.size()-1; i >= 1; i--) {
    sum = B[i];
    for (j = i+1; j < A.size(); j++){
      sum -= A[i][j] * X[j];
    }
    X[i] = sum / A[i][i];
  }

}




//with partial pivoting with lookup vector...
void forward_elim(matrix& A, col& B, col& X , icol& L ) {
  cout<<"Matrix A forward elimination pp with lookup vector"<<endl;
  for (int i=1;i<L.size();i++) L[i]=i;
  
  int i, j, k,l,mpos, row;
  int tmpl,tmpk;
  float xmult,max,tmp;
  

  int matrix_size = A.size()-1;
  
  for (k = 1; k < matrix_size; k++) {

    max=A[L[k]][k ];
    mpos=k;
    for(int l=k+1;l<= matrix_size;l++){
      tmpl=L[l];
      tmp=A[tmpl][k];
      tmp=mabs(tmp);
      if( tmp > max ){
        max=tmp;
        mpos=l;
      }
    }
    
    if(k!=mpos) swap(L,k,mpos);
    
    if( max == 0.0 ){
      cerr << "Matrix is singular!" << endl;
      exit(0);
    }
        
        
    for (i = k+1; i <= matrix_size; i++) { // i = row index in L
      tmpk=L[k];      //spil pos
      row=L[i];       //row to change
      xmult = A[row][k] / A[tmpk][k];
      A[row][k]=0; //set 0 on spil
      for (j = k+1 ; j <= matrix_size; j++){ //col
        A[row][j] -= xmult * A[tmpk][j];
      }
      
      B[row] -= xmult * B[tmpk];      
    }
  }


}


// with partial pivoting with lookup
void back_subst (matrix& A, col& B, col& X, icol& L) {
  int i, j,tmpi;
  float sum;

  tmpi=L[A.size()-1];
  X[tmpi] = B[tmpi] / A[tmpi][tmpi];
  
  for (i = A.size()-1; i >= 1; i--) {
    tmpi= L[i];
    sum = B[tmpi];
    for (j = i+1; j < A.size(); j++)
      sum -= A[tmpi][j] * X[L[j]];
    X[tmpi] = sum / A[tmpi][i];
  }

}





 //other tests ...
 
 /* matrix_size=3;
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
  */


  /*
  //a singular matrix
  A[1][1]=1.0;  A[1][2]=3.0;  A[1][3]=1.0;  A[1][4]=1.0;
  A[2][1]=1.0;  A[2][2]=3.0; A[2][3]=-1.0;    A[2][4]=2.0; 
  
  A[3][1]=1.0;  A[3][2]=3.0;  A[3][3]=8.0;  A[3][4]=3.0;
  A[4][1]=1.0;  A[4][2]=3.0;  A[4][3]=5.0;    A[4][4]=4.0;
    
  B[1]=1.0;    B[2]=2.0;    B[3]=3.0;    B[4]=4.0; 
  */

/*
  A[1][1]=2.0;  A[1][2]=4.0;  A[1][3]=1.0;  A[1][4]=4.0;
  A[2][1]=7.0;  A[2][2]=3.0;  A[2][3]=-5.0; A[2][4]=8.0; 
  A[3][1]=9.0;  A[3][2]=6.0;  A[3][3]=2.0;  A[3][4]=1.0;
  A[4][1]=1.0;  A[4][2]=-3.0; A[4][3]=4.0;  A[4][4]=6.0;
    
  B[1]=15.0;    B[2]=13.0;    B[3]=6.0;     B[4]=2.0; 
*/


void init(matrix& A, col& B){
  A[1][1]=2.0;  A[1][2]=3.0;  A[1][3]=-4.0; A[1][4]=1.0;
  A[2][1]=1.0;  A[2][2]=-1.0; A[2][3]=0;    A[2][4]=-2.0; 
  A[3][1]=3.0;  A[3][2]=3.0;  A[3][3]=4.0;  A[3][4]=3.0;
  A[4][1]=4.0;  A[4][2]=1.0;  A[4][3]=0;    A[4][4]=4.0;

  B[1]=-2.0;    B[2]=-7.0;    B[3]=20.0;    B[4]=13.0; 

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






  //=============test with partial pivoting and lookup =====================
  init(A,B);  
  
  forward_elim( A, B, X , L);   
  show(A);
  
  
  back_subst ( A, B, X , L );
  cout<<"Matrix A after back subst"<<endl;
  show(A);
  
  cout<<" Vector X = "<<endl;  
  show(X,L); //with right order using offsets L
  //show(X); 
  
  cout<<"L="<<endl;
  show(L);


  //====================test with partial pivoting ======================== 
  init(A,B);
  
  forward_elim_pp( A, B, X );   
  show(A);
  
  
  back_subst ( A, B, X );
  cout<<"Matrix A after back subst"<<endl;
  show(A);
  
  cout<<" Vector X = "<<endl;  
  show(X); 
  
  
  //========================= no partial pivoting ==========================
  cout<<endl<<endl;
  init(A,B);
  show(A);cout<<endl;
  
  forward_elim( A, B, X );   
  show(A);
  
  
  back_subst ( A, B, X );
  cout<<"Matrix A after back subst"<<endl;
  show(A);
  
  cout<<" Vector X = "<<endl;  
  show(X); 
  
  
  return 0;
}


