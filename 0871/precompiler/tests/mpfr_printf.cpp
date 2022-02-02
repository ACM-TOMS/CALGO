#include <stdio.h>


int main(){
  float a[4][4]={
    1,2,3,4,
    5,6,7,8,
    9,10,11,12,
    13,14,15,16
  };
  
  a[1][1]=0.0;
  
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      printf("a[%i,%i]=%f\n",i,j,a[i][j]);
    }
  }
}


