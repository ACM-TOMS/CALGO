#include <stdio.h>

float func(float& a){
  float k;
  for(int i=1;i<5;i++)
    for(int j=1;j<3;j++){
      char k=32;
      if(i+j>5) k='y';
      if(k=='y') printf("k=%c, i+j=%i\n",k,i+j);
    }
  k=2*a;
  return 3*a;
}

int main(){
  float b=0.1,c=0;
  double d[3][3]=
    {
      1.0, 2.0, 3,
      4, 5, 6,
      7, 8, 9
    };
  c=3+26;
  b=func(c);
  printf("b=%f\n",b);
  printf("d[1][2]=%f\n",d[1][2]);
  return 0;
}

