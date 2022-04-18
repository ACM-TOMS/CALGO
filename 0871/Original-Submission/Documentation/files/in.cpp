#include <stdio.h>

double globD=2.2;

int main(){
  float b;
  for(float a=0;a<10;a=a+0.5){
    b=a*2;    
    printf("a=%f",a);   
  }
  return 0;
}

