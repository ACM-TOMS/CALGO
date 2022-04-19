#include <stdio.h>
int f=1;
int pow(int a,unsigned int n){
  if(n<0) return 0;
  for(int i=1;i<n;i=i+1){
    a*=a;
  }
  return a+1;
}

int main(){
  int a=3,b=1,c=3;
  if(a+b == 3){
    char c;
    c=32; //c is a char here!
  }
  else{
    c=32+pow(a,b); //c is a long int
  }
  printf("c=%i\n",c);
  return 0;
}

