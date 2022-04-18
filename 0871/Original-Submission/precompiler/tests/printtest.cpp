#include <stdio.h>

float square( float a){
  return a*a;
}

int main(){
  float a=2.0;
  printf("a=%f,square(3)=%f\n", a, square(3.) );
  return 0;

}

