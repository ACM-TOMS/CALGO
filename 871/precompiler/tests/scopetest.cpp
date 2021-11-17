int f(int a){
  return a+1;
}


int main(){
  int a=2,b=1,c=3;
  if(a+b == 3){
    char c;
    c=32; //c is a char here!
  }
  else{
    c=32; //c is a long int
  }
  return 0;
}

