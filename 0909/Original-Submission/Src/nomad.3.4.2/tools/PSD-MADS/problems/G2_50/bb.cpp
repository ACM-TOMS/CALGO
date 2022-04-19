/*--------------*/
/*  problem G2  */
/*--------------*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <sstream>
using namespace std;

#define N 50

int main ( int argc , char ** argv ) {

  long double z = 1e+20 , g1 = 1e+20 , g2 = 1e+20;

  if ( argc != 2 ) {
    cout << g1 << " " << g2 << " " << z << endl;
    return 1;
  }
  
  long double x[N];

  ifstream in ( argv[1] );
  if ( in.fail() ) {
    cout << g1 << " " << g2 << " " << z << endl;
    return 1;
  }
  
  int i;
 
  for ( i = 0 ; i < N ; i++ )
    in >> x[i];

  if ( in.fail() ) {
    cout << g1 << " " << g2 << " " << z << endl;
    in.close();
    return 1;
  }
  
  in.close();

  /*-------------------------------------------------------------------*/

  long double sum1 = 0.0 , sum2 = 0.0 , sum3 = 0.0 , prod1 = 1.0 , prod2 = 1.0;
  
  for ( i = 0 ; i < N ; i++ ) {
    sum1  += pow ( cos(x[i]) , 4 );
    sum2  += x[i];
    sum3  += (i+1)*x[i]*x[i];
    prod1 *= pow ( cos(x[i]) , 2 );
    prod2 *= x[i];
  }


  g1 = -prod2+0.75;
  g2 = sum2 -7.5*N;

  z  = - fabs ( ( sum1 - 2 * prod1 ) / sqrt(sum3) );

  cout << setprecision(16);
  cout << g1 << " " << g2 << " " << z << endl;

//   /*-----------------------------------*/
//   /*           for observer            */
//   /*-----------------------------------*/
//   ostringstream s;
  
//   string path = "/home/sebld/DOCTORAT/GPS/NOMADr/problems/Hedar_Fukushima/G2/50/observer/";

//   unsigned t   = time(0);
//   unsigned pid = getpid();

//   s << path << "bb_" << t << "_" << pid << ".txt";
//   ofstream out ( s.str().c_str() );

//   out << g1 << " " << g2 << " " << z << endl;
//   out.close();

//   ostringstream s2;
//   s2 << path << "ready_bb_" << t << "_" << pid << ".txt";
//   out.open ( s2.str().c_str() );


//   out.close();

  return 0;
}
