
#include <iostream.h>

void functie( double f( double ) )
{
  cout << f( 0 ) << endl;
}


double dezeFunctie( double d )
{
  return d + 1;
}

int main()
{
  functie( dezeFunctie );
  return 0;
}


