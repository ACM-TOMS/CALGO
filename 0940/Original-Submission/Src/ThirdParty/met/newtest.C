#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char * argv[]) {
  if (argc <= 0) {
    cerr << "Usage : " << argv[0] << " size [size2] ... " << endl;
    exit(1);
  }

  for (int k=1; k<argc; k++) {
    int size = atoi(argv[k]);
    
    double *v1 = new double[size], 
    *v2 = new double[size],
    *v3 = new double[size],
    *v4 = new double[size];

    cout << "v2 -v1 :" << v2 - v1 << endl;
    cout << "v3 -v2 :" << v3 - v2 << endl;
    cout << "v4 -v3 :" << v4 - v3 << endl;
    
    assert(v2 - v1 > size );
    assert(v3 - v2 > size );
    assert(v4 - v3 > size );

    for (int i=0; i<size; i++) {
      v1[i] = 0.1;
      v2[i] = 0.2;
      v3[i] = 0.4;
    }
    for (int i=0; i<size; i++) v4[i] = v1[i]*3.0 + 2.0*v2[i] + v3[i];
  }
  return 0;
}
