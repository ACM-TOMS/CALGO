#ifndef VBF_vec_GF2E__H
#define VBF_vec_GF2E__H

#include <NTL/vec_GF2E.h>
#include <NTL/fileio.h>

NTL_CLIENT

class vbf_vec_GF2E: public vec_GF2E {

  void conv(vec_GF2E& x, const mat_GF2& A);

  void conv(mat_GF2& A, const vec_GF2E& x);

}; // end class vbf_vec_GF2E

void conv(vec_GF2E& x, const mat_GF2& A)
{
   ofstream output;
   OpenWrite(output, "mat_GF2_to_vec_GF2E.tmp");
   output << A;
   output.close();

   ifstream input;
   OpenRead(input, "mat_GF2_to_vec_GF2E.tmp");
   input >> x;
   input.close();
}

void conv(mat_GF2& A, const vec_GF2E& x)
{
   ofstream output;
   OpenWrite(output, "mat_GF2_to_vec_GF2E.tmp");
   output << x;
   output.close();

   ifstream input;
   OpenRead(input, "mat_GF2_to_vec_GF2E.tmp");
   input >> A;
   input.close();
}
#endif
