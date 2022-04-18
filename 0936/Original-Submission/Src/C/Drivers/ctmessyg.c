// A test driver illustrating us different precisions
#include "cmessy.h"


int main(void)
{
  struct smessy_ty e_s;
  struct dmessy_ty e_d;
  struct qmessy_ty e_q;
  float r_s[2];
  double r_d[2];
  long double r_q[2];
  
  allocate_smessy_interface_(2,1);// Must be called before any other interface calls.
  allocate_dmessy_interface_(2,1);
  allocate_qmessy_interface_(2,1);
  get_smessy_defaults_(e_s); //Illustrates good practice before ?messy call
  get_dmessy_defaults_(e_d); //Illustrates good practice before ?messy call
  get_qmessy_defaults_(e_q); //Illustrates good practice before ?messy call
#ifdef MSWin  // This part a total mystery to Krogh
  system("rd /s /q Results");
  system("md Results"); 
// Needed to capture output in a file for comparison using checkctmessy.
//  OPEN_CMESSY_FILES_(ep, 3, "Results/" OUT_);
#endif
  r_s = { 4.0F * atanf(1.0F), expf(1.0F) };
  r_d = { 4.0 * atan(1.0), exp(1.0) };
  r_q = { 4.0L * atanx(1.0L), expx(1.0L) };
  smessy(e_s, "$NSingle: $NPi=$R$N e=$R", r_s, 0);
  dmessy(e_d, "$NSingle: $NPi=$R$N e=$R", r_d, 0);
  qmessy(e_q, "$NSingle: $NPi=$R$N e=$R", r_q, 0);
}
