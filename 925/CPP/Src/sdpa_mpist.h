/*--------------------------------------------------
  sdpa_mpist.h
--------------------------------------------------*/

#ifndef __sdpa_mpist_h__
#define __sdpa_mpist_h__

#include <cstdio>

#define MPICH_IGNORE_CXX_SEEK
#include <mpi.h>

namespace sdpa {

class MpiSt
{
public:
  static int iam;
  static int nprocs;
  static int ictxt;
  static int myrow;
  static int mycol;
  static int nprow;
  static int npcol;
  // static int mb;
  static int ictxt2;
  static int myrow2;
  static int mycol2;
  static int nprow2;
  static int npcol2;
  static int mb2;
  static void initialize(int& argc, char**& argv);
  static void finalize();
  static void display(FILE* fpOut=stdout);
  static void barrier();
};

// rMpiSt::mb2 defines block size of
// two cyclic distribution of ScaLAPACK

}; // end of namespace 'sdpa'

#endif // __sdpa_mpist_h__
