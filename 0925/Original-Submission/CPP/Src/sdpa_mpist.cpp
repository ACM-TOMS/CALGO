/*--------------------------------------------------
  rsdpa_mpist.cpp
--------------------------------------------------*/

#include "sdpa_include.h"
#include "sdpa_scalapack.h"
#include <cmath>


namespace sdpa {

int MpiSt::iam    = 0;
int MpiSt::ictxt  = 0;
int MpiSt::myrow  = 0;
int MpiSt::mycol  = 0;
int MpiSt::nprow  = 0;
int MpiSt::npcol  = 0;
int MpiSt::nprocs = 0;
  // int MpiSt::mb     = 1;
int MpiSt::ictxt2 = 0;
int MpiSt::myrow2 = 0;
int MpiSt::mycol2 = 0;
int MpiSt::nprow2 = 0;
int MpiSt::npcol2 = 0;
#ifdef SCALAPACK_BLOCK
  int MpiSt::mb2    = SCALAPACK_BLOCK;
#else
  int MpiSt::mb2    = 48;
  /*---------------------
    mb2 is a block size of Two-Dimensional Block-Cyclic Distribution
    The defalut value '48' is decided by refering
    MUMPS default block size KEEP(51)
    at Lines 968-976 of dmumps_part2.F
    ---------------------*/
#endif


void MpiSt::initialize(int& argc, char** & argv)
{
  MPI_Init(&argc,&argv);
  Cblacs_pinfo(&iam, &nprocs);
  // rMessage("iam = " << iam << " nprocs = " << nprocs);
  Cblacs_get(-1,0,&ictxt);
  Cblacs_gridinit(&ictxt,(char*)"Column-major",1,nprocs);
  Cblacs_gridinfo(ictxt,&nprow,&npcol,&myrow,&mycol);
  Cblacs_get(-1,0,&ictxt2);
  nprow2 = (int) sqrt((double)nprocs);
  npcol2 = nprocs/nprow2;
  while (nprow2*npcol2 != nprocs) {
    nprow2 --;
    npcol2 = nprocs/nprow2;
  }
  Cblacs_gridinit(&ictxt2,(char*)"Row-major",nprow2,npcol2);
  Cblacs_gridinfo(ictxt2,&nprow2,&npcol2,&myrow2,&mycol2);
}

void MpiSt::finalize()
{
  Cblacs_barrier(ictxt,(char*)"All");
  Cblacs_gridexit(ictxt);
  Cblacs_exit(0);
}

void MpiSt::display(FILE* fpOut)
{
  if (fpOut != NULL) {
    fprintf(fpOut," iam = %d , nprocs = %d , ictxt = %d \n",
	    iam , nprocs,ictxt);
    fprintf(fpOut," myrow  = %d , mycol  = %d , nprow  = %d , npcol  = %d \n",
	    myrow,mycol,nprow,npcol);
    fprintf(fpOut," myrow2 = %d , mycol2 = %d , nprow2 = %d , npcol2 = %d \n",
	    myrow2,mycol2,nprow2,npcol2);
  }
}

void MpiSt::barrier()
{
  Cblacs_barrier(ictxt,(char*)"All");
}

} // end of namespace 'sdpa'
