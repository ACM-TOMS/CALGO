/* -----------------------------------------------------------
*  sdpa_rdpotf2.f
*  modification of pdpotf2.f (ScaLAPACK 1.7)
*  to call the function 'rdpotrfL'.
*
*  modified by Makoto Yamshita 2002.07.11
* 
* -----------------------------------------------------------
*  rsdpa_rdpotrf.f
*  modification of pdpotrf.f
*  to call the function 'rpdpotf2'.
*
*  modified by Makoto Yamshita 2002.07.11 
*
* -----------------------------------------------------------
*
*     Further modify to remove gfortran warning
*     by Makoto Yamashita 2008.08.30
*
*     Convert to c++
*     by Makoto Yamashita 2008.08.31
*
* -----------------------------------------------------------
*/

#include "sdpa_tool.h"
#include "sdpa_algebra.h"
#include "sdpa_scalapack.h"
#include "sdpa_dpotrf.h"
#include "sdpa_mpist.h"
#include <cmath> //for ceil
using namespace std;
using namespace sdpa;

#define DLEN1_            9
#define DTYPE1_           0
#define CTXT1_            1
#define M1_               2
#define N1_               3
#define MB1_              4
#define NB1_              5
#define RSRC1_            6
#define CSRC1_            7
#define LLD1_             8

extern "C" {
  
int rpdpotf2_(char* uplo, int* n, double* A,
	      int* ia, int* ja, int desca[DLEN1_], int* info)
{
  char colbtop, rowbtop;
  int  iacol, iarow, ictxt, idiag, iia, jja, lda;
  int  mycol, myrow, npcol, nprow;

  ictxt = desca[CTXT1_];
  Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
  #if 0
  rMessage("ictxt = " << ictxt
	   << " nprow = " << nprow
	   << " npcol = " << npcol
	   << " myrow = " << myrow
	   << " mycol = " << mycol);
  MpiSt::barrier();
  #endif
  *info = 0;
  if (*n == 0) {
    return 0;
  }
  // Do not use PB_Cinfog2l
  // because PB_Cinfog2l and infog2l use different desca
  infog2l_(ia,ja, desca, &nprow, &npcol, &myrow, &mycol,
	      &iia, &jja, &iarow, &iacol);
  #if 0
  rMessage("ia = " << *ia << " ja = " << *ja
	   << " iia = " << iia << " jja = " << jja
	   << " iarow = " << iarow << " iacol = " << iacol);
  #endif
  PB_topget_(&ictxt,(char*)"BroadCast",(char*)"Rowwise", &rowbtop);
  PB_topget_(&ictxt,(char*)"BroadCast",(char*)"Columnwise", &colbtop);
  if (mycol == iacol) {
    if (myrow == iarow) {
      lda = desca[LLD1_];
      idiag = (iia-1) + (jja-1)*lda;
      #if 0
      rMessage("idiag = " << idiag
	       << " iia = " << iia << " jja = " << jja);
      #endif
      rdpotrfl_(n, &A[idiag], &lda, info);
      Cigebs2d(ictxt,(char*)"Columnwise",&colbtop,1,1,(char*)info,1);
    }
    else {
      Cigebr2d(ictxt,(char*)"Columnwise",&colbtop, 1, 1, (char*)info,
	       1, iarow, mycol);
    }
    Cigebs2d(ictxt,(char*)"Rowwise", &rowbtop, 1, 1, (char*)info, 1 );
  }
  else {
    Cigebr2d(ictxt,(char*)"Rowwise", &rowbtop, 1, 1, (char*)info, 1,
	     myrow, iacol );
  }
  return 0;
}

void recoverTopology(int ictxt, char rowbtop, char colbtop)
{
  PB_topset_(&ictxt,(char*) "Broadcast", (char*) "Rowwise",
	     &rowbtop);
  PB_topset_(&ictxt, (char*) "Broadcast", (char*)"Columnwise",
	     &colbtop);
}


int rpdpotrf_(char* uplo, int* n, double* A,
	      int* ia, int* ja,int desca[DLEN1_],int* info)
{
  // This routine is only for lower triangular
  if (uplo[0] != 'L' && uplo[1] != 'l') {
    *info = -1;
    return -1;
  }

  char colbtop, rowbtop;
  int  i, icoff, ictxt, iroff, j, jb, jn, mycol, myrow, npcol, nprow;
  int  idum1[1], idum2[1];

  int one = 1;
  int two = 2;
  int six = 6;
  double done =   1.0;
  double dmone = -1.0;

  ictxt = desca[CTXT1_];
  Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
  *info = 0;
  if (nprow == -1) {
    *info = -(600+CTXT1_+1);
  }
  else {
    chk1mat_(n,&two,n,&two,ia,ja,desca,&six,info);
    if (*info == 0) {
      iroff = (*ia-1) % desca[MB1_];
      icoff = (*ja-1) % desca[NB1_];
      if (iroff != 0) {
	*info = -4;
      }
      else if (icoff != 0) {
	*info = -5;
      }
      else if (desca[MB1_] != desca[NB1_]) {
	*info = -(600 + NB1_ + 1);
      }
    }

    idum1[0] = 'L';
    idum2[0] = 1;
    pchk1mat_(n, &two, n, &two, ia, ja, desca, &six, &one,
	      idum1, idum2, info);
  }

  if (*info != 0) {
    *info = - (*info);
    pxerbla_(&ictxt, (char*) "pdpotrf", info, 7);
    return 0;
  }
  // Quick return if possible
  if (*n == 0) {
    return 0;
  }
  
  PB_topget_(&ictxt, (char*)"Broadcast", (char*)"Rowwise",
	     &rowbtop);
  PB_topget_(&ictxt, (char*)"Broadcast", (char*)"Columnwise",
	     &colbtop);
  
  // A is lower triangular, compute Cholesky factorization A = L*L'
  PB_topset_(&ictxt, (char*)"Broadcast", (char*)"Rowwise",
	     (char*)"S-ring");
  PB_topset_(&ictxt, (char*)"Broadcast", (char*)"Columnwise",
	     (char*)" ");
  //(right-looking)
  //      Handle the first block of columns separately

  jn = iceil_(ja,&desca[NB1_]) * desca[NB1_];
  if (*ja + *n -1 < jn) {
    jn = *ja + *n -1;
  }
  jb = jn - *ja + 1;

  // Perform unblocked Cholesky factorization on JB block
  rpdpotf2_((char*)"L", &jb, A, ia, ja, desca, info);
  if (*info != 0) {
    recoverTopology(ictxt,rowbtop,colbtop);
    return -1;
  }

  if (jb + 1 <= *n) {

    int n_m_jb  = *n - jb;
    int ia_p_jb = *ia + jb;
    int ja_p_jb = *ja + jb;
    
    // Form the column panel of L using the triangular solver
    pdtrsm_((char*)"Right",(char*)"Lower",
	    (char*)"Tranpose",(char*)"Non-Unit",
	    &n_m_jb, &jb, &done, A, ia, ja, desca,
	    A, &ia_p_jb, ja, desca);
    pdsyrk_((char*)"Lower", (char*)"Notranspose",
	    &n_m_jb, &jb, &dmone,
	    A, &ia_p_jb, ja, desca, &done,
	    A, &ia_p_jb, &ja_p_jb, desca);
    
  }

  for (j = jn+1; j<= *ja + *n -1; j+=desca[NB1_]) {
    // Computing MIN
    jb = *n - j + *ja;
    if (desca[NB1_] < jb) {
      jb = desca[NB1_];
    }
    i = *ia + j - *ja;
    // Perform unblocked Cholesky factorization on JB block
    rpdpotf2_((char*)"Lower", &jb, A, &i, &j, desca, info);
    
    if (*info != 0) {
      *info = *info + j - *ja;
      recoverTopology(ictxt,rowbtop,colbtop);
      return -1;
    }

    if (j - *ja + jb + 1 <= *n) {
      // Form the column panel of L using the triangular solver
      int n_m_j_m_jb_p_ja = *n - j - jb + *ja;
      int i_p_jb = i+jb;
      int j_p_jb = j+jb;
      pdtrsm_((char*)"Right",(char*)"Lower",
	      (char*)"Transpose",(char*)"Non-Unit",
	      &n_m_j_m_jb_p_ja, &jb, &done, A, &i, &j, desca,
	      A, &i_p_jb, &j, desca);
      // Update the trailing matrix, A = A - L*L'
      pdsyrk_((char*)"Lower", (char*)"No Transpose",
	      &n_m_j_m_jb_p_ja, &jb, &dmone,
	      A, &i_p_jb, &j, desca, &done,
	      A, &i_p_jb, &j_p_jb, desca);
    }
  }
  
  recoverTopology(ictxt,rowbtop,colbtop);
  return 0;
}

}; // end of extern "C"
