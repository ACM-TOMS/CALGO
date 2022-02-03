#include "cmessy.h"
void MESSY_(struct CMESSY_TY_* e, char* const s, ...)
{
  static int flag, j;
  static int toobig=100000;
  char * enumnames[8]={"m_rdat ", "m_rmat ", "m_zdat ", "m_zmat ",
                       "m_idat ", "m_imat ", "m_ix   ", "m_ptext"};
  static int msize[20];
  /* Values passed for sizes of arrays */
  static int nd,nda[2],ni,nia[2],nz,nza[2],nx,plen,slen;
  /* Pointers to arrays and a char that may be used within messy.*/
  static CTYPE_* D;
  static CTYPE_* DA;
  static ck* Z;
  static ck* ZA;
  static int* II;
  static int* IA;
  static int* ix;
  static char* ptext;
  static char* ptextt=" "; // Pointed to when ptext is not present.
  /* Holds error data */
  static int itemp[4];
  va_list ap;
#if numt_ > 1
#pragma omp threadprivate(flag, j, msize,nd,nda,ni,nia,nz,nza,nx,plen,slen)
#pragma omp threadprivate(D,DA,Z,ZA,II,IA,ix,ptext,ptextt,itemp)
#endif
  /* Make check that an initialized copy of ?messy_ty is used. */
  if(e->cinit != 1234565) /* This arbitrary value is set as
                             a default within ?messy_ty, in Fortran.
                             If a C program did not make a call to
                             get_?messy_defaults(), then this value
                             will not be there.  Hence the warning.*/
    {
      GET_CMESSY_DEFAULTS_(e); // Intializes and sets e->cinit to 1234565
      strcpy(e -> ename, MESSY_NAME_);
      D=NULL;
      itemp[0]=3;
      //       Recursive call to cmessy that prints out an error message.
      MESSY_(e, "$E18 No initialization call to '" GET_CMESSY_DEFAULTS_NAME_
              "'for the " MESSY_NAME_ "_ty structure.$NAn initialization call"
              "was made within " MESSY_NAME_ ".$NNote: Any settings of struct"
              MESSY_NAME_ "_ty in user code may be overwritten.",
              m_idat,1,itemp,0);
    }
  va_start(ap, s);
  flag=-1;
  /* Give valid values to sizes and pointers.*/
  nd=0;nda[0]=0;nda[1]=0; D=NULL;DA=NULL;
  ni=0;nia[0]=0;nia[1]=0; II=NULL;IA=NULL;
  nz=0,nza[0]=0;nza[1]=0; Z=NULL;ZA=NULL;
  nx=0; ix=NULL; ptext=ptextt;
  for(j=0;j<20;j++) msize[j]=0; /*Record dimensions for finding errors */
  j=0;
  while(flag != 0) /* Process "optional" argument list past the fixed parameters (e, s, ...) */
    {
      flag=va_arg(ap,int); /* Get next enumerated flag for setting up messy */
      switch(flag)
        {	
        case m_rdat:
          nd= va_arg(ap,int); /* Get array size */
          msize[0]=nd;
          D = va_arg(ap, CTYPE_*);/* Get pointer to rank1 assumed size array*/
          j=0;
          break;

        case m_rmat:
          nda[0]=va_arg(ap,int); /* Get array row size */
          nda[1]=va_arg(ap,int); /* Get array col size */
          msize[2]=nda[0];msize[3]=nda[1];
          DA = va_arg(ap, CTYPE_*);   /* Get pointer to rank2 assumed size array*/
          j=1;
          break;

        case m_zdat:
          nz=va_arg(ap,int); /* Get array size */
          msize[6]=nz;
          Z = va_arg(ap, ck*);
          j=3;
          break;

        case m_zmat:
          nza[0]=va_arg(ap,int); /* Get array row size */
          nza[1]=va_arg(ap,int); /* Get array col size */
          msize[8]=nza[0]; msize[9]=nza[1];
          ZA = va_arg(ap, ck*);   /* Get pointer to rank2 assumed size array*/
          j=4;
          break;

        case m_idat:
          ni=va_arg(ap,int); /* Get array size */
          msize[12]=ni;
          II = va_arg(ap, int*); /* Get pointer to rank1 assumed size array*/
          j=6;
          break;

        case m_imat:
          nia[0]=va_arg(ap,int); /* Get array row size */
          nia[1]=va_arg(ap,int); /* Get array col size */
          msize[14]=nia[0]; msize[15]=nia[1];
          IA = va_arg(ap, int*);   /* Get pointer to rank2 assumed size array*/
          j=7;
          break;
        case m_ix:
          nx = va_arg(ap,int); /* Get array size */
          ix = va_arg(ap, int*); /* Get pointer to rank1 assumed size array*/
          msize[18]=nx;
          j=9;
          break;

        case m_ptext:
          ptext =  va_arg(ap, char*); // Pointer to char arguments inserted
          j=10;
          break; 

        case 0:
          break; //Signals that there are no more optional arguments.

        default: // Fatal error condition. No valid stop flag (0) on
                 // the groups.
          itemp[0]=1; itemp[1]=flag;
          strcpy(e -> ename, MESSY_NAME_);
          MESSY_(e, "$E99Undefined option flag with value = $I when calling"
                 MESSY_NAME_ "().$NWas the required trailing 0 or arguments"
                 " neglected?$NThis may be either near the beginning of"
                 " optional groups or it could be near group $P$N",
                 m_idat,2,itemp, m_ptext, enumnames[j], 0);
        } /* end of switch */
    } /* end while */
  va_end(ap); /* End of processing variable argument list */
  /* Check the dimensions of the arrays.  They should be all
     nonnegative and of reasonable size , <= toobig*/
  for(j=0;j<10;j++)
    {
      if( (msize[2*j]      >= 0 && msize[2*j]   <= toobig)
          && (msize[2*j+1] >= 0 && msize[2*j+1] <= toobig)) continue;

      itemp[0]=2; itemp[1]=msize[2*j];itemp[2]=msize[2*j+1]; itemp[3]=j+1;
      MESSY_(e,"$E99Dimensions of array are too big: [$I, $I].$N"
             "Probably forgotten or misplaced arguments using group $I,"
             " with name $P$N",m_idat,4,itemp, m_ptext,
             strlen(enumnames[j]), enumnames[j], 0);
    }
  /* Call Fortran wrapper routine that calls messy().
     Every optional argument is provided here but only some are usually
     activated with non-zero sizes.*/
  slen = strlen(s);
  plen = strlen(ptext);

  CALLMESSY_(e, s, slen,
             nd, D,
             nda[0], nda[1], DA,
             nz, Z,
             nza[0], nza[1], ZA,
             ni, II,
             nia[0], nia[1], IA,
             nx, ix,
             ptext,plen);
}
