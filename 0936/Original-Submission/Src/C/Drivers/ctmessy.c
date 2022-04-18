// 
// Code written by Fred T. Krogh with Richard J. Hanson adding MS Windows support.

/* Note that this was written to mirror the Fortran version as much as possible.
   This tends to put the C code in a bad light.  Things that were easy to state
   in Fortran code are awkward in the C code.  If a driver designed to be easy
   to write in C, were to be converted to Fortran, one may see a slightly
   similar effect, but probably not to the same degree.
*/
#undef FFMT
#if plet_=='d'
#define CTYPE_ double
#define MESSY_ dmessy
#define CALLMESSY_ calldmessy
#define SAMPLE_ dsample
#define ALLOCATE_CMESSY_INTERFACE_ allocate_dmessy_interface
#define DEALLOCATE_CMESSY_INTERFACE_ deallocate_dmessy_interface
#define GET_CMESSY_DEFAULTS_ get_dmessy_defaults
#define OPEN_CMESSY_FILES_ open_dmessy_files
#define CMESSY_TY_ dmessy_ty
#define SAMPLE_TY_ dsample_ty
#define DIG DBL_DIG
#define FFMT "%7.2f"
#define OUT_ "result.d"
#elif plet_=='s'
# define CTYPE_ float
#define MESSY_ smessy
#define CALLMESSY_ callsmessy
#define SAMPLE_ ssample
#define ALLOCATE_CMESSY_INTERFACE_ allocate_smessy_interface
#define DEALLOCATE_CMESSY_INTERFACE_ deallocate_smessy_interface
#define GET_CMESSY_DEFAULTS_ get_smessy_defaults
#define OPEN_CMESSY_FILES_ open_smessy_files
#define CMESSY_TY_ smessy_ty
#define SAMPLE_TY_ ssample_ty
#define  DIG FLT_DIG
#define FFMT "%7.2f"
#define OUT_ "result.s"
#elif plet_=='q'
#define CTYPE_ long double
#define MESSY_ qmessy
#define CALLMESSY_ callqmessy
#define SAMPLE_ qsample
#define ALLOCATE_CMESSY_INTERFACE_ allocate_qmessy_interface
#define DEALLOCATE_CMESSY_INTERFACE_ deallocate_qmessy_interface
#define GET_CMESSY_DEFAULTS_ get_qmessy_defaults
#define OPEN_CMESSY_FILES_ open_qmessy_files
#define CMESSY_TY_ qmessy_ty
#define SAMPLE_TY_ qsample_ty
#define DIG LDBL_DIG
#define FFMT "%7.2Lf"
#define OUT_ "result.q"
#endif
#define rk CTYPE_
#include "csample.h"
#include <complex.h>

int main(void)
{
  struct CMESSY_TY_ e;
  struct CMESSY_TY_ * const ep =&e;
  struct SAMPLE_TY_ s;
  struct SAMPLE_TY_ * const sp=&s;
  int idat[6] = { 2, 1, 1, 1, 6, 6};
  int i, ix[1], j, k, m;
  CTYPE_ c0;
  CTYPE_ rdat[5] = {(rk)0.,(rk)2.328306436538696E-10, (rk)4.440892097466651E+4,
		(rk)1.e-20, (rk)1.421085471189328E-14};
  int idat1[6] = {0, -0, 17, -17, 8922, -8922 };
  rk rdat1[6] = { (rk)0.0, (rk)0.0, (rk)17.123456789, (rk)-.17123456789,
                  (rk)4.123456789e6, (rk)3.123456789e-4 };
  rk rdat2[5] = { (rk)9.994, (rk)9.9996, (rk).99994999, (rk).999999999999999,
                  (rk)9.99949999 };

  int idat2[16] = {0, -0, 17, -17, 8922, -8922, 7, 123, 5, -12, 16, -26,
                   13, 267, 48, 9 };
  int imat[15][5];
  rk rmat[15][5];
  rk rdat3[15];
  int idat3[9];
  int idat4[6] = {13, 0, 15, 6, 7, 7};
  rk rdat4[3] = {(rk)6.636e-2, (rk)2.789e-2, (rk)0.0};
  int idat5[4] = {13, 2, 1, 8};
  rk rdat5[6] = {(rk)6.63568416E-02, (rk)2.78940E-02, (rk)1.59E-01, (rk)2.71E-02,
                 (rk)7.91E+01, (rk)1.42500};
  rk rdat6[11] = {(rk)7.1E-03, (rk)1.0E+00, (rk)1.0E-09, (rk)9.9332068E-01,
                  (rk)2.392E-06, (rk)-2.223E-08, (rk)-8.231E-09, (rk)-3.193E-09,
                  (rk)1.1E-08, (rk)0.00E+00, (rk)4.36246E+00};
  rk rdat7[11] = {(rk)2.8E-04, (rk)8.1E-03, (rk)1.0E-09, (rk)-9.3828367E-02,
                  (rk)-1.307E-07, (rk)1.103E-09, (rk)2.654E-10, (rk)2.654E-10,
                  (rk)4.4E-09, (rk)0.00E+00, (rk)6.09012E+00};
  rk rdat8[11] = {(rk)2.4E-04, (rk)6.9E-03, (rk)1.0E-09, (rk)-4.7748202E-02,
                  (rk)-9.342E-08, (rk)8.125E-10, (rk)2.265E-10, (rk)2.265E-10,
                  (rk)7.3E-09, (rk)0.00E+00, (rk)6.09012E+00};
  rk rdat9[6] = { (rk)9.42508693E-02, (rk)3.13808E-02, (rk)6.57E-02, (rk)9.48E-02,
                  (rk)3.29E+01, (rk)1.42500E+00};
  rk rdat10[11] = {(rk)5.1E-03, (rk)7.3E-02, (rk)1.0E-09, (rk)9.8817022E-01,
                   (rk)-8.535E-08, (rk)-2.466E-08, (rk)4.624E-09, (rk)4.624E-09,
                   (rk)9.3E-09, (rk)0.00E+00, (rk)3.55720E+00};
  rk rdat11[11] = {(rk)2.6E-03, (rk)3.7E-02, (rk)1.0E-09, (rk)-1.2463414E-01,
                   (rk)-3.498E-07, (rk)6.896E-09, (rk)2.973E-09, (rk)1.733E-09,
                   (rk)3.0E-08, (rk)0.00E+00, (rk)3.55720E+00};
  rk rdat12[11] = {(rk)1.5E-03, (rk)2.1E-02, (rk)1.0E-09, (rk)-6.3317119E-02,
                   (rk)-2.503E-07, (rk)4.775E-09, (rk)1.885E-09, (rk)8.272E-10,
                   (rk)3.4E-08, (rk)0.00E+00, (rk)3.55720E+00};
  rk rdat13[5] = {(rk).123, (rk).124, (rk).125, (rk).126, (rk).127};
  ck zdat[50]; ck zmat[10][5];
  char* const diva1=
    "$N$D4T=$R  H=$R, NS=$I NSE=$I$D3 E=$R ND=$I KQ=$W$B";
  char* const diva2="KSTEP=$I$D9 T=$R $D6H=$R LSC=$I EIMIN=$D3$R"
    " EAVE=$R KSC=$I SIGMA($I) RQ=$D6$R$B";
  char* const diva3="$H|5RI|3RKQ|3RLI|8CE|8CEI|8CEPS|CF|"
    "44CHigh Order Predicted Differences|9C Rnoise|10CStiff|13CBeta|";
  char* const diva4="$FI5$2FI3$3FE8.1$FD$4FE11.3$FE9.1$FE10.2$FE13.5";
  char* const diva5="$H|5RI|3RKQ|3RLI|8CE|8CEI|8CEPS|CF$Y|"
    "11CD$I$y|9C Rnoise|10CStiff|13CBeta|";
  char* const diva6="$FI5$2FI3$3FE8.1$FD$Y$FE11.3$y$FE9.1$FE10.2$FE13.5";
  const int bitvec[8] ={ // Note only 31 bits per integer
    0x2AAA,0x55555555,0x1999,0x4CCCCCCC,0x2468,0x56789ABC,0x1DDD,0x6EEEEEEE};
  char tmp_out[12]; // For showing how C output looks.

#ifdef NAG
  f90_io_init(void);
#endif
  /* *********************************************************************
*********************** Start of examples *****************************
**************************************ddd*********************************/
  // void allocate_?messy_interface(int maxcstructs, int maxcthreads);
  //     The value of maxcstructs must be no smaller than the number
  //     of different C structures ?messy_ty  used in a C program.

  //     The value of maxcthreads must be no smaller than the
  //     number of parallel OpenMP threads used in a C program.

  ALLOCATE_CMESSY_INTERFACE_(2,1);// Must be called before any other interface calls.

  GET_CMESSY_DEFAULTS_(ep); //Illustrates good practice before cmessy call

#ifdef MSWin
  system("rd /s /q Results");
  system("md Results"); 
// Needed to capture output in a file for comparison using checkctmessy.
  OPEN_CMESSY_FILES_(ep, 3, "Results/" OUT_);
#endif
  MESSY_(ep, "$L75If you don't make a call to get_$Pmessy_defaults, before the"
         " first call to $Pmessy you will get an error message telling you that"
         " a call to get_$Pmessy_defaults was made within cmessy, and settings"
         " in struct $Pmessy_ty may be overwritten.",m_ptext,qplet_,0);

  strcpy(ep->ename,"ctmessy"); // Setting the package name for error messages.
  ep->fpprec = e.fpprec>32?32:e.fpprec; // So NAG compiler will match others

  ep->line_len = 75;
  s.se = e;
  MESSY_(ep, "$NDefault line lengths have been changed from 128 to 75."
         "This change makes it easier to show some of the features of $Pmessy"
         " that would not be illustrated as nicely with the longer line lengths." 
         " Below is an example of a non-stopping error in the integration program"
         " diva.", m_ptext,qplet_, 0);
  MESSY_(ep,"$E25At: TN=$R, KSTEP=$I, with H=$R, Error tolerances are too"
         " small.  (Estimated Error) / (Requested Error) for equation $I is $R."
         "  Tolerance $I for F($I) = $R.  Replacing F($I) with $R.",
         m_idat,6,idat, m_rdat,5,rdat, 0);
  c0 = (rk)(idat1[2] + idat1[3]); // Used to show how a -0 will print
  rdat1[1] = -c0; MESSY_(ep,"$NYou might get the next line, using an i0 format, but there is"
         " no similar format for the following line of reals with 8 sig. digits."
         "  The next output is for the same data when 4 digits after the decimal"
         " are requested.$N$Nidat(1:6) is: i1=$I i2=$I i3=$I i4=$I i5=$I i6=$I$N"
         "$D8rdat(1:6) is: r1=$R r2=$R r3=$R r4=$R r5=$R r6=$R$B",
         m_idat,6,idat1, m_rdat,6,rdat1, 0);
  MESSY_(ep,"$D-4rdat(1:6) is: r1=$R r2=$R r3=$R r4=$R r5=$R r6=$R$B",
         m_rdat,6,rdat1, 0);

  MESSY_(ep, "rdat(1:6) is: r1=$FF8.4 r2=$G r3=$G, r4=$G r5=$G r6=$G",
         m_rdat,6,rdat1, 0);

  sprintf(tmp_out, "%3d " FFMT, idat1[1], rdat1[1]); // FFMT needed for quad p.
  /* This looks weird, but to guarantee the C output goes in the proper place
     we need the output to go through the Fortran I/O units. */
  MESSY_(ep, "$NTo show how C handles cases with (a possibly) negative 0,"
         " here are$Nresults for the same idat(2) and rdat(2) with a C"
         " printf statement$N$P$Nmessy avoids printing the \"-\", C may not.",
         m_ptext,tmp_out, 0);
  MESSY_(ep,"$N$D8rdat as a vector with a starting index of -3:$V$O-3$B",
         m_idat,1,idat1, m_rdat,6,rdat1, 0);

  MESSY_(ep,"$NNext we check the close calls on the rounding.  In order to"
         " protect against a close call giving a bad format we sometimes get more"
         " digits and/or a wider field that would be absolutely necessary."
         "  We are asking for 4 significant figures.  Imagine what you would want"
         " to see for the numbers: 9.994, 9.9996, .99994999, .999999999999999,"
         " 9.99949999$N$N$D4rdat(1:5) is:  r1=$R r2=$R r3=$R r4=$R r5=$R",
         m_rdat,5,rdat2, 0);

  MESSY_(ep, "$N$D4An idat as a vector:$WThen last rdat:$V$B",
         m_idat,16,idat2, m_rdat,5,rdat2,0);
  MESSY_(ep, "Same with no offset:$W",m_idat, 16, idat2, 0);
  MESSY_(ep, "Same with offset:$W$O10",m_idat, 16, idat2, 0);
  
  
  MESSY_(ep,"$NThe last vectors could be printed like this as well:$N"
         "$D4An idat as a vector:$N$WThen last rdat:$N$V$B",
         m_idat,16,idat2, m_rdat,5,rdat2, 0);

  k = 1;
  for (j=0; j<15; j++) {
    for (i=0; i<5; i++) { // Want numbers that match on different machines.
      rmat[j][i] = (rk)(-12.9375) + ((rk)k)/(rk)(1<<(k%4));
      k++;
      if (rmat[j][i] >= (rk)0.0) {
        m = (int)(rmat[j][i]+(rk).5);
      }
      else {
        m = (int)(rmat[j][i]-(rk).5);
      }
      imat[j][i] = m;
    }
  }
  MESSY_(ep, "$Nimat:$M", m_imat, 5,15,imat, 0);
  MESSY_(ep,"$NThis would be more compact with shorter column labels."
         "  Thus:$N$Nimat:$M$O<1C$O$O",  m_imat,5,15,imat, 0);
  MESSY_(ep,"$NOriginal, but for the transposed matrix:$M$T", m_imat,5,15,imat, 0);
  MESSY_(ep,"$NHow about just some of this data with fancier labels?$N"
         "$Nimat:$M$O 5Earth  Air FireWater$O 8HydrogenCopper  Iron    Gold    "
         "Lead    $O",   m_imat,5,4,imat,0);

  MESSY_(ep,"$N$D6Transpose real matrix printing (6 significant digits)."
         "$Nrmat^T:$A$T$B", m_rmat, 5, 15, rmat, 0);

  MESSY_(ep,"$N$D6Just default real matrix printing (6 significant digits)."
         "$Nrmat:$A$B", m_rmat, 5, 15, rmat, 0);
  // Note rmat has 5 rows and 15 columns.  (Row index is last index.)
  MESSY_(ep,"$N$D6Same, but testing strange labels and indexes.$Nrmat:"
         "$A$O-2|11<>$O-3<2R:$O$B", m_rmat,5,15,rmat, 0);

  ix[0] = 2;
  for (j=0; j<15; j++) {
    rdat3[j] = rmat[j][1];
  }
  MESSY_(ep,"$N$D6Row $X of rmat:$V$B", m_rdat,15,rdat3, m_ix,1,ix, 0);



  idat3[0]=e.fpprec; idat3[1]=e.line_len; idat3[2]=e.maxerr; idat3[3]=e.lstop;
  idat3[4]=e.lprint; idat3[5]=e.errcnt; idat3[6]=e.dblev;
  MESSY_(ep,"$NUsing tabs to display portable public integers in messy_ty:"
         "$14T$Nfpprec=$I$Tline_len=$I$Tmaxerr=$I$Tlstop=$I$Tlprint=$I$T"
         "errcnt=$I$Tdblev=$I", m_idat,7,idat3, 0);

  MESSY_(ep,"$NData entered by hand to show what output from the 'DIVA'"
         " integration program would look like.  Even though this output was"
         " designed for 128 character lines, it can still be understandable with"
         " shorter lines.", 0);

  // e.line_len=128 ! Uncomment this to see what output should look like.
  MESSY_(ep,diva1, m_idat,6,idat4, m_rdat,3,rdat4, 0);

  MESSY_(ep,diva2, m_idat,4,idat5, m_rdat,6,rdat5,0);
  e.kdf=(DIG >= 9) ? 9:DIG;
  MESSY_(ep,diva3,0);
  idat5[0]=1; idat5[1]=6; idat5[2]=1;

  MESSY_(ep, diva4, m_idat,3,idat5, m_rdat,11,rdat6,0);
  idat5[0]=2; idat5[1]=7; idat5[2]=1;
  MESSY_(ep, diva4, m_idat,3,idat5, m_rdat,11,rdat7, 0);
  idat5[0]=3; idat5[1]=7; idat5[2]=1;
  MESSY_(ep, diva4, m_idat,3,idat5, m_rdat,11,rdat8, 0);
  MESSY_(ep, "$B", 0); 
  rdat4[0]=(rk)9.425E-02; rdat4[1]=(rk)3.138E-02;
  idat4[0]=14; idat4[2]=16; idat4[3]=7;
  MESSY_(ep,diva1,m_idat,6,idat4, m_rdat,3,rdat4, 0);

  idat5[0]=14; idat5[1]=2; idat5[2]=1; idat5[3]=8;
  MESSY_(ep,diva2, m_idat,4,idat5, m_rdat,6,rdat9, 0);
  e.kdf=(DIG >= 9) ? 9:DIG;
  MESSY_(ep, "The following is just to illustrate $$Y...$$y actions.", 0);
  idat5[0]=6; idat5[1]=7; idat5[2]=8; idat5[3]=9;
  idat4[0]=4;
  MESSY_(ep,diva5, m_ix,1,idat4, m_idat,4,idat5, 0);
  idat5[0]=1; idat5[1]=7; idat5[2]=0;
  MESSY_(ep, diva6, m_idat,3,idat5, m_rdat,11,rdat10, m_ix,1,idat4, 0);
  idat5[0]=2; idat5[1]=7; idat5[2]=0;
  MESSY_(ep, diva6, m_idat,3,idat5, m_rdat,11,rdat11, m_ix,1,idat4, 0);
  idat5[0]=3; idat5[1]=7; idat5[2]=0;
  MESSY_(ep, diva6, m_idat,3,idat5, m_rdat,11,rdat12, m_ix,1,idat4, 0);
  MESSY_(ep, "$B",0);

  MESSY_(ep, "$NTesting subtle issues with 'F' format.$Ntest1: $D4$V",
          m_rdat,5,rdat13, 0);
  rdat13[0] = (rk).0942408;
  MESSY_(ep, "test2: t=$D4$R", m_rdat,1,rdat13, 0);
  idat5[0]=3; idat5[1]=9; idat5[2]=23; idat5[3]=40;
  rdat3[0]=(rk)7.46; rdat3[1]=(rk).946; rdat3[2]=(rk)1.78; rdat3[3]=(rk)-9.40;
  idat4[0]=7;
  MESSY_(ep, "Test printing column 7 of a sparse matrix.$NCol $X:$S",
         m_idat,4,idat5, m_rdat,4,rdat3, m_ix,1,idat4, 0);

  MESSY_(ep, "$L50Best to put $$N before what you want at start of a line,"
         " Strange things happen with long preamble+short line. Col $X:$S$B",
         m_idat,4,idat5, m_rdat,4,rdat3, m_ix,1,idat4, 0);

  for (i=0; i<20; i++) zdat[i] = (rk)1.0-(rk).0625*(rk)i + I*((rk).0625*(rk)i);
  MESSY_(ep, "$N$D5Complex data: zdat(1)=$ZR, zdat(2)=$ZR", m_zdat,2,zdat, 0);
  MESSY_(ep, "$N$D5Complex vector zdat:$N$ZV", m_zdat,20,zdat, 0);

  k = 0;
  for (j=0; j<10; j++) {
    for (i=0; i<5; i++) { // Want numbers that match on different machines.
      zmat[j][i]=(rk)1.0-(rk).0625*(rk)k + I*((rk).0625*(rk)k);
      k++;
      if (rmat[j][i] >= (rk)0.0) {
        m = (int)(rmat[j][i]+(rk).5);
      }
      else {
        m = (int)(rmat[j][i]-(rk).5);
      }
      imat[j][i] = m;
    }
  }
  MESSY_(ep,"$N$D5Complex matrix zmat:$N$ZA", m_zmat,5,10,zmat, 0);

  zdat[0] = (rk)1.e12 -(rk)7.3*I; zdat[1]=zdat[0]; zdat[2]=zdat[0];
  zdat[3] = (rk)7.36 + (rk)0.0*I;
  MESSY_(ep, "$NA user formatted complex number with same format for real"
         " and imaginary parts, then one with both formats specified, then two"
         " numbers printed with a default of 4 significant digits.  z1=$ZFE10.3"
         " z2=$ZFE9.3,E10.3 $D4z3=$ZR z4=$ZR", m_zdat,4,zdat, 0);

  MESSY_(ep,"$NJust printing part of text.$K2  And some more.$K4"
         "  But not all of it!", 0);

  MESSY_(ep,"$NA repeated message, with name=$P so we can change names.",
         m_ptext,"\"this name\"", 0);

  idat4[0]=45;
  MESSY_(ep,"$NAn example of binary output: nodec=$QBI$H noinc=$QBI"
         "$H nodecx=$QBI$H noincx=$QBI$H", m_idat,8,bitvec, m_ix,1,idat4, 0);
  idat4[0]=-45;
  MESSY_(ep, "$NThe same no leading 0's and hex output:  nodec=$QZI$H"
         " noinc=$QZI$H nodecx=$QZI$H noincx=$QZI$H",
         m_idat,8,bitvec, m_ix,1,idat4, 0);

  MESSY_(ep, "$NThe same no leading 0's and octal output:  nodec=$QOI$H"
         " noinc=$QOI$H nodecx=$QOI$H noincx=$QOI$H",
         m_idat,8,bitvec, m_ix,1,idat4, 0);
  idat4[0]=45;
  MESSY_(ep, "$NPrinting a bit string vector$Nbitvec:$QBV",
         m_idat,8,bitvec, m_ix,1,idat4, 0);
  MESSY_(ep, "$L50Lining up output with vectors longer than the line length"
         " is tricky.$Nbitvec:$QBV",m_idat,8,bitvec, m_ix,1,idat4, 0);
  MESSY_(ep, "$L40Etc. bitvec:$QBV$B",m_idat,8,bitvec, m_ix,1,idat4, 0);
  MESSY_(ep, "Etc., without ix (number of bits <= 31) and longer line.$N"
         "bitvec:$QBV", m_idat,8,bitvec, 0);

  sp->what = setup_sample; // Code to show various interactions.
  SAMPLE_(sp);
  sp->what = partial_message;

  SAMPLE_(sp); 
  sp->what = finish_message;
  SAMPLE_(sp);
  idat4[0]=sp->se.errcnt; idat4[1]=sp->se.maxerr;
  MESSY_(ep, "Finished first error message for $Psample,  count=$I maxerr=$I",
         m_idat,2,idat4,m_ptext,qplet_, 0);
  sp->what = 17;
  sp->se.lstop = 5; 
  SAMPLE_(sp);

  idat4[0] = 999;

  MESSY_(ep,"$E88This error message illustates how an error message ends things,"
         "ctmessy is done.",  m_idat,1,idat4, 0);

  printf("Done processing\n");
#ifdef NAG
  f90 io finish(void); 
#endif
  
}
