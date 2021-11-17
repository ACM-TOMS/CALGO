#include "mex.h"
#include "dsdp5.h"
#include <math.h>

/* #define K_IN prhs[0] */
#define CA_IN prhs[1]
#define B_IN prhs[0]
#define PARS_IN prhs[2]
#define Y_IN prhs[3]
#define STAT_OUT plhs[0]
#define Y_OUT plhs[1]
#define X_OUT plhs[2]
static int DSDPPrintStats2(DSDP, void*);
static int CountNonzeroMatrices(double*,int,int,int,int*,int*,int*);
static int CheckForConstantMat(double*,int,int);

static int printlevel=10;
#define CHKERR(a)  { if (a){mexWarnMsgTxt("DSDP Numerical Error"); } }

#define mexErrMsgTxt(a); mexPrintf("Error: "); mexWarnMsgTxt(a); return;

/*! \file dsdp.c
  \brief Call DSDP from the Matlab environment
 */
/*!
\fn void mexFunction(int nlhs,mxArray *plhs[],int nrhs, const mxArray *prhs[]){
\brief Call DSDP from the Matlab environment.
\param nlhs is the number of output arguments
\param plhs are the output arguments
\param nrhs is the number of input arguments
\param prhs are the input arguments
\ingroup Examples
\note Must be called from Matlab
*/
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
  
  mxArray  *CA_cell_pr,*X_cell_pr;
  const mxArray  *OPTIONS_FIELD;
  mxArray  *STAT_FIELD;
  int     i,j,k,ii,itmp,index,info;
  int    *air, *ajc, *air2, *ajc2, *str,*iptr;
  int   nvars,nb,mC,nC,m1,n1,m2,n2;
  int   nsubs=2, subs[2];
  int   iscellCA;
  int   its,reuse=4,print_info=1,printsummary=0;
  int   ijnnz,spot,n,nn=0,nzmats,vecn;
  int it1,it2;
  int nsdpblocks=0,sdpblockj=0,sdpnmax=1,lpnmax=1,stat1=1,xmaker=0;
  int sspot,nsubblocks,blockj;
  int jj,tnnz,tnnz2;
  int maxit=1000,fastblas=1,rpos=0,drho=1,iloginfo=0,aggressive=0;
  double penalty=1e8,rho=4,zbar=1e10,cc=0,r0=-1,mu0=-1,ylow,yhigh,gaptol=1e-6,pnormtol=1e30;
  double maxtrust=1e30,steptol=0.01,inftol=1e-8,lpb=1.0,dbound=1e20,infptol=1e-4;
  double dtmp,pstep,dstep,pnorm,mu;
  double *blockn,datanorm[3];
  double *aval,*aval2,*bval,*yout,*y0=0,*xout,*stat;
  double pobj,dobj,dinf;
  DSDP dsdp;
  SDPCone  sdpcone=0;
  DSDPTerminationReason reason;
  DSDPSolutionType pdfeasible;
  LPCone   lpcone=0;
  BCone bcone=0;
  char conetype[30];
  int nfields=25;
  const char *fnames[25]={"stype","obj","pobj","dobj","stopcode","termcode","iterates","r","mu",
			  "pstep","dstep","pnorm","gaphist","infeashist","errors",
			  "datanorm","ynorm","boundy","penalty","tracex","reuse","rho","xy","xdy","xmu"};
  
  if (nrhs < 2){
    mexErrMsgTxt("Two input arguments required.  See help for details. ");}
  if (nrhs > 4){
    mexErrMsgTxt("Fewer input arguments required.  See help for details. ");}
  if (nlhs < 2){
    mexErrMsgTxt("Two output arguments required.  See help for details. ");}
  if (nlhs > 3){
    mexErrMsgTxt("Fewer output arguments required.  See help for details. ");}

  if (!mxIsDouble(B_IN) || mxIsSparse(B_IN)){
    mexErrMsgTxt("DSDP: 1ST input must be a dense vector of doubles"); }
  nvars = mxGetM(B_IN);
  nb = mxGetN(B_IN);    
  if (nb > 1){
    mexErrMsgTxt("DESP: 1ST input must be a column vector"); }
  
  iscellCA = mxIsCell(CA_IN); 
  if (!iscellCA){
    mexErrMsgTxt("DSDP: 2ND input must be a cell array"); }
  mC = mxGetM(CA_IN);
  nC = mxGetN(CA_IN); 
  if (nC != 3){
    mexErrMsgTxt("DSDP: dimension of 2ND cell array is p x 3");}

  if (nrhs >2){
    if(!mxIsStruct(PARS_IN)){
      mexErrMsgTxt("3RD input `OPTIONS' should be a structure.");
    }
  }

  if (nrhs>3){
    if (!mxIsDouble(Y_IN) || mxIsSparse(Y_IN)){
      mexErrMsgTxt("DSDP: 4TH input must be a dense vector of doubles"); }
    m1 = mxGetM(Y_IN);
    n1 = mxGetN(Y_IN);
    if (m1 != nvars || n1 != nb){
      mexErrMsgTxt("DSDP: dimensions of 1ST and 4TH input not compatible");}
    y0=mxGetPr(Y_IN);
    if (!y0){
       mexErrMsgTxt("DSDP: Cannot read 4TH argument");}
  } else y0=0;
  
  
  /* Check data */
  for (j=0; j<mC; j++){
    subs[0] = j; subs[1] = 0;
    index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
    CA_cell_pr = mxGetCell(CA_IN,index); 
    if (!CA_cell_pr){
      mexPrintf("??? DSDP DATA ERROR: Cell: %d, Column: %d\n",j+1,1);
      mexErrMsgTxt("DSDP: Empty Cell. Missing String"); }
    if (!mxIsChar(CA_cell_pr)){
      mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
      mexErrMsgTxt("DSDP: First column of cells in 2ND argument are a string to determine cone type."); }
    mxGetString(CA_cell_pr,conetype,20); 
    
    subs[0] = j; subs[1] = 1;
    index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
    CA_cell_pr = mxGetCell(CA_IN,index); 
    if (!CA_cell_pr){
      mexPrintf("??? DSDP DATA ERROR: Cell: %d, Column: %d\n",j+1,2);
      mexErrMsgTxt("DSDP: Empty Cell. Provide dimension of block."); }
    if (mxIsSparse(CA_cell_pr)){
      mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
      mexErrMsgTxt("DSDP: Second column in 2ND argument must be a dense array of scalars that specify dimension."); }
    if (!mxIsDouble(CA_cell_pr)){
      mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
      mexErrMsgTxt("DSDP: Second column in 2ND argument must specify dimension."); }
    aval=mxGetPr(CA_cell_pr);

    subs[0] = j; subs[1] = 2;
    index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
    CA_cell_pr = mxGetCell(CA_IN,index); 
    if (!CA_cell_pr){
      mexPrintf("??? DSDP DATA ERROR: Cell: %d, Column: %d\n",j+1,3);
      mexErrMsgTxt("DSDP: Empty Cell. Provide sparse data matrix."); }
    if (!mxIsSparse(CA_cell_pr) || !mxIsDouble(CA_cell_pr)){ 
      mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
      mexErrMsgTxt("DSDP: Third column in 2ND argument must be a real sparse data matrix."); }
    
    if (strcmp(conetype,"SDP")==0){ 
      subs[0] = j; subs[1] = 1;
      index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
      CA_cell_pr = mxGetCell(CA_IN,index); 
      aval=mxGetPr(CA_cell_pr);
      it1 = mxGetM(CA_cell_pr);
      it2 = mxGetN(CA_cell_pr); 
      if (it1!=1 && it2!=1){
	mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
 	mexErrMsgTxt("DSDP: Use a dense row vector in the second column in 2ND of the argument.");}

      subs[0] = j; subs[1] = 2;
      index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
      CA_cell_pr = mxGetCell(CA_IN,index); 
      m1 = mxGetN(CA_cell_pr)-1;
      if (m1 != nvars){
	mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
 	mexErrMsgTxt("DSDP: The matrix in the third column in 2ND argument must have number of columns equal to number of variables+1.");}
      vecn = mxGetM(CA_cell_pr); 
      for (tnnz=0,i=0;i<it1*it2;i++){n=(int)aval[i]; tnnz=tnnz+n*(n+1)/2;}
      if ( tnnz != vecn){
	mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
 	mexErrMsgTxt("DSDP: Check Dimensions:  The columns of A and C cannot be converted into square matrices");}
      nsdpblocks=nsdpblocks+it1*it2;
    } else if (strcmp(conetype,"LP")==0){ 

      subs[0] = j; subs[1] = 1;
      index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
      CA_cell_pr = mxGetCell(CA_IN,index); 
      it1 = mxGetM(CA_cell_pr);
      it2 = mxGetN(CA_cell_pr); 
      /* THe sum of the dimensions should equal the number of constraints */
      
      subs[0] = j; subs[1] = 2;
      index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
      CA_cell_pr = mxGetCell(CA_IN,index); 
      if (!mxIsSparse(CA_cell_pr)){ 
	mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
	mexErrMsgTxt("DSDP: Matrices in the third column in 2ND argument must be sparse."); }     
      m1 = mxGetN(CA_cell_pr)-1;
      if (m1 != nvars){
	mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
 	mexErrMsgTxt("DSDP: The matrix in the third column in 2ND argument must have number of columns equal to number of variables+1.");}
    } else if (strcmp(conetype,"LB")==0 || strcmp(conetype,"UB")==0 ){ 
      
      subs[0] = j; subs[1] = 2;
      index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
      CA_cell_pr = mxGetCell(CA_IN,index); 
      if (mxIsSparse(CA_cell_pr)){ 
	mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
	mexErrMsgTxt("DSDP: Row vector in the third column in 2ND argument must be full."); }     
      it1 = mxGetM(CA_cell_pr);
      it2 = mxGetN(CA_cell_pr); 
      if (it2 != 1){
	mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
 	mexErrMsgTxt("DSDP: The matrix in the third column in 2ND argument must have a single column of bounds.");}
      if (it1 != nvars){
	mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
 	mexErrMsgTxt("DSDP: The column matrix in the third column in 2ND argument must contain a bound for each y variable.");}

    } else if (strcmp(conetype,"FIXED")==0){ 
      int dim1;
      subs[0] = j; subs[1] = 1;
      index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
      CA_cell_pr = mxGetCell(CA_IN,index); 
      if (mxIsSparse(CA_cell_pr)){ 
	mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
	mexErrMsgTxt("DSDP: Vector in the third column in 2ND argument must be full."); }     
      it1 = mxGetM(CA_cell_pr);
      it2 = mxGetN(CA_cell_pr); 
      dim1=it1*it2;
      if (it1 != 1){
	mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
 	mexErrMsgTxt("DSDP: Third column in 2ND argument must have 1 row,");}
      subs[0] = j; subs[1] = 2;
      index = mxCalcSingleSubscript(CA_IN,nsubs,subs);
      CA_cell_pr = mxGetCell(CA_IN,index); 
      if (mxIsSparse(CA_cell_pr)){ 
	mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
	mexErrMsgTxt("DSDP: Vector in the third column in 2ND argument must be full."); }
      it1 = mxGetM(CA_cell_pr);
      it2 = mxGetN(CA_cell_pr); 
      if (it1 != 1){
	mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
 	mexErrMsgTxt("DSDP: Third column in 2ND argument must have 1 row,");}
      if (it2 != dim1){
	mexPrintf("??? DSDP DATA ERROR: Cell: %d \n",j+1);
 	mexErrMsgTxt("DSDP: Secord and third column must have same dimension,");}
    } else {
      mexPrintf("??? DSDP DATA ERROR: Cell: %d, Conetype: %s \n",j+1,conetype);
      mexErrMsgTxt("DSDP: Unknown Cone type in 2ND argument. Try 'SDP' or 'LP' or 'Bounds'. ");
    }
  }
  
  
  /* Create output arrays */
  if (nlhs>2){
    if (X_OUT != NULL) mxDestroyArray(X_OUT) ;
    X_OUT = mxCreateCellMatrix(mC,1);
  }
  if (Y_OUT != NULL) mxDestroyArray(Y_OUT) ;
  Y_OUT = mxCreateDoubleMatrix(nvars, 1, mxREAL) ;
  
  /* Create the Solver */
  info = DSDPCreate(nvars,&dsdp); CHKERR(info);
  info = DSDPCreateSDPCone(dsdp,nsdpblocks,&sdpcone); CHKERR(info);

  /* Set Dual Objective Vector */
  bval=mxGetPr(B_IN);
  if (!bval){ mexErrMsgTxt("DSDP: Problems with 1ST argument");}
  for (i=0;i<nvars;i++){info=DSDPSetDualObjective(dsdp,i+1,bval[i]); CHKERR(info);}
  
  /* Set Matrix Data */
  for (j=0; j<mC; j++){  /* Begin Block */
    subs[0] = j; subs[1] = 0;
    index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
    CA_cell_pr = mxGetCell(CA_IN,index); 
    mxGetString(CA_cell_pr,conetype,20); 
    if (strcmp(conetype,"SDP")==0){ 

      subs[0] = j; subs[1] = 1;
      index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
      CA_cell_pr = mxGetCell(CA_IN,index); 
      it1 = mxGetM(CA_cell_pr);
      it2 = mxGetN(CA_cell_pr); 
      blockn=mxGetPr(CA_cell_pr);
      nsubblocks=it1*it2;

      subs[0] = j; subs[1] = 2;
      index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
      CA_cell_pr = mxGetCell(CA_IN,index); 
      aval=mxGetPr(CA_cell_pr); air =mxGetIr(CA_cell_pr); ajc =mxGetJc(CA_cell_pr);
      if (!aval||!air||!ajc)
	{ mexErrMsgTxt("DSDP: Problems with 2ND argument");}
      for (tnnz=0,jj=0;jj<nsubblocks;jj++){n=(int)blockn[jj];tnnz+=n*(n+1)/2;}
      if (nlhs>2){
	subs[0] = j; subs[1] = 0;
	index = mxCalcSingleSubscript(X_OUT,nsubs,subs);
	X_cell_pr = mxCreateDoubleMatrix(tnnz,1,mxREAL);
	mxSetCell(X_OUT,index,X_cell_pr);
	xout=mxGetPr(X_cell_pr);
	if (tnnz>0 && !xout){ mexErrMsgTxt("DSDP: Cannot create array. Out of Memory");}
      }
      
      for (ii=0; ii<=nvars; ii++){ /* Begin Variable matrix constraints */
	i=ii+1;
	if (i==nvars+1) i=0;
	tnnz=0; spot=ajc[ii]; blockj=sdpblockj;
	for (jj=0;jj<nsubblocks;jj++){
	  n=(int)blockn[jj];
	  if (sdpnmax<n) sdpnmax=n;
	  if (ii==0){
	    nn+=n;
	    info=SDPConeSetBlockSize(sdpcone,blockj,n); CHKERR(info);
	    info=SDPConeUsePackedFormat(sdpcone,blockj); CHKERR(info);

	    if (nlhs>2){info=SDPConeSetXArray(sdpcone,blockj,n,xout+tnnz,(n*n+n)/2); CHKERR(info);} 
	    info=CountNonzeroMatrices(blockn,nsubblocks,jj,nvars, air, ajc, &nzmats); CHKERR(info);
	    info=SDPConeSetSparsity(sdpcone,blockj,nzmats); CHKERR(info);
	    if (stat1<nzmats)stat1=nzmats;
	  }
	  for (tnnz2=tnnz+n*(n+1)/2,ijnnz=0;ijnnz<ajc[ii+1]-spot && air[spot+ijnnz]<tnnz2;ijnnz++){}
	  if ( ijnnz==0 ){     /*  info=DSDPSetZeroMat(dsdp,sdpblockj,i,n); */
	  } else if(ijnnz==n*(n+1)/2){ /* check for dense matrix  */
	    if (CheckForConstantMat(aval+spot,ijnnz,n)){
	      info=SDPConeSetConstantMat(sdpcone,blockj,i,n,aval[spot]); CHKERR(info);
	    } else {
	      info=SDPConeSetADenseVecMat(sdpcone,blockj,i,n,1.0,aval+spot,ijnnz); CHKERR(info);
	    }
	  } else {
	    info=SDPConeSetASparseVecMat(sdpcone,blockj,i,n,1.0,tnnz,air+spot,aval+spot,ijnnz); CHKERR(info);
	  }
	  tnnz+=n*(n+1)/2; spot+=ijnnz; blockj++;
	}
      }  /* End Matrices in block */
      sdpblockj=sdpblockj+nsubblocks;

    } else if (strcmp(conetype,"LB")==0 || strcmp(conetype,"UB")==0 ){
      subs[0] = j; subs[1] = 2;
      index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
      CA_cell_pr = mxGetCell(CA_IN,index);
      aval=mxGetPr(CA_cell_pr);
      info=DSDPCreateBCone(dsdp,&bcone);  CHKERR(info);
      info=BConeAllocateBounds(bcone,nvars); CHKERR(info);
      for (i=0;i<nvars;i++){
	if (strcmp(conetype,"LB")==0){
	  info=BConeSetLowerBound(bcone,i+1,aval[i]); CHKERR(info);
	} else {
	  info=BConeSetUpperBound(bcone,i+1,aval[i]); CHKERR(info);
	} 
      }
      if (nlhs>2){
	subs[0] = j; subs[1] = 0;
	index = mxCalcSingleSubscript(X_OUT,nsubs,subs); 
	X_cell_pr = mxCreateDoubleMatrix(nvars,1,mxREAL);
	mxSetCell(X_OUT,index,X_cell_pr);
	aval2=mxGetPr(X_cell_pr);
	info=BConeSetXArray(bcone,aval2,nvars); CHKERR(info);
      }
    } else if (strcmp(conetype,"LP")==0){
      subs[0] = j; subs[1] = 2;
      index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
      CA_cell_pr = mxGetCell(CA_IN,index); 
      n = mxGetM(CA_cell_pr); 
      if (lpnmax<n) lpnmax=n;
      nn+=n;
      aval=mxGetPr(CA_cell_pr); air =mxGetIr(CA_cell_pr); ajc =mxGetJc(CA_cell_pr);
      if (!aval||!air||!ajc)
	{ mexErrMsgTxt("DSDP: Problems with 2ND argument");}
      info=DSDPCreateLPCone(dsdp,&lpcone); CHKERR(info);
      info=LPConeSetData2(lpcone,n,ajc,air,aval); CHKERR(info);
      if (nlhs>2){
	subs[0] = j; subs[1] = 0;
	index = mxCalcSingleSubscript(X_OUT,nsubs,subs); 
	X_cell_pr = mxCreateDoubleMatrix(n,1,mxREAL);
	mxSetCell(X_OUT,index,X_cell_pr);
	xout=mxGetPr(X_cell_pr);
	if (n>0 && !xout){ mexErrMsgTxt("DSDP: Cannot create array. Out of Memory");}
	info=LPConeSetXVec(lpcone,xout,n); CHKERR(info);
      }
    
    } else if (strcmp(conetype,"FIXED")==0){
      int vari;
      subs[0] = j; subs[1] = 1;
      index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
      CA_cell_pr = mxGetCell(CA_IN,index);
      aval=mxGetPr(CA_cell_pr);
      it1 = mxGetM(CA_cell_pr);
      it2 = mxGetN(CA_cell_pr); 
      subs[0] = j; subs[1] = 2;
      index = mxCalcSingleSubscript(CA_IN,nsubs,subs); 
      CA_cell_pr = mxGetCell(CA_IN,index);
      aval2=mxGetPr(CA_cell_pr);
      if (nlhs>2){
	subs[0] = j; subs[1] = 0;
	index = mxCalcSingleSubscript(X_OUT,nsubs,subs); 
	X_cell_pr = mxCreateDoubleMatrix(it1*it2,1,mxREAL);
	mxSetCell(X_OUT,index,X_cell_pr);
	xout=mxGetPr(X_cell_pr);
	if (it1*it2>0 && !xout){ mexErrMsgTxt("DSDP: Cannot create array. Out of Memory");}
      }  else {xout=0;}
      info=DSDPSetFixedVariables(dsdp,aval,aval2,xout,it1*it2); CHKERR(info);
      for (i=0;i<it1*it2;i++){
	/*
	vari=(int)aval[i];
	printf("FixedVariable %d to %4.4e\n",vari,aval2[i]);
       	info=DSDPSetFixedVariable(dsdp,vari,aval2[i]); CHKERR(info); 
	*/
      }
    }    
  } /* End Block */

  /* Set initial point */
  if (y0){
    for (i=0;i<nvars;i++){ info = DSDPSetY0(dsdp,i+1,y0[i]); CHKERR(info);}
  } 

  reuse=(nvars-2)/sdpnmax; 
  if (nvars<50 && reuse==0) reuse=1;
  if (reuse>=1) reuse++;
  reuse=reuse*reuse;
  if (nvars<2000 && reuse>10) reuse=10;
  if (reuse>12) reuse=12;
  info=DSDPReuseMatrix(dsdp,reuse); CHKERR(info);

  info=DSDPGetPPObjective(dsdp,&zbar); CHKERR(info);
  info=DSDPGetR(dsdp,&r0); CHKERR(info);
  info=DSDPGetMaxIts(dsdp,&maxit);  CHKERR(info);
  info=DSDPGetPenaltyParameter(dsdp,&penalty);  CHKERR(info);
  info=DSDPGetPotentialParameter(dsdp,&rho);  CHKERR(info);
  info=DSDPGetDualBound(dsdp,&dbound);  CHKERR(info);
  info=DSDPGetGapTolerance(dsdp,&gaptol); CHKERR(info);
  info=DSDPGetRTolerance(dsdp,&inftol);  CHKERR(info);
  info=DSDPGetBarrierParameter(dsdp,&mu0); CHKERR(info);
  info=DSDPGetMaxTrustRadius(dsdp,&maxtrust);  CHKERR(info);
  info=DSDPGetStepTolerance(dsdp,&steptol);  CHKERR(info);
  info=DSDPGetPTolerance(dsdp,&infptol);  CHKERR(info);
  info=DSDPGetDataNorms(dsdp, datanorm); CHKERR(info);
  info=DSDPGetYBounds(dsdp,&ylow,&yhigh); CHKERR(info);
  if (datanorm[0]==0){info=DSDPSetYBounds(dsdp,-1.0,1.0);CHKERR(info);}
  if (nrhs >2){
    if(!mxIsStruct(PARS_IN)){
      mexErrMsgTxt("Fifth Parameter `OPTIONS' should be a structure.");}

      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"maxit") ){
	maxit=(int) mxGetScalar(OPTIONS_FIELD);
	info=DSDPSetMaxIts(dsdp,maxit); CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"fastblas") ){
	fastblas= (int) mxGetScalar(OPTIONS_FIELD);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"print") ){
	print_info= (int) mxGetScalar(OPTIONS_FIELD);
	printlevel=print_info;
      }
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"bigM") ){
	rpos=(int)mxGetScalar(OPTIONS_FIELD);
	info=DSDPUsePenalty(dsdp,rpos); CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"penalty") ){
	penalty=mxGetScalar(OPTIONS_FIELD);
	info=DSDPSetPenaltyParameter(dsdp,penalty); CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"rho") ){
	rho= mxGetScalar(OPTIONS_FIELD);
	info=DSDPSetPotentialParameter(dsdp,rho); CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"dynamicrho") ){
	drho= (int)mxGetScalar(OPTIONS_FIELD);
	info=DSDPUseDynamicRho(dsdp,drho); CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"zbar") ){
	zbar= mxGetScalar(OPTIONS_FIELD);
	info=DSDPSetZBar(dsdp,zbar); CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"dual_bound") ){
	dbound= mxGetScalar(OPTIONS_FIELD);
	info=DSDPSetDualBound(dsdp,dbound); CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"reuse") ){
	reuse= (int)mxGetScalar(OPTIONS_FIELD);
	info=DSDPReuseMatrix(dsdp,reuse); CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"gaptol") ){
	gaptol= mxGetScalar(OPTIONS_FIELD);
	info=DSDPSetGapTolerance(dsdp,gaptol);CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"lp_barrier") ){
	lpb= mxGetScalar(OPTIONS_FIELD);
	if (lpb<0.1) lpb=0.1;
	if (lpcone){info=LPConeScaleBarrier(lpcone,lpb); CHKERR(info);}}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"lpb") ){
	lpb= mxGetScalar(OPTIONS_FIELD);
	if (lpb<0.1) lpb=0.1;
	if (lpcone){info=LPConeScaleBarrier(lpcone,lpb); CHKERR(info);}}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"cc") ){
	cc= mxGetScalar(OPTIONS_FIELD);
	info=DSDPAddObjectiveConstant(dsdp,cc); CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"inftol") ){
	inftol= mxGetScalar(OPTIONS_FIELD);
	info=DSDPSetRTolerance(dsdp,inftol); CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"infptol") ){
	infptol= mxGetScalar(OPTIONS_FIELD);
	info=DSDPSetPTolerance(dsdp,infptol);CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"pnormtol") ){
	pnormtol= mxGetScalar(OPTIONS_FIELD);
	info=DSDPSetPNormTolerance(dsdp,pnormtol);CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"boundy") ){
	yhigh= fabs(mxGetScalar(OPTIONS_FIELD)); ylow=-yhigh;
	info=DSDPSetYBounds(dsdp,ylow,yhigh);CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"r0") ){
	r0= mxGetScalar(OPTIONS_FIELD); 
	info=DSDPSetR0(dsdp,r0); CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"mu0") ){
	mu0=mxGetScalar(OPTIONS_FIELD);
	info=DSDPSetBarrierParameter(dsdp,mu0);CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"maxtrustradius")){
	maxtrust= mxGetScalar(OPTIONS_FIELD);
	info=DSDPSetMaxTrustRadius(dsdp,maxtrust); CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"steptol") ){
	steptol = mxGetScalar(OPTIONS_FIELD);
	info=DSDPSetStepTolerance(dsdp,steptol); CHKERR(info);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"dobjmin") ){
	dtmp= mxGetScalar(OPTIONS_FIELD);
	info = DSDPSetDualLowerBound(dsdp,dtmp);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"dloginfo") ){
	iloginfo= (int) mxGetScalar(OPTIONS_FIELD);
	info=DSDPLogInfoAllow(iloginfo,0);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"logtime") ){
	printsummary= (int) mxGetScalar(OPTIONS_FIELD);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"printproblem") ){
	info=DSDPPrintData(dsdp,sdpcone,lpcone);}
      if ( OPTIONS_FIELD = mxGetField(PARS_IN,0,"xmaker") ){
	xmaker = (int) mxGetScalar(OPTIONS_FIELD);}
      /*
	if( OPTIONS_FIELD = mxGetField(PARS_IN,0,"maxlanczos") ){
	itmp= (int) mxGetScalar(OPTIONS_FIELD);
	info=DSDPSetLanczosIterations(dsdp,itmp);}
      */
      
  }

  info = DSDPSetMonitor(dsdp,DSDPPrintStats2,0); CHKERR(info);

  info = DSDPSetup(dsdp);
  if (info){
    mexErrMsgTxt("DSDP: Setup Error, Probably out of memory");}


  info = DSDPSolve(dsdp); CHKERR(info);
  info = DSDPStopReason(dsdp,&reason); 
  if (reason!=DSDP_INFEASIBLE_START){
    info=DSDPComputeX(dsdp);CHKERR(info);
  }
  info = DSDPGetSolutionType(dsdp,&pdfeasible);  CHKERR(info);

  if (printsummary){ DSDPEventLogSummary();}

  if (info){
    mexErrMsgTxt("DSDP: Numerical error");}

  if ( reason == DSDP_INFEASIBLE_START){
    mexErrMsgTxt("DSDP Terminated Due to Infeasible Starting Point\n");
  } else if (print_info){

    if (reason == DSDP_CONVERGED)
      mexPrintf("DSDP Converged. \n");
    else if ( reason == DSDP_UPPERBOUND )
      mexPrintf("DSDP Converged: Dual Objective exceeds its bound\n");
    else if ( reason == DSDP_SMALL_STEPS )
      mexPrintf("DSDP Terminated Due to Small Steps\n");
    else if ( reason == DSDP_MAX_IT)
      mexPrintf("DSDP Terminated Due Maximum Number of Iterations\n");
    else if ( reason == DSDP_USER_TERMINATION)
      mexPrintf("DSDP Terminated By User\n");
    else if ( reason == DSDP_INFEASIBLE_START)
      mexPrintf("DSDP Terminated Due to Infeasible Starting Point\n");
    else 
      mexPrintf("DSDP Finished.\n");
  }

  if (pdfeasible==DSDP_UNBOUNDED){
    mexPrintf("DSDP: Dual Unbounded, Primal Infeasible\n");
  } else if (pdfeasible==DSDP_INFEASIBLE){
    mexPrintf("DSDP: Primal Unbounded, Dual Infeasible\n");
  }

  /* Set the dual solution */  
  yout=mxGetPr(Y_OUT);
  info = DSDPGetY(dsdp,yout,nvars); CHKERR(info);
  if (info){
    mexErrMsgTxt("DSDP: Numerical error");}
  
  
  /* Output statistics */
  if (STAT_OUT != NULL) mxDestroyArray(STAT_OUT) ;
  subs[0] = 1;  subs[1] = 1;
  STAT_OUT = mxCreateStructArray(2,subs,nfields,fnames);
  info= DSDPGetDObjective(dsdp,&dobj);  CHKERR(info);
  info= DSDPGetPObjective(dsdp,&pobj);  CHKERR(info);
  info= DSDPGetR(dsdp,&dinf);  CHKERR(info);
  info= DSDPStopReason(dsdp,&reason);  CHKERR(info);
  info= DSDPGetIts(dsdp,&its);  CHKERR(info);
  info= DSDPGetBarrierParameter(dsdp,&mu); CHKERR(info);
  info= DSDPGetStepLengths(dsdp,&pstep,&dstep); CHKERR(info);
  info= DSDPGetPnorm(dsdp,&pnorm); CHKERR(info);
  info= DSDPGetBarrierParameter(dsdp,&mu); CHKERR(info);
  info= DSDPGetYBounds(dsdp,&ylow,&yhigh); CHKERR(info);

  if (pdfeasible==DSDP_UNBOUNDED){
    STAT_FIELD = mxCreateString("Unbounded");
    mxSetField(STAT_OUT,0,fnames[0],STAT_FIELD);
  } else if (pdfeasible==DSDP_INFEASIBLE){
    STAT_FIELD = mxCreateString("Infeasible");
    mxSetField(STAT_OUT,0,fnames[0],STAT_FIELD);
  } else {
    STAT_FIELD = mxCreateString("PDFeasible");
    mxSetField(STAT_OUT,0,fnames[0],STAT_FIELD);
  }

  STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); stat[0]=dobj;
  mxSetField(STAT_OUT,0,fnames[1],STAT_FIELD);

  STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); stat[0]=pobj;
  mxSetField(STAT_OUT,0,fnames[2],STAT_FIELD);
  
  STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); stat[0]=dobj;
  mxSetField(STAT_OUT,0,fnames[3],STAT_FIELD);

  if (reason==DSDP_CONVERGED){
    STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
    stat=mxGetPr(STAT_FIELD); stat[0]=0;
    mxSetField(STAT_OUT,0,fnames[4],STAT_FIELD);
  } else {
    STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
    stat=mxGetPr(STAT_FIELD); stat[0]=(double)reason;
    if (stat[0]==0) stat[0]=-1;
    mxSetField(STAT_OUT,0,fnames[4],STAT_FIELD);
  }
  
  if (reason==DSDP_CONVERGED){ /* For YALMIP */
    STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
    stat=mxGetPr(STAT_FIELD); stat[0]=0;
    if (pdfeasible==DSDP_UNBOUNDED){ stat[0]=1;} 
    if (pdfeasible==DSDP_INFEASIBLE){ stat[0]=2;}
    mxSetField(STAT_OUT,0,fnames[5],STAT_FIELD);
  } else {
    STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
    stat=mxGetPr(STAT_FIELD); stat[0]=-1;
    if (reason==DSDP_INFEASIBLE_START)stat[0]=-6;
    if (reason==DSDP_UPPERBOUND)stat[0]=-5;
    if (reason==DSDP_MAX_IT)stat[0]=-3;
    mxSetField(STAT_OUT,0,fnames[5],STAT_FIELD);
  }

  STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); stat[0]=(double) its;
  mxSetField(STAT_OUT,0,fnames[6],STAT_FIELD);
  
  STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); stat[0]=dinf;
  mxSetField(STAT_OUT,0,fnames[7],STAT_FIELD);
    
  STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); stat[0]=mu;
  mxSetField(STAT_OUT,0,fnames[8],STAT_FIELD);
  
  STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); stat[0]=pstep;
  mxSetField(STAT_OUT,0,fnames[9],STAT_FIELD);
  
  STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); stat[0]=dstep;
  mxSetField(STAT_OUT,0,fnames[10],STAT_FIELD);
  
  STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); stat[0]=pnorm;
  mxSetField(STAT_OUT,0,fnames[11],STAT_FIELD);
  
  itmp=100; if (its < itmp) itmp=its;
  STAT_FIELD = mxCreateDoubleMatrix(itmp, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); 
  info= DSDPGetGapHistory(dsdp,stat,itmp);  CHKERR(info);
  mxSetField(STAT_OUT,0,fnames[12],STAT_FIELD);
  
  STAT_FIELD = mxCreateDoubleMatrix(itmp, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); 
  info= DSDPGetRHistory(dsdp,stat,itmp);  CHKERR(info);
  mxSetField(STAT_OUT,0,fnames[13],STAT_FIELD);

  STAT_FIELD = mxCreateDoubleMatrix(6, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); 
  info = DSDPGetFinalErrors(dsdp,stat);  CHKERR(info);
  mxSetField(STAT_OUT,0,fnames[14],STAT_FIELD);

  STAT_FIELD = mxCreateDoubleMatrix(3, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); 
  info = DSDPGetDataNorms(dsdp,stat);  CHKERR(info);
  dtmp=stat[0];stat[0]=stat[1];stat[1]=dtmp;
  mxSetField(STAT_OUT,0,fnames[15],STAT_FIELD);

  STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); 
  info = DSDPGetYMaxNorm(dsdp,stat);  CHKERR(info);
  mxSetField(STAT_OUT,0,fnames[16],STAT_FIELD);

  STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); 
  stat[0]=yhigh;
  mxSetField(STAT_OUT,0,fnames[17],STAT_FIELD);

  STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); 
  info = DSDPGetPenaltyParameter(dsdp,stat);  CHKERR(info);
  mxSetField(STAT_OUT,0,fnames[18],STAT_FIELD);

  STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); 
  info = DSDPGetTraceX(dsdp,stat);  CHKERR(info);
  mxSetField(STAT_OUT,0,fnames[19],STAT_FIELD);

  STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
  info = DSDPGetReuseMatrix(dsdp,&its);  CHKERR(info);
  stat=mxGetPr(STAT_FIELD); stat[0]=(double) its;
  mxSetField(STAT_OUT,0,fnames[20],STAT_FIELD);
  
  STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
  stat=mxGetPr(STAT_FIELD); 
  info = DSDPGetPotentialParameter(dsdp,stat);  CHKERR(info);
  mxSetField(STAT_OUT,0,fnames[21],STAT_FIELD);

  if (xmaker){
    STAT_FIELD = mxCreateDoubleMatrix(nvars+1, 1, mxREAL) ;
    stat=mxGetPr(STAT_FIELD); 
    info = DSDPGetYMakeX(dsdp,stat,nvars+1);  CHKERR(info);
    mxSetField(STAT_OUT,0,fnames[22],STAT_FIELD);

    STAT_FIELD = mxCreateDoubleMatrix(nvars+1, 1, mxREAL) ;
    stat=mxGetPr(STAT_FIELD); 
    info = DSDPGetDYMakeX(dsdp,stat,nvars+1);  CHKERR(info);
    mxSetField(STAT_OUT,0,fnames[23],STAT_FIELD);

    STAT_FIELD = mxCreateDoubleMatrix(1, 1, mxREAL) ;
    stat=mxGetPr(STAT_FIELD); 
    info = DSDPGetMuMakeX(dsdp,stat);  CHKERR(info);
    mxSetField(STAT_OUT,0,fnames[24],STAT_FIELD);
  }
  /* Free internal data structure */

  info = DSDPDestroy(dsdp); CHKERR(info);

  return;
} /* main */
  


#undef __FUNCT__
#define __FUNCT__ "CheckForConstantMat"
static int CheckForConstantMat(double*v,int nnz, int n){
  int i;double vv;
  if (n<=1){ return 0; }
  if (nnz!=(n*n+n)/2){ return 0; }
  for (vv=v[0],i=1;i<nnz;i++){
    if (v[i]!=vv){ return 0;}
  }
  return 1;
}

#undef __FUNCT__
#define __FUNCT__ "CountNonzeroMatrices"
 static int CountNonzeroMatrices(double *blocksize, int nblocks, int block, int nvars, int *indd, int*nnnz, int *nnzmats){
   int i,j,n,nzmats=0;
   int marker1=0,marker2=0;
  for (i=0;i<block;i++){n=(int)blocksize[i]; marker1=marker1+n*(n+1)/2;}
  n=(int)blocksize[block];  
  marker2=marker1+n*(n+1)/2;
  for (i=0;i<nvars;i++){
    j=nnnz[i];
    while (indd[j]<marker1 && j< nnnz[i+1]){j++;}
    if (j<nnnz[i+1] && indd[j]<marker2){ nzmats++;}
  }
  *nnzmats=nzmats;
  return 0;
}

/* ---------------------------------------------------------- */
 static int DSDPPrintStats2(DSDP dsdp, void* dummy){
   
  int    iter,info;
  double pobj,dobj,pstp=0,dstp,mu,res,pnorm,pinfeas;
  DSDPTerminationReason reason;

  if(printlevel<=0) return(0);

  info = DSDPStopReason(dsdp,&reason);
  info = DSDPGetIts(dsdp,&iter);

  if( (reason!=CONTINUE_ITERATING) || ((iter % printlevel)==0)){

    info = DSDPGetDDObjective(dsdp,&dobj);
    info = DSDPGetPPObjective(dsdp,&pobj);
    info = DSDPGetR(dsdp,&res);
    info = DSDPGetPInfeasibility(dsdp,&pinfeas);
    info = DSDPGetStepLengths(dsdp,&pstp,&dstp);
    info = DSDPGetBarrierParameter(dsdp,&mu);
    info = DSDPGetPnorm(dsdp,&pnorm);
    
    
    if (iter==0){
      mexPrintf("Iter   PP Objective      DD Objective     PInfeas   DInfeas     Nu     StepLength   Pnrm\n")
	;
      mexPrintf("----------------------------------------------------------------------------------------\n")
	;
    }
    mexPrintf("%-3d %16.8e  %16.8e  %9.1e %9.1e %9.1e",iter,pobj,dobj,pinfeas,res,mu);
    mexPrintf("  %4.2f  %4.2f",pstp,dstp);
    if (pnorm>1.0e3){
      mexPrintf("  %1.0e \n",pnorm);
    } else {
      mexPrintf("  %5.2f \n",pnorm);
    }
    fflush(NULL);
  }
  return(0);
}
