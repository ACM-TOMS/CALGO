#include "dsdp5.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "pdsdp5scalapack.h"

/*! \file readsdpa.c 
  \brief Read SDPA data files, pass data into DSDP solver, and print solution.
*/

/*
static char help[]="\n";
*/
static char help[]="\
DSDP Usage: dsdp5 filename <sdpa format file> \n\
                  -print <10> - print information at each k iteration\n\
                  -save <Solution File in SDPA format> \n\
                  -fout <filename> to print standard monitor to a file\n\
                  -y0 <initial solution file> \n\
                  -benchmark <file containing names of SDPA files> \n\
                    -directory <directory containing benchmark SDPA files> \n\
                    -suffix  <add to each benchmark problem name> \n\
                  -dloginfo <0> - print more information for higher numbers \n\
                  -dlogsummary <1> - print timing information \n\
                  -help  for this help message\n";

static int qusort(int[],int[],int[],double[],int,int);
static int partition(int[],int[],int[],double[],int, int);
static int Parseline(char *,int *,int *,int *,int *,double *, int *);
static int ReadInitialPoint(char*, int, double[]);
static int TCheckArgs0(DSDP,SDPCone,int,int,char *[]);
static int TCheckArgs(DSDP,SDPCone,int,int,char *[]);
static int CheckForConstantMat(double[],int, int);
static int CountNonzeroMatrices(int, int[],int[], int*);

typedef struct{
  char sformat;
  int blocksize;
} DBlock;

typedef struct{
  int *block,*constraint,*matind;
  double*nnz;
  char *sformat;
  int totalnonzeros;
  double *dobj,*y0;
  char *conetypes;
  int *blocksizes; 
  int m; int n; int nblocks;
  int lpn,lpspot,lpblock,lpnnz;
  int *lpi,*lui,*cmap;
  double cnorm;
  double fixedvari;
  double fixedvard,xout;
} DSDPData;

static int ReadSDPA2(char*,DSDPData*);
static int GetMarkers(int, int, int*, int*, int*);
static int ComputeY0(DSDP,DSDPData);
static int rank=0;
int ReadSDPAFile(int argc,char *argv[]);

#define CHKDATA(a,b,c)  { if (c){ printf("Possible problem in variable %d, block %d. \n",a+1,b+1);} }


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char *argv[]){
  int info;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  DSDPSetRank(rank);
  info=ReadSDPAFile(argc,argv);
  MPI_Finalize( );
  return info;
}

#undef __FUNCT__
#define __FUNCT__ "ReadSDPAFile"
/*!
\fn int ReadSDPAFile(int argc,char *argv[]);
\brief Read SDPA formatted file and solve the semidefinite program.
\param argc number of command line arguments
\param argv command line arguments
\ingroup Examples
*/
int ReadSDPAFile(int argc,char *argv[]){

  int      i,j,m,n,np,its,info;
  int      spot,ijnnz,nzmats,sdpnmax,sdpn,stat1,printsummary=1;
  int      runbenchmark=0,saveit=0,justone=1,fileout=0;
  double   t1,t2,t3,t4,t5,dd,yhigh; 
  double   derr[6],dnorm[3];
  double   ddobj,ppobj;
  char     problemname[100],thisline[100], filename[300],savefile[100];
  char     directory[100]="/home/benson/sdpexamples/sdplib/";
  char     outputfile[50]="",suffix[20]=".dat-s", tablename[20]="results-dsdp-5.7";
  char     success='f',sformat;
  FILE     *fp1=0,*fp2=0,*fout;
  DSDPData dddd;
  DSDP     dsdp;
  DSDPTerminationReason reason;
  DSDPSolutionType pdfeasible;
  SDPCone  sdpcone=0;
  LPCone   lpcone=0;
  int *ittt,sspot;

  if (argc<2){ printf("%s",help); DSDPPrintOptions(); return(0); }
  for (i=0; i<argc; i++){ if (strncmp(argv[i],"-help",5)==0){
    printf("%s",help); DSDPPrintOptions(); return(0);}}
  for (i=1; i<argc; i++){ printf("%s ",argv[i]); } printf("\n");
  for (i=1; i<argc-1; i++){   /* Are we reading a file with a lot of problem names? */
    if (strncmp(argv[i],"-benchmark",8)==0){
      strncpy(thisline,argv[i+1],90); fp1=fopen(thisline,"r");runbenchmark=1; justone=0;
    };
    if (strncmp(argv[i],"-directory",8)==0){strncpy(directory,argv[i+1],90);}
    if (strncmp(argv[i],"-table",4)==0){strncpy(tablename,argv[i+1],90);};
    if (strncmp(argv[i],"-suffix",4)==0){strncpy(suffix,argv[i+1],20);};
    if (strncmp(argv[i],"-save",5)==0){ strncpy(savefile,argv[i+1],40);saveit=1;};
    if (strncmp(argv[i],"-dlogsummary",8)==0){printsummary=atoi(argv[i+1]);}
    if (rank==0&&strncmp(argv[i],"-fout",5)==0){ strncpy(outputfile,argv[i+1],45);fileout=1;};
  }


  while ((runbenchmark && !feof(fp1)) || justone==1){
    justone=0;
    if (runbenchmark){ /* Get filename with the data */
      fgets(thisline,100,fp1); if (sscanf(thisline,"%s",problemname)<1) break;
      strncpy(filename,directory,90); strncat(filename,problemname,90);strncat(filename,suffix,90);
      printf("%s\n",problemname);
    } else {
      strncpy(filename,argv[1],90);
      strncpy(problemname,argv[1],90);
    }
    
    if (rank==0 && fileout){
      dsdpoutputfile=fopen(outputfile,"a");
      fprintf(dsdpoutputfile,"%s\n",problemname);
    } else {dsdpoutputfile=0;fileout=0;}
    DSDPTime(&t1);

    info=ReadSDPA2(filename, &dddd);  m=dddd.m;
    if (info){  printf("Problem reading SDPA file\n"); return 1;}
    DSDPTime(&t2);
    if (rank==0){
      printf("\nVariables %d \n",dddd.m);
      printf("Matrix Blocks: %d, ",dddd.nblocks);
      printf("Total Number of Constraints: %d \n",dddd.n);
      printf("Nonzeros in Constraints: %d\n\n",dddd.totalnonzeros);
      printf("Read Data File into Buffer:      %4.3e seconds\n",t2-t1);
    }
    
    for (i=0; i<argc-1; i++){ 
      if (strncmp(argv[i],"-dloginfo",8)==0){ info=DSDPLogInfoAllow(atoi(argv[i+1]),0);};
    }

    info = DSDPCreate(dddd.m,&dsdp); 
    info = DSDPCreateSDPCone(dsdp,dddd.nblocks,&sdpcone);
    /* Set Dual objective vector */
    for (i=0;i<m;i++){info = DSDPSetDualObjective(dsdp,i+1,dddd.dobj[i]);}
    
    /* Set  initial point */
    for (i=0; i<m; i++)
      if (dddd.dobj[i]> 0.0) dddd.y0[i]=-0.0; else dddd.y0[i]=0.0;
    for (i=0; i<m; i++){info = DSDPSetY0(dsdp,i+1,dddd.y0[i]);}
    info=ComputeY0(dsdp,dddd);
    if (dddd.fixedvari){
      info = DSDPSetY0(dsdp,(int)dddd.fixedvari,dddd.fixedvard);
      printf("Fixed: %2.0f %4.2f ?\n",dddd.fixedvari,dddd.fixedvard);
      DSDPSetFixedVariables(dsdp,&dddd.fixedvari,&dddd.fixedvard,&dddd.xout,1);
    }

    spot=0;ijnnz=0;np=0;sdpnmax=1;sdpn=0;stat1=1;
    /* Insert the SDP data */
    for (j=0;j<dddd.nblocks; j++){
      if (dddd.conetypes[j]=='S'){
	n=dddd.blocksizes[j];
	sformat=dddd.sformat[j];
	info=CountNonzeroMatrices(j+1,dddd.block+spot,dddd.constraint+spot,&nzmats);
	info=SDPConeSetBlockSize(sdpcone,j,n); DSDPCHKERR(info);
	info=SDPConeSetSparsity(sdpcone,j,nzmats); DSDPCHKERR(info);
	info=SDPConeSetStorageFormat(sdpcone,j,sformat); DSDPCHKERR(info);
	np+=n; sdpn+=n;
	if (sdpnmax<n) sdpnmax=n;
	if (stat1<nzmats) stat1=nzmats;
	for (i=0; i<=m; i++){
	  info=GetMarkers(j+1,i,dddd.block+spot,dddd.constraint+spot,&ijnnz);	  
	  if (0==1){
	  } else if ( ijnnz==0 ){  /* info=DSDPSetZeroMat(dsdp,j,i,n); */
	  } else if (CheckForConstantMat(dddd.nnz+spot,ijnnz,n)){
	    info=SDPConeSetConstantMat(sdpcone,j,i,n,dddd.nnz[spot+1]);CHKDATA(i,j,info);
	    if(sformat=='P'){info=SDPConeSetXArray(sdpcone,j,n,dddd.nnz+spot,n*(n+1)/2);}
	  } else if (sformat=='P' && ijnnz==n*(n+1)/2 ){     /* check for dense matrix  */
	    info=SDPConeSetADenseVecMat(sdpcone,j,i,n,1.0,dddd.nnz+spot,ijnnz);CHKDATA(i,j,info);
	  } else {     /* sparse matrix  */
	    info=SDPConeSetASparseVecMat(sdpcone,j,i,n,1.0,0,dddd.matind+spot,dddd.nnz+spot,ijnnz);CHKDATA(i,j,info);
	  }
	  if (0==1){ info=SDPConeViewDataMatrix(sdpcone,j,i);}
	  spot+=ijnnz;
	  /*	  SDPConeScaleBarrier(sdpcone,j,j+1.0); */
	}
	if (0==1){info=SDPConeView2(sdpcone);}
      }  else if (dddd.conetypes[j]=='L'){
	info=DSDPCreateLPCone(dsdp,&lpcone); sformat='P';
	info=SDPConeSetStorageFormat(sdpcone,j,sformat); DSDPCHKERR(info);
	n=dddd.blocksizes[j];
	np+=n;
	sspot=spot;
	DSDPCALLOC2(&ittt,int,(m+2),&info);
	for (i=0;i<=m;i++){ittt[i]=0;}
	for (i=0;i<=m;i++){
	  info=GetMarkers(j+1,i,dddd.block+spot,dddd.constraint+spot,&ijnnz);
	  ittt[i+1]=ijnnz; spot+=ijnnz;
	}
	for (i=1;i<=m;i++)ittt[i+1]+=ittt[i];
	info=LPConeSetData(lpcone,n,ittt,dddd.matind+sspot,dddd.nnz+sspot);CHKDATA(i,0,info);
	if (0==1){info=LPConeView(lpcone);}
	if (0==1){info=LPConeView2(lpcone);}
	/*	info=DSDPFree(&ittt); */
      }
    }
    if (0==1){
      BCone bcone;
      info=DSDPCreateBCone(dsdp, &bcone); 
      info=BConeAllocateBounds(bcone,2*m);
      for (i=0;i<m;i++){
        info=BConeSetUpperBound(bcone,i+1,10);
      }
      for (i=0;i<m;i++){
        info=BConeSetLowerBound(bcone,i+1,-10);
      }
    }

    DSDPTime(&t3);
    if (rank==0){printf("DSDP Set Data:                   %4.3e seconds\n",t3-t2);}

    its=(m-2)/sdpnmax;
    if (np<100 && its==0) its=1;
    if (dddd.lpn>=m && its==0) its=1; 
    if (its>=1) its++;
    its=its*its;
    if (m<2000 && its>10) its=10;
    if (its>12) its=12;

    info=DSDPReuseMatrix(dsdp,its);

    DSDPFREE(&dddd.blocksizes,&info);
    DSDPFREE(&dddd.sformat,&info);
    DSDPFREE(&dddd.dobj,&info);
    DSDPFREE(&dddd.y0,&info);
    DSDPFREE(&dddd.conetypes,&info);
    DSDPFREE(&dddd.constraint,&info);
    DSDPFREE(&dddd.block,&info);
 
    info=DSDPGetDataNorms(dsdp, dnorm);
    if (dnorm[0]==0){
      info=DSDPSetR0(dsdp,np); 
      info=DSDPSetGapTolerance(dsdp,1e-3);
      info=DSDPSetYBounds(dsdp,-1.0,1.0);
    } else {
    }
    info = TCheckArgs0(dsdp,sdpcone,dddd.m,argc,argv);
    info = TCheckArgs(dsdp,sdpcone,dddd.m,argc,argv);

    info=PDSDPUseSCALAPACKLinearSolver(dsdp);
    info = DSDPSetup(dsdp); if (info){ printf("\nProblem Setting problem.  Likely insufficient memory\n");  return 1;}
    if (0==1){info=SDPConeCheckData(sdpcone);}


    DSDPTime(&t4);
    if (rank==0){
      printf("DSDP Process Data:               %4.3e seconds\n\n",t4-t3);
      printf("Data Norms: C: %4.2e, A: %4.2e, b: %4.2e\n",dnorm[0],dnorm[1],dnorm[2]);
      info=DSDPGetScale(dsdp,&ddobj);
      printf("Scale C: %4.2e\n\n",ddobj);
      info=DSDPGetPotentialParameter(dsdp,&ddobj); 
      printf("Potential Parameter: %4.2f\n",ddobj);
      info=DSDPGetReuseMatrix(dsdp,&its);
      printf("Reapply Schur matrix: %d\n\n",its);
    }
    if (0==1){info=DSDPPrintData(dsdp,sdpcone,lpcone);}

    info = DSDPSolve(dsdp); 
    if (info){ printf("\nNumerical errors encountered in DSDPSolve(). \n");}
    
    info=DSDPStopReason(dsdp,&reason);
    if (reason!=DSDP_INFEASIBLE_START){
      info=DSDPComputeX(dsdp);DSDPCHKERR(info);
    }
    info=DSDPStopReason(dsdp,&reason);
    info=DSDPGetSolutionType(dsdp,&pdfeasible);

    DSDPTime(&t5);

    if (rank==0){

      if (reason == DSDP_CONVERGED){
	printf("DSDP Converged. \n"); 
	success='s';
      } else if ( reason == DSDP_UPPERBOUND ){
	printf("DSDP Terminated Because Dual Objective Exceeded its Bound\n");
	success='s';
      } else if ( reason == DSDP_SMALL_STEPS ){
	printf("DSDP Terminated Due to Small Steps\n");
	success='f';
      } else if ( reason == DSDP_MAX_IT){
	printf("DSDP Terminated Due Maximum Number of Iterations\n");
	success='f';
      } else if ( reason == DSDP_INFEASIBLE_START){
	printf("DSDP Terminated Due to Infeasible Starting Point\n");
	success='f';
      } else if ( reason == DSDP_INDEFINITE_SCHUR_MATRIX){
	printf("DSDP Terminated Due to Indefinite Schur Complement\n");
	success='f';
      } else {
	printf("DSDP Finished\n");
	success='f';
      }

      if (pdfeasible == DSDP_UNBOUNDED ){
	printf("DSDP Dual Unbounded, Primal Infeasible\n"); 
      } else if ( pdfeasible == DSDP_INFEASIBLE ){
	printf("DSDP Primal Unbounded, Dual Infeasible\n");
      }

      
      info=DSDPGetDObjective(dsdp,&ddobj);
      info=DSDPGetPObjective(dsdp,&ppobj);
      printf("\nP Objective  : %16.8e \n",ppobj);
      printf("DSDP Solution: %16.8e \n\n",ddobj);
      printf("DSDP Solve Time:                     %4.3e seconds\n",t5-t4);
      printf("DSDP Preparation and Solve Time:     %4.3e seconds\n\n",t5-t3);

      info=DSDPGetFinalErrors(dsdp,derr);
      info=DSDPGetIts(dsdp,&its);
      if (1 || runbenchmark){
	fp2=fopen(tablename,"a");
	if (pdfeasible==DSDP_UNBOUNDED){
	  fprintf(fp2," %-18s & %4d & %4d &    infeasible &     unbounded & %4.0e & %c & %3d & %6.2f \\\\\n",problemname,m,np,derr[0],success,its,t5-t3);
	} else if (pdfeasible==DSDP_INFEASIBLE){
	  fprintf(fp2," %-18s & %4d & %4d &     unbounded &    infeasible & %4.0e & %c & %3d & %6.2f \\\\\n",problemname,m,np,derr[0],success,its,t5-t3);
	} else { 
	  fprintf(fp2," %-18s & %4d & %4d & %13.7f & %13.7f & %4.0e & %c & %3d & %6.2f \\\\\n",problemname,m,np,-ppobj,-ddobj,derr[0],success,its,t5-t3);
	}
	fclose(fp2);
      }
      
      /*      info=DSDPComputeMinimumXEigenvalue(dsdp,&derr[1]); */
      printf("\nP Infeasible: %8.2e \n",derr[0]);
      printf("D Infeasible: %8.2e \n",derr[2]);
      printf("Minimal P Eigenvalue: %6.2e \n",derr[1]);
      printf("Minimal D Eigenvalue: 0.00 \n");
      printf("Relative P - D Objective values: %4.2e \n",derr[4]);
      printf("Relative X Dot S: %4.2e \n",derr[5]);
      info=DSDPGetYBounds(dsdp,&dd,&yhigh);
      info=DSDPGetYMaxNorm(dsdp,&dd);
      printf("\nMax Y: %10.8e,  Bounded by %6.1e\n",dd,yhigh);
      info=DSDPGetTraceX(dsdp,&dd);
      printf("Trace X: %4.8e,   ",dd);
      info=DSDPGetPenaltyParameter(dsdp,&dd);
      printf("Bounded by Penalty Parameter: %4.1e \n\n",dd);
      
      if (printsummary){ DSDPEventLogSummary();}
      printf("--- DSDP Finished ---\n\n");

      if (saveit){
	fout=fopen(savefile,"w");
	/*	fprintf(fout,"** %s \n",filename); Deleted */
	info= DSDPPrintSolution(fout,dsdp,sdpcone,lpcone);
	if (dddd.fixedvari){
	  sspot=dddd.nblocks+1,dd=dddd.xout;
	  fprintf(fout,"1 %d 1 1 1.0e-11\n1 %d 2 2 1.0e-11\n",sspot,sspot);
	  fprintf(fout,"2 %d 1 1 %12.8e\n",sspot,DSDPMax(1.0e-10,dd));
	  fprintf(fout,"2 %d 2 2 %12.8e\n",sspot,DSDPMax(1e-10,-dd));
	}
	fclose(fout);
      }
    }
    
    info = DSDPDestroy(dsdp);
    DSDPFREE(&dddd.matind,&info);
    DSDPFREE(&dddd.nnz,&info);

    if (fileout){fclose(dsdpoutputfile);}
    if (0){ DSDPMemoryLog();}

  }
  if (runbenchmark){ fclose(fp1);}

  return 0;
} /* main */



#define BUFFERSIZ 4000
#undef __FUNCT__
#define __FUNCT__ "ReadSDPA2"
static int ReadSDPA2(char *filename, DSDPData*ddd){
  FILE*fp;
  char ctmp,refline[BUFFERSIZ]="*",thisline[BUFFERSIZ]="*";
  int info,tline,line=0;
  int i,k,m,n;
  /* int spot,nzmark, */
  int bigint=1000000;
  int i1,nblk,nmat,col,row;
  int np=0,nblocks;
  int nargs,nonzero;
  double val;

  fp=fopen(filename,"r");
  if (!fp){
    printf("Cannot open file %s !",filename); return(1);
  }

  /* Read comments */
  while(!feof(fp) && (thisline[0] == '*' || thisline[0] == '"') ){
    fgets(thisline,BUFFERSIZ,fp); line++;
  }
  /* Read number of constraints */
  if (sscanf(thisline,"%d",&m)<1){
    printf("Error: line %d.  Number of constraints not given.\n",line);
    return(1);
  }
  /* Read number of blocks */
  fgets(thisline,BUFFERSIZ,fp); line++;
  if (sscanf(thisline,"%d",&nblocks)!=1){
    printf("Error: line %d.  Number of blocks not given.\n",line);
    return(1);
  }
  ddd->lpn=0;ddd->lpspot=0;ddd->lpblock=0;ddd->cnorm=0;
  /* Read block sizes */
  DSDPCALLOC2(&ddd->sformat,char, (nblocks+1),&info);
  DSDPCALLOC2(&ddd->blocksizes,int, (nblocks+1),&info);
  DSDPCALLOC2(&ddd->conetypes,char, (nblocks+1),&info );
  line++;
  for (i=0;i<nblocks; i++){
    if (fscanf(fp,"{")==1 || fscanf(fp,"(")==1 || fscanf(fp,",")==1 ){
      i--;
    } else if (fscanf(fp,"%d",&col)==1){ 
      if (col>0) {  ddd->blocksizes[i]=col;  np+=col; ddd->conetypes[i]='S';
      } else if (col>0){ddd->blocksizes[i]=-col;  np+=-col; ddd->conetypes[i]='S';
      } else if (col<0){ddd->blocksizes[i]=-col; np += -col; ddd->conetypes[i]='L';ddd->lpn=-col;ddd->lpblock=i;
      } else { ddd->blocksizes[i]=0; ddd->conetypes[i]='N';}
      if (ddd->blocksizes[i]<10){ddd->sformat[i]='U';} else {ddd->sformat[i]='P';}
    }
    else{ printf("Error block sizes, line %d",line); return(1);}
  }
  if (ddd->blocksizes[nblocks-1]==0) nblocks--;
  fgets(thisline,BUFFERSIZ,fp); 
  
  /* Read objective vector */
  DSDPCALLOC2(&ddd->y0,double,m,&info);
  DSDPCALLOC2(&ddd->dobj,double,m,&info);
  line++;
  for (i=0;i<m;i++){
    if (fscanf(fp,",")==1){
      i--;
      continue;
    }
    while (fscanf(fp,"%lg",&val)!=1){
      fscanf(fp,"%c",&ctmp);
      if (ctmp=='\n'){
	printf("Constraints: %d, Blocks: %d\n",m,nblocks);
	printf("Error reading objective, line %d, i=%d \n",line,i); return 1;
      }
    }
    ddd->dobj[i]=val;
  }
  fgets(thisline,BUFFERSIZ,fp); 
  tline=line;

  nargs=5;  nonzero=0;
  while(!feof(fp)){
    thisline[0]='\0';
    nmat=-1; nblk=-1; row=-1; col=-1; val=0.0;
    fgets(thisline,BUFFERSIZ,fp); line++;
    info = Parseline(thisline,&nmat,&nblk,&row,&col,&val,&nargs); 
    if (!feof(fp)&&nargs!=5&&nargs>0){
      printf("Error: line %d \n%s\n",line,thisline);return 1;}
    if (nargs==5 && val!=0.0){
      nonzero++;
      i1=row*(row+1)/2 + col;
      if (row >= ddd->blocksizes[nblk-1] || col >= ddd->blocksizes[nblk-1] ) { 
	printf("Data Error in line: %d.  Row %d or col %d > blocksize %d\n%s",line,row+1,col+1,ddd->blocksizes[nblk-1],thisline);
	return 1;
      }
      if (row<0 || col<0){
	printf("Data Error in line: %d.  Row %d or col %d <= 0 \n%s",line,row+1,col+1,thisline);
	return 1;
      }
      if (nmat>m || nmat<0){
	printf("Data Error in line: %d.  Is Var  0 <= %d <= %d \n%s",line,nmat,m,thisline);
	return 1;
      }
      if (nblk>nblocks || nblk<0){
	printf("Data Error in line: %d.  Is Block  0 <= %d <= %d \n%s",line,nmat,m,thisline);
	return 1;
      }
    } else if (nargs==5 && val==0.0){
    } else if (nargs==0){
    } else {
      printf("Problem Reading SDPA file at line %d:  %s\n",line, thisline);
    }
  }

  /* Allocate memory for the data */
  nonzero++;
  DSDPCALLOC2(&ddd->matind,int,nonzero,&info);
  DSDPCALLOC2(&ddd->nnz,double,nonzero,&info);
  DSDPCALLOC2(&ddd->block,int,nonzero,&info);
  DSDPCALLOC2(&ddd->constraint,int,nonzero,&info);
  nonzero--;

  fseek(fp,0,SEEK_SET);
  line=0;
  for (i=0;i<tline;i++){ctmp='*'; while (ctmp!='\n') fscanf(fp,"%c",&ctmp); line++;}

  nargs=5;k=0;
  while(!feof(fp) /* && nargs==5 */){  
    thisline[0]='\0';
    fgets(thisline,BUFFERSIZ,fp);
    if (k==0){strncpy(refline,thisline,BUFFERSIZ-1); }
    info = Parseline(thisline,&nmat,&nblk,&row,&col,&val,&nargs); 
    if (!feof(fp)&&nargs!=5&&nargs<0){
      /* if (k>=nonzero && !feof(fp) ){ */
      printf("Problem Reading SDPA file at line %d:  %s\n",line, thisline);
      printf("Problem could be earlier in file \n"); 
      printf("The first recorded matix nonzero in the file occured at line %d: \n %s",tline,refline); 
      printf(" Please check data file\n"); 
      return 1;
    }
    if (nargs==5 && val!=0.0){
      if (row>col){
	printf("Warning: Line: %d Row < Column. %s \n",line,thisline);
      }
      i=row;row=col;col=i;
      n=ddd->blocksizes[nblk-1];
      if (nmat==0) {val=-val;}
      if (ddd->conetypes[nblk-1]=='S'){
	/*	if (row==col) val/=2; */
	ddd->matind[k]=row*(row+1)/2 + col;
	if (ddd->sformat[nblk-1]=='U'){ddd->matind[k]=row*n + col;}
      } else {
	ddd->matind[k]=col;
      }
      ddd->block[k]=nblk;
      ddd->constraint[k]=nmat;
      ddd->nnz[k]=val;
      k++;
    } else if (nargs==5 && val==0.0){
    } else if (nargs==0){
    } else {
      printf("Problem Reading SDPA file at line %d:  %s\n",line, thisline);
      printf("Problem could be earlier in file \n"); 
      printf("The first recorded matix nonzero in the file occured at line %d: \n %s",tline,refline); 
      printf(" Please check data file\n"); 
      return 1;
    }
  }
  ddd->block[k]=nblocks+1;  ddd->constraint[k]=m+2;
  ddd->matind[k]=10000000;  ddd->nnz[k]=0.0;

  qusort(ddd->block,ddd->constraint,ddd->matind,ddd->nnz,0,nonzero-1);

  for (i=0;i<nonzero-1; i++){
    while (i<nonzero-1 && ddd->matind[i]==ddd->matind[i+1] && 
	   ddd->constraint[i]==ddd->constraint[i+1] && 
	   ddd->block[i]==ddd->block[i+1] ){
      printf("DSDPError:   Reading Input File:\n");
      printf("Possible problem with data input file: Double Entry: \n");
      printf("   %d %d %d %2.10e\n",
	     ddd->constraint[i],ddd->block[i],ddd->matind[i]+1,ddd->nnz[i]);
      printf("   %d %d %d %2.10e\n\n",
	     ddd->constraint[i+1],ddd->block[i+1],ddd->matind[i+1]+1,ddd->nnz[i+1]);
      for (k=i+1;k<nonzero-1;k++){
	ddd->constraint[k]=ddd->constraint[k+1]; ddd->block[k]=ddd->block[k+1];
	ddd->matind[k]=ddd->matind[k+1];ddd->nnz[k]=ddd->nnz[k+1];
      }
      ddd->constraint[nonzero-1]=bigint;ddd->nnz[nonzero-1]=0;
      nonzero--;
    }
  }

  ddd->fixedvari=0;ddd->fixedvard=0;
  if (ddd->lpblock>0){ 
    int spot;
    if (ddd->blocksizes[ddd->lpblock]==2){
      i=0;k=0;
      while (ddd->block[i]<=ddd->lpblock && i<nonzero){ i++;} spot=i;
      while (ddd->block[i]==ddd->lpblock+1 && i<nonzero){ i++;k++;}
      if (k==4){
	if (ddd->constraint[spot]==ddd->constraint[spot+1] && 
	    ddd->constraint[spot+2]==ddd->constraint[spot+3] &&
	    ddd->matind[spot]==ddd->matind[spot+2] && 
	    ddd->matind[spot+1]==ddd->matind[spot+3] &&
	    fabs(ddd->nnz[spot+2])==1.0 && fabs(ddd->nnz[spot+3])==1.0 &&
	    fabs(ddd->nnz[spot] + ddd->nnz[spot+1]) <=1e-6   ){
	  ddd->fixedvari=ddd->constraint[spot+2];
	  ddd->fixedvard=ddd->nnz[spot]/ddd->nnz[spot+2];
	  nblocks--;ddd->lpblock=0;
	}
      }
    }
  }
  
  ddd->totalnonzeros=nonzero;
  for (ddd->n=0,i=0;i<nblocks;i++)  ddd->n += ddd->blocksizes[i];
  ddd->m=m;
  ddd->nblocks=nblocks;
  fclose(fp);
  return 0;
}


#undef __FUNCT__  
#define __FUNCT__ "Parseline"
static int Parseline(char thisline[],int *nmat,int *nblk,int *row,
		     int *col,double *value, int *nargs){

  int temp;
  int rtmp,coltmp;

  *nmat=-1;*nblk=-1;rtmp=-1;coltmp=-1;*value=0.0;
  temp=sscanf(thisline,"%d %d %d %d %lg",nmat,nblk,&rtmp,&coltmp,value);
  if (temp==5) *nargs=5;
  else *nargs=0;
  *row=rtmp-1; *col=coltmp-1;

  return(0);
}


static int partition(int list1[], int list2[], int list3[], double list5[], int lstart, int lend){
  int k=lend;
  int pivot1=list1[k],pivot2=list2[k],pivot3=list3[k];
  double pivot5 = list5[k];
  int bottom = lstart-1, top = lend;                                 
  int done = 0;
  int ordered=1;
  while (!done){                    
    
    while (!done) {
      bottom = bottom+1;     
      
      if (bottom == top){ 
	done = 1;                     
	break;
      }
      if ( list1[bottom] > pivot1 || 
	   (list1[bottom] == pivot1 && list2[bottom] > pivot2) ||
	   (list1[bottom] == pivot1 && list2[bottom] == pivot2 && 
	    list3[bottom] > pivot3) ){
	list1[top] = list1[bottom];
	list2[top] = list2[bottom];
	list3[top] = list3[bottom];
	list5[top] = list5[bottom];
	ordered=0;
	break;                     
      }
    }
    while (!done){ 
      top = top-1; 
      
      if (top == bottom){
	done = 1;                 
	break;
      }
      if ( list1[top] < pivot1 || 
	   (list1[top] == pivot1 && list2[top] < pivot2) ||
	   (list1[top] == pivot1 && list2[top] == pivot2 && list3[top] < pivot3)){
	list1[bottom] = list1[top];
	list2[bottom] = list2[top];
	list3[bottom] = list3[top];
	list5[bottom] = list5[top];
	ordered=0;
	break;
      }
    }
  }
  list1[top] = pivot1;
  list2[top] = pivot2;
  list3[top] = pivot3;
  list5[top] = pivot5;
  
  ordered=0;
  
  if (bottom==lend){
    ordered=1;
    for (k=lstart;k<lend-1;k++){
      if ( (list1[k] > list1[k+1]) || 
	   (list1[k] == list1[k+1] && list2[k] > list2[k+1]) ||
	   (list1[k] == list1[k+1] && list2[k] == list2[k+1] && list3[k] > list3[k+1]) ){
	ordered=0;
	break;
      }
    }
  }
  if (ordered && lend-lstart>2){
    top=(lend+lstart)/2;
  }
  return top;
}



static int qusort(int list1[], int list2[], int list3[], double list5[], int lstart, int lend){
  int split;
  if (lstart < lend){
    split = partition(list1, list2, list3, list5,lstart, lend);
    qusort(list1, list2, list3, list5, lstart, split-1);
    qusort(list1, list2, list3, list5, split+1, lend);
  }   else {
    return 0;
  }
  return 0; 
}



#undef __FUNCT__
#define __FUNCT__ "GetMarkers"
static int GetMarkers(int block, int constraint, int blockn[], int constraintn[], 
		       int*m3){
  int i=0;
  while (blockn[i]==block && constraintn[i]==constraint){ i++;}
  *m3=i;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "CountNonzeroMatrices"
static int CountNonzeroMatrices(int block, int blockn[], int constraintn[], int*m3){
  int i=0,cvar=-1,nnzmats=0;
  while (blockn[i]==block){
    if (constraintn[i]>cvar){ cvar=constraintn[i];nnzmats++;}
    i++;
  }
  *m3=nnzmats;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "CheckForConstantMat"
static int CheckForConstantMat(double v[],int nnz, int n){
  int i; double vv;
  if (n<=1){ return 0; }
  if (nnz!=(n*n+n)/2){ return 0; }
  for (vv=v[0],i=1;i<nnz;i++){
    if (v[i]!=vv){ return 0;}
  }
  return 1;
}

static int ComputeY0(DSDP dsdp,DSDPData dddd){
  int i,ii,info,ijnnz=0,spot=0,ddiag=0,diag=0,n=dddd.blocksizes[0],m=dddd.m;
  double aa,bb=0,ddmax=0,dd=0,cnorm=0;
  char sformat=dddd.sformat[0];
  if (dddd.nblocks>1) return 0;
  if (dddd.fixedvari) return 0;

  info=GetMarkers(1,0,dddd.block+spot,dddd.constraint+spot,&ijnnz);	  
  for (i=0;i<ijnnz;i++){if (cnorm<fabs(dddd.nnz[i])) cnorm=fabs(dddd.nnz[i]);}
  spot+=ijnnz;
  
  for (i=1;i<=m;i++,spot+=ijnnz){
    info=GetMarkers(1,i,dddd.block+spot,dddd.constraint+spot,&ijnnz);	  
    ddiag=0;
    if (ijnnz==n){
      dd=dddd.nnz[spot]; ddiag=1;bb=dddd.dobj[i-1];
      for (ii=0;ii<n;ii++){
	if (sformat=='U'){
	  if (dddd.matind[spot+ii] != ii*n+ii || dddd.nnz[spot+i]!=dd){ ddiag=0;}
	} else if (sformat=='P'){
	  if (dddd.matind[spot+ii] != ii*(ii+1)/2+ii || dddd.nnz[spot+i]!=dd){ ddiag=0;}
	}
      }
    }
    if (ddiag){
      info=DSDPSetY0(dsdp,i,-100*n*cnorm*dddd.dobj[i-1]);
      info=DSDPSetR0(dsdp,0); 
      info=DSDPSetZBar(dsdp,100*dd/bb*cnorm); 
      info=DSDPSetPotentialParameter(dsdp,5.0);
    }
  }
  
  if (m==n){
    spot=0; info=GetMarkers(1,0,dddd.block+spot,dddd.constraint+spot,&ijnnz); spot+=ijnnz;
    dd=dddd.nnz[spot]; bb=dddd.dobj[0];
    for (diag=1,i=0;i<n;i++,spot+=ijnnz){
      info=GetMarkers(1,i+1,dddd.block+spot,dddd.constraint+spot,&ijnnz);      
      dd=dddd.nnz[spot];bb=dddd.dobj[i];
      if (ddmax<bb/dd) ddmax=bb/dd;
      if (sformat=='P'){
	if (ijnnz!=1 || bb!=dddd.dobj[i] || dd!=dddd.nnz[spot] || dddd.matind[spot]!= ((i+1)*(i+2))/2-1) { diag=0; }
      } else if (sformat=='U'){
	if (ijnnz!=1 || bb!=dddd.dobj[i] || dd!=dddd.nnz[spot] || dddd.matind[spot]!= ((i)*n)+i) { diag=0; }
      }
    }
    if (diag && cnorm*sqrt(1.0*n)<1e6){
      for (i=0;i<n;i++){info = DSDPSetY0(dsdp,i+1,-10*sqrt(1.0*n)*cnorm);} 
      info=DSDPSetR0(dsdp,0); info=DSDPSetZBar(dsdp,100*ddmax*n*cnorm); info=DSDPSetPotentialParameter(dsdp,5.0); info=DSDPReuseMatrix(dsdp,0);
    }
  }
  if (m==n+1){
    spot=0; info=GetMarkers(1,0,dddd.block+spot,dddd.constraint+spot,&ijnnz); spot+=ijnnz;
    dd=dddd.nnz[spot]; bb=dddd.dobj[0];diag=1;
    info=GetMarkers(1,1,dddd.block+spot,dddd.constraint+spot,&ijnnz); 
    if (CheckForConstantMat(dddd.nnz+spot,ijnnz,n)){aa=dddd.nnz[0]; spot+=ijnnz;ii=2;} else {ii=1;}
    for (i=0;i<n;i++,spot+=ijnnz){
      info=GetMarkers(1,i+ii,dddd.block+spot,dddd.constraint+spot,&ijnnz);
      dd=dddd.nnz[spot];bb=dddd.dobj[i+ii-1];
      if (ddmax<bb/dd) ddmax=bb/dd;
      if (sformat=='U'){
	if (ijnnz!=1 || bb!=dddd.dobj[ii-1] || dd!=dddd.nnz[spot] || dddd.matind[spot]!= ((i)*(n))+i ) { diag=0; }
      } else if (sformat=='P'){
	if (ijnnz!=1 || bb!=dddd.dobj[ii-1] || dd!=dddd.nnz[spot] || dddd.matind[spot]!= ((i+1)*(i+2))/2-1) { diag=0; }
      }
    } 
    if (ii==1 && diag==1){
      info=GetMarkers(1,m,dddd.block+spot,dddd.constraint+spot,&ijnnz);
      if (CheckForConstantMat(dddd.nnz+spot,ijnnz,n)){aa=dddd.nnz[spot];} else {diag=0;}
    }
    if (diag && cnorm*sqrt(1.0*n)<1e5){  
      /*
      if (ii=2){info = DSDPSetY0(dsdp,1,-10000*aa);} else {info = DSDPSetY0(dsdp,m,-10000*aa);}
      for (i=0;i<n;i++){info = DSDPSetY0(dsdp,i+ii,-100*sqrt(1.0*n)*cnorm);} 
      */
      /* info=DSDPSetR0(dsdp,cnorm); info=DSDPSetZBar(dsdp,n*n*n*ddmax*cnorm); */ info=DSDPSetPotentialParameter(dsdp,5.0); info=DSDPReuseMatrix(dsdp,0);
    }
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TCheckArgs0"
static int TCheckArgs0(DSDP dsdp,SDPCone sdpcone, int m,int nargs,char *runargs[]){

  int kk,info,iloginfo=0;
  for (kk=1; kk<nargs-1; kk++){
    if (0){
    } else if (strncmp(runargs[kk],"-dloginfo",8)==0){
      iloginfo=atoi(runargs[kk+1]);
    }
  }
  info=DSDPLogInfoAllow(iloginfo,0);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TCheckArgs"
static int TCheckArgs(DSDP dsdp,SDPCone sdpcone, int m,int nargs,char *runargs[]){

  int i,kk, info;
  int printlevel=10;
  double *yy0;

  for (kk=1; kk<nargs-1; kk++){
    if (strncmp(runargs[kk],"-y0",7)==0){
      DSDPCALLOC2(&yy0,double,m,&info);
      for (i=0;i<m;i++) yy0[i]=0.0;
      info = ReadInitialPoint(runargs[kk+1],m,yy0);
      for (i=0;i<m;i++){info = DSDPSetY0(dsdp,i+1,yy0[i]);}
      DSDPFREE(&yy0,&info);
    } else if (strncmp(runargs[kk],"-params",6)==0){
      info=DSDPReadOptions(dsdp,runargs[kk+1]);
    } else if (strncmp(runargs[kk],"-print",6)==0){
      printlevel=atoi(runargs[kk+1]);
    } 
  }

  info=DSDPSetOptions(dsdp,runargs,nargs);
  /*  if (lpcone){info=LPConeScaleBarrier(lpcone,lpb); */
  if (rank==0){info=DSDPSetStandardMonitor(dsdp,printlevel);}
  if (rank==0){info=DSDPSetFileMonitor(dsdp,printlevel);}
  return(0);
}

#undef __FUNCT__
#define __FUNCT__ "ReadInitialPoint"
static int ReadInitialPoint(char* filename, int m, double yy0[])
{
  FILE   *fp;
  int i,count;

  fp = fopen(filename,"r");
  if(!fp)
  { printf("\n Cannot open file containing initial dual solution %s",filename);
    return 0;
  }

  for(count=0,i=0;(i < m) &&(!feof(fp));i++){
    if(fscanf(fp,"%lf",&yy0[i] ) ){ count++; }
  } 

  if (count < m){
    printf("Warning: reading initial y vector: \n");
    printf("   Expecting %d entries but read only %d entries \n",m,count);
  }
  fclose(fp);
  return 0;
} /* ReadUserY */
