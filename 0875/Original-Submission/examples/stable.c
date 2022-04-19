#include "dsdp5.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

/*! \file stable.c
  \brief Read graph from file, formulate the Maximum Stable Set problem, and solve using DSDP.
 */

char help[]="\n\
Compute the stable set for a graph. \n\n\
DSDP Usage: stable <graph filename> \n";

typedef struct {
  double v[3];
  int indd[3];
} EdgeMat;

static int ReadGraphFromFile(char*,int*, int*, EdgeMat*[]); 
int SetStableSetData(DSDP, SDPCone, int, int, EdgeMat[]);
int StableSet(int argc,char *argv[]);
int StableRandomized(SDPCone,int, int, EdgeMat[]);
#define CHKERR(a)  { if (a){printf("DSDP Numerical Error or Memory Problem"); return 2;} }

int main(int argc,char *argv[]){
  int info;
  info=StableSet(argc,argv);
  return 0;
}

/*!
\fn int StableSet(int argc,char *argv[]);
\param argc number of command line arguments
\param argv command line arguments
\brief Formulate and solve the maximum Stable Set problem.
\ingroup Examples
\sa SetStableSetData()
 */
int StableSet(int argc,char *argv[]){

  int info,kk,nedges,nnodes;
  EdgeMat *Edges;
  DSDP dsdp;
  SDPCone  sdpcone;

  if (argc<2){ printf("%s",help); return(1); }

  info = ReadGraphFromFile(argv[1],&nnodes,&nedges,&Edges);
  if (info){ printf("Problem reading file\n"); return 1; }

  info = DSDPCreate(nedges+nnodes+1,&dsdp);CHKERR(info); 
  info = DSDPCreateSDPCone(dsdp,1,&sdpcone);CHKERR(info); 
  info = SDPConeSetBlockSize(sdpcone,0,nnodes+1);CHKERR(info); 
  info = SDPConeUsePackedFormat(sdpcone,0);CHKERR(info); 
  info = SDPConeSetSparsity(sdpcone,0,nedges+nnodes+1);CHKERR(info); 
  info = SetStableSetData(dsdp, sdpcone, nnodes, nedges, Edges);
  if (info){ printf("Problem setting data\n"); return 1; }

  info = DSDPSetGapTolerance(dsdp,0.0001);CHKERR(info); 
  info = DSDPSetZBar(dsdp,1e10*nnodes+1);CHKERR(info); 
  info = DSDPReuseMatrix(dsdp,10);CHKERR(info); 

  for (kk=1; kk<argc-1; kk++){
    if (strncmp(argv[kk],"-dloginfo",8)==0){
      info=DSDPLogInfoAllow(DSDP_TRUE,0);CHKERR(info); 
    } else if (strncmp(argv[kk],"-params",7)==0){
      info=DSDPReadOptions(dsdp,argv[kk+1]);CHKERR(info); 
    } else if (strncmp(argv[kk],"-help",5)==0){
      printf("%s\n",help);
    } 
  }
  info=DSDPSetOptions(dsdp,argv,argc);CHKERR(info); 

  info = DSDPSetStandardMonitor(dsdp,1);CHKERR(info); 

  info = DSDPSetup(dsdp);CHKERR(info); 
  if (0==1){info=SDPConeCheckData(sdpcone);}
  
  info=DSDPSolve(dsdp); CHKERR(info);
  info=StableRandomized(sdpcone,nnodes,nedges,Edges);

  info=DSDPComputeX(dsdp);CHKERR(info);
  
  if (0==1){ /* Look at the solution */
    int n; double *xx;
    info=SDPConeGetXArray(sdpcone,0,&xx,&n);CHKERR(info); 
  }

  info = DSDPDestroy(dsdp);CHKERR(info); 
  free(Edges);
  
  return 0;
} /* main */

/*!
\fn int SetStableSetData(DSDP dsdp, SDPCone sdpcone, int nodes, int edges, EdgeMat Edge[]);
\param dsdp the solver
\param sdpcone the semidefinite cone
\param nodes number of nodes in graph
\param edges number of edges in graph
\param Edge edges in graph
\brief Given a graph, formulate maximum Stable Set problem and place data into solver.
\ingroup Examples
\sa StableSet
 */
int SetStableSetData(DSDP dsdp, SDPCone sdpcone, int nodes, int edges, EdgeMat Edge[]){

  int i,ii,info,nnodes=nodes+1;
  int *iptr,*iptr2;
  double *diag;

  /* The c matrix has all elements equal to 1.0 */
  diag=(double*)malloc(nnodes*sizeof(double));
  iptr=(int*)malloc(nnodes*sizeof(int));
  iptr2=(int*)malloc(nnodes*sizeof(int));

  ii=nodes*(nodes+1)/2;
  for (i=0;i<nnodes;i++){
    diag[i]=1.0;
    iptr[i]=i*(i+1)/2+i; 
    iptr2[i]=i;
  }
  info = SDPConeSetASparseVecMat(sdpcone,0,0,nnodes,-0.50,0,iptr,diag,nodes);CHKERR(info); 
  info = SDPConeAddASparseVecMat(sdpcone,0,0,nnodes,-0.25,-ii,iptr2,diag,nodes);CHKERR(info); 
  if (0){info=SDPConeViewDataMatrix(sdpcone,0,0);}
  /* Diagonal elements must equal 1 */
  for (i=0;i<nnodes;i++){
    info = DSDPSetDualObjective(dsdp,i+1,1.0);CHKERR(info); 
    info = DSDPSetY0(dsdp,i+1,0.0);CHKERR(info); 
    info = SDPConeSetASparseVecMat(sdpcone,0,i+1,nnodes,1.0,0,iptr+i,diag+i,1);CHKERR(info); 
  }

  /* 
     For each edge connecting nodes i and j, 
     X(i,i)+X(j,j)+X(i,j)+X(j,i)+X(i,n)+X(n,i)+X(j,n)+X(n,j)+X(n,n) = 1
     where nodes i,j numbered 0 ... n-1.
  */
  for (i=0; i<edges; i++){
    info = SDPConeAddARankOneMat(sdpcone,0,i+nnodes+1,nnodes,1.0,0,Edge[i].indd,Edge[i].v,3);
    if (0==1){info = SDPConeViewDataMatrix(sdpcone,0,i+nnodes+1);CHKERR(info);}
    info = DSDPSetDualObjective(dsdp,i+nnodes+1,1.0);CHKERR(info); 
    info = DSDPSetY0(dsdp,i+nnodes+1,0.0);CHKERR(info); 
  }

  /* The initial point y is feasible and near the central path */
  /*
  info = DSDPSetR0(dsdp,0);
  */  
  return(0);
}

/*!
int StableRandomized(SDPCone sdpcone,int nodes, int edges, EdgeMat Edge[]);
\brief Apply a randomized procedure to find feasible stable sets.
\param sdpcone the SDP cone
\param nodes number of nodes in the graph
\param edges number of edges in the graph
\param Edge Array of edges
\note This routine is an example! It is not part of the solver library.
\ingroup Examples
\sa  MaxCutRandomized()
*/
int StableRandomized(SDPCone sdpcone,int nodes, int edges, EdgeMat Edge[]){
  int i,j,derror,info,nnodes=nodes+1;
  double dd,scal=RAND_MAX,*vv,*tt,*cc,ymin=0;
  int e0,e1,e2;

  vv=(double*)malloc(nnodes*sizeof(double));
  tt=(double*)malloc(nnodes*sizeof(double));
  cc=(double*)malloc((edges+nnodes+2)*sizeof(double));
  info=SDPConeComputeXV(sdpcone,0,&derror);
  for (i=0;i<nnodes;i++){
      for (j=0;j<nnodes;j++){dd = (( rand())/scal - .5); vv[j] = tan(3.1415926*dd);}
      info=SDPConeXVMultiply(sdpcone,0,vv,tt,nnodes);
      for (j=0; j<edges; j++){
	e0=Edge[j].indd[0];e1=Edge[j].indd[1];e2=Edge[j].indd[2];
	if (tt[e0] * tt[e2] > 0 && tt[e1]*tt[e2] >0){
	  if ( fabs(tt[e0]-tt[e2]) > fabs(tt[e1]-tt[e2]) ){
	    tt[e0]*=-1;
	  } else {
	    tt[e1]*=-1;
	  }
	} 
      }
      for (j=0;j<nnodes;j++){if (tt[j]<0) tt[j]=-1; else tt[j]=1;}
      for (j=0;j<edges+nodes+1;j++){cc[j]=0;}
      info=SDPConeAddXVAV(sdpcone,0,tt,nnodes,cc,edges+nnodes+2);
      if (cc[0]<ymin) ymin=cc[0];
  }
  printf("Stable Set Size: %4.0f\n",-ymin);
  free(vv); free(tt); free(cc);

  return(0);
}


#define BUFFERSIZ 100

#undef __FUNCT__  
#define __FUNCT__ "ParseGraphline"
static int ParseGraphline(char * thisline,int *row,int *col,double *value, 
			  int *gotem){

  int temp;
  int rtmp,coltmp;

  rtmp=-1, coltmp=-1, *value=0.0;
  temp=sscanf(thisline,"%d %d %lf",&rtmp,&coltmp,value);
  if (temp==3 && rtmp>0 && coltmp>0) *gotem=3;
  else if (temp==2 && rtmp>0 && coltmp>0){ *value = 1.0; *gotem=3;}
  else *gotem=0;
  *row=rtmp-1; *col=coltmp-1;
  if (*gotem && (*col < 0 || *row < 0)){
    printf("Node Number must be positive.\n, %s\n",thisline);
  }
  return(0);
}

#undef __FUNCT__  
#define __FUNCT__ "ReadGraphFromFile"
static int ReadGraphFromFile(char* filename,int *nnodes, int *nedges, EdgeMat**EE){

  FILE*fp;
  char thisline[BUFFERSIZ]="*";
  int i,k=0,line=0,nodes,edges,gotone=3;
  int info,row,col;
  double value;
  EdgeMat *E;

  fp=fopen(filename,"r");
  if (!fp){ printf("Cannot open file %s !",filename); return(1); }

  while(!feof(fp) && (thisline[0] == '*' || thisline[0] == '"') ){
    fgets(thisline,BUFFERSIZ,fp); line++; }

  if (sscanf(thisline,"%d %d",&nodes, &edges)!=2){
    printf("First line must contain the number of nodes and number of edges\n");
    return 1;
  }

  E=(EdgeMat*)malloc(edges*sizeof(EdgeMat)); *EE=E;
  for (i=0; i<edges; i++){ 
    E[i].v[0]=0; E[i].v[1]=0; E[i].v[2]=0; 
    E[i].indd[0]=0; E[i].indd[1]=0; E[i].indd[2]=0; 
  }

  while(!feof(fp) && gotone){
    thisline[0]='\0';
    fgets(thisline,BUFFERSIZ,fp); line++;
    info = ParseGraphline(thisline,&row,&col,&value,&gotone);
    if (gotone && k<edges && 
	col < nodes && row < nodes && col >= 0 && row >= 0){
      if (row > col){i=row;row=col;col=i;}
      if (row == col){}
      else {
	E[k].indd[0]=row; E[k].indd[1]=col; E[k].indd[2]=nodes; 
	E[k].v[0]=1.0;    E[k].v[1]=1.0;    E[k].v[2]=1.0; 
	k++;
      }
    } else if (gotone &&  k>=edges) {
      printf("To many edges in file.\nLine %d, %s\n",line,thisline);
      return 1;
    } else if (gotone&&(col >= nodes || row >= nodes || col < 0 || row < 0)){
      printf("Invalid node number.\nLine %d, %s\n",line,thisline);
      return 1;
    }
  }

  *nnodes=nodes; *nedges=k;
						 
  return 0;
}

