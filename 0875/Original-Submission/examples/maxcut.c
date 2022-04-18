#include "dsdp5.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*! \file maxcut.c
  \brief Most Basic Example: read graph from file, formulate the SDP relaxation of maximum cut
problem, solve using DSDP, and apply randomized algorithm to generate
approximate solutions.
 */

char help[]="\
DSDP Usage: maxcut <graph filename> \n\
                   -gaptol <relative duality gap: default is 0.001> \n\
                  -maxit <maximum iterates: default is 200> \n";

static int ReadGraph(char*,int *, int *,int**, int **, double **); 
static int TCheckArgs(DSDP,int,char **);
static int ParseGraphline(char*,int*,int*,double*,int*);
int MaxCutRandomized(SDPCone,int);
int MaxCut(int,int, int[], int[], double[]);


int main(int argc,char *argv[]){
  int info;
  int *node1,*node2,nedges,nnodes;
  double *weight;

  if (argc<2){ printf("%s",help); return(1); }

  info = ReadGraph(argv[1],&nnodes,&nedges,&node1,&node2,&weight);
  if (info){ printf("Problem reading file\n"); return 1; }

  MaxCut(nnodes,nedges,node1,node2,weight);

  free(node1);free(node2);free(weight);
  return 0;
}

/*! 
\fn int MaxCut(int nnodes,int nedged, int node1[], int node2[], double weight[]);
\brief Formulate and solve the SDP relaxation of the Maximum Cut problem
\ingroup Examples
\param nnodes number of nodes in graph
\param nedges number of edges in graph
\param node1 first node of each edge
\param node2 second node of each edge
\param weight weight of each edge
\note This routine is an example! It is not part of the solver library.
*/
int MaxCut(int nnodes,int nedges, int node1[], int node2[], double weight[]){
  
  int i,info;
  int *indd,*iptr;
  double *yy,*val,*diag,tval=0;
  DSDPTerminationReason reason;
  SDPCone sdpcone;
  DSDP dsdp;
  
  info = DSDPCreate(nnodes,&dsdp); 
  info = DSDPCreateSDPCone(dsdp,1,&sdpcone);

  if (info){ printf("Out of memory\n"); return 1; }

  info = SDPConeSetBlockSize(sdpcone,0,nnodes);


  /* Formulate the problem from the data */
  /* 
     Diagonal elements equal 1.0 
     Create Constraint matrix A_i for i=1, ..., nnodes.
     that has a single nonzero element.
  */
  diag=(double*)malloc(nnodes*sizeof(double));
  iptr=(int*)malloc(nnodes*sizeof(int));
  for (i=0;i<nnodes;i++){
    iptr[i]=i*(i+1)/2+i; 
    diag[i]=1.0;
  }

  for (i=0;i<nnodes;i++){
    info = DSDPSetDualObjective(dsdp,i+1,1.0);
    info = SDPConeSetASparseVecMat(sdpcone,0,i+1,nnodes,1.0,0,iptr+i,diag+i,1);
    if (0==1){
      printf("Matrix: %d\n",i+1);
      info = SDPConeViewDataMatrix(sdpcone,0,i+1);
    }
  }

  /* C matrix is the Laplacian of the adjacency matrix */
  /* Also compute a feasible initial point y such that S>=0 */ 
  yy=(double*)malloc(nnodes*sizeof(double));
  for (i=0;i<nnodes;i++){yy[i]=0.0;}
  indd=(int*)malloc((nnodes+nedges)*sizeof(int));
  val=(double*)malloc((nnodes+nedges)*sizeof(double));
  for (i=0;i<nnodes+nedges;i++){indd[i]=0;}
  for (i=0;i<nnodes;i++){indd[nedges+i]=i*(i+1)/2+i;}
  for (i=0;i<nnodes+nedges;i++){val[i]=0;}
  for (i=0;i<nedges;i++){
    indd[i]=(node1[i])*(node1[i]+1)/2 + node2[i];
    tval+=fabs(weight[i]);
    val[i]=weight[i]/4;
    val[nedges+node1[i]]-=weight[i]/4;
    val[nedges+node2[i]]-=weight[i]/4;
    yy[node1[i]]-= fabs(weight[i]/2);
    yy[node2[i]]-= fabs(weight[i]/2);
  }

  if (0){
    info = SDPConeSetASparseVecMat(sdpcone,0,0,nnodes,1.0,0,indd,val,nedges+nnodes);
  } else { /* Equivalent */
    info = SDPConeSetASparseVecMat(sdpcone,0,0,nnodes,1.0,0,indd,val,nedges);
    info = SDPConeAddASparseVecMat(sdpcone,0,0,nnodes,1.0,0,indd+nedges,val+nedges,nnodes);
  } 
  if (0==1){ info = SDPConeViewDataMatrix(sdpcone,0,0);}
  
  /* Initial Point */
  info = DSDPSetR0(dsdp,0.0);
  info = DSDPSetZBar(dsdp,10*tval+1.0);
  for (i=0;i<nnodes; i++){ 
    info = DSDPSetY0(dsdp,i+1,10*yy[i]);
  }
  if (info) return info;

  /* Get read to go */
  info=DSDPSetGapTolerance(dsdp,0.001);
  info=DSDPSetPotentialParameter(dsdp,5);
  info=DSDPReuseMatrix(dsdp,0);
  info=DSDPSetPNormTolerance(dsdp,1.0);
  /*
  info = TCheckArgs(dsdp,argc,argv);
  */

  if (info){ printf("Out of memory\n"); return 1; }
  info = DSDPSetStandardMonitor(dsdp,1);

  info = DSDPSetup(dsdp);
  if (info){ printf("Out of memory\n"); return 1; }

  info = DSDPSolve(dsdp);
  if (info){ printf("Numerical error\n"); return 1; }
  info = DSDPStopReason(dsdp,&reason); 

  if (reason!=DSDP_INFEASIBLE_START){ /* Randomized solution strategy */
    info=MaxCutRandomized(sdpcone,nnodes);
    if (0==1){ /* Look at the solution */
      int n; double *xx,*y=diag;
      info=DSDPGetY(dsdp,y,nnodes);
      info=DSDPComputeX(dsdp);DSDPCHKERR(info);
      info=SDPConeGetXArray(sdpcone,0,&xx,&n);
    }
  }
  info = DSDPDestroy(dsdp);

  free(iptr);
  free(yy);
  free(indd);
  free(val);
  free(diag);
  
  return 0;
} /* main */



/*!
int MaxCutRandomized(SDPCone sdpcone,int nnodes);
\brief Apply the Goemens and Williamson randomized cut algorithm to the SDP relaxation of the max-cut problem
\param sdpcone the SDP cone
\param nnodes number of nodes in the graph
\note This routine is an example! It is not part of the solver library.
\ingroup Examples
\sa  MaxCut()
*/
int MaxCutRandomized(SDPCone sdpcone,int nnodes){
  int i,j,derror,info;
  double dd,scal=RAND_MAX,*vv,*tt,*cc,ymin=0;

  vv=(double*)malloc(nnodes*sizeof(double));
  tt=(double*)malloc(nnodes*sizeof(double));
  cc=(double*)malloc((nnodes+2)*sizeof(double));
  info=SDPConeComputeXV(sdpcone,0,&derror);
  for (i=0;i<nnodes;i++){
      for (j=0;j<nnodes;j++){dd = (( rand())/scal - .5); vv[j] = tan(3.1415926*dd);}
      info=SDPConeXVMultiply(sdpcone,0,vv,tt,nnodes);
      for (j=0;j<nnodes;j++){if (tt[j]<0) tt[j]=-1; else tt[j]=1;}
      for (j=0;j<nnodes+2;j++){cc[j]=0;}
      info=SDPConeAddXVAV(sdpcone,0,tt,nnodes,cc,nnodes+2);
      if (cc[0]<ymin) ymin=cc[0];
  }
  printf("Best integer solution: %4.2f\n",ymin);
  free(vv); free(tt); free(cc);

  return(0);
}

static int TCheckArgs(DSDP dsdp,int nargs,char *runargs[]){

  int kk, info;
  
  info=DSDPSetOptions(dsdp,runargs,nargs);
  for (kk=1; kk<nargs-1; kk++){
    if (strncmp(runargs[kk],"-dloginfo",8)==0){
      info=DSDPLogInfoAllow(DSDP_TRUE,0);
    } else if (strncmp(runargs[kk],"-params",7)==0){
      info=DSDPReadOptions(dsdp,runargs[kk+1]);
    } else if (strncmp(runargs[kk],"-help",7)==0){
      printf("%s\n",help);
    } 
  }

  return 0;
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

  return(0);
}


#undef __FUNCT__  
#define __FUNCT__ "ReadGraph"
int ReadGraph(char* filename,int *nnodes, int *nedges,
	      int**n1, int ** n2, double **wght){

  FILE*fp;
  char thisline[BUFFERSIZ]="*";
  int i,k=0,line=0,nodes,edges,gotone=3;
  int *node1,*node2;
  double *weight;
  int info,row,col;
  double value;

  fp=fopen(filename,"r");
  if (!fp){printf("Cannot open file %s !",filename);return(1);}

  while(!feof(fp) && (thisline[0] == '*' || thisline[0] == '"') ){
    fgets(thisline,BUFFERSIZ,fp); line++;
  }

  if (sscanf(thisline,"%d %d",&nodes, &edges)!=2){
    printf("First line must contain the number of nodes and number of edges\n");
    return 1;
  }

  node1=(int*)malloc(edges*sizeof(int));
  node2=(int*)malloc(edges*sizeof(int));
  weight=(double*)malloc(edges*sizeof(double));

  for (i=0; i<edges; i++){ node1[i]=0;node2[i]=0;weight[i]=0.0;}
  
  while(!feof(fp) && gotone){
    thisline[0]='\0';
    fgets(thisline,BUFFERSIZ,fp); line++;
    info = ParseGraphline(thisline,&row,&col,&value,&gotone);
    if (gotone && value!=0.0 && k<edges && 
	col < nodes && row < nodes && col >= 0 && row >= 0){
      if (row<col){info=row;row=col;col=info;}
      if (row == col){}
      else {
	node1[k]=row;        node2[k]=col;
	weight[k]=value;     k++;
      }
    } else if (gotone &&  k>=edges) {
      printf("To many edges in file.\nLine %d, %s\n",line,thisline);
      return 1;
    } else if (gotone&&(col >= nodes || row >= nodes || col < 0 || row < 0)){
      printf("Invalid node number.\nLine %d, %s\n",line,thisline);
      return 1;
    }
  }
  *nnodes=nodes; *nedges=edges;
  *n1=node1; *n2=node2; *wght=weight;
  return 0;
}
