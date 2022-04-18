#include "dsdp5.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*! \file color.c
  \brief Second Basic Example: Read graph from file, formulate the SDP relaxation of k-coloring
problem, solve using DSDP, and apply randomized algorithm to generate
approximate solutions.
 */
char help2[]="\nA positive semidefinite relaxation of the\n\
graph k coloring problem can be rewritten as\n\n\
   Find       X>=0 \n\
   such that  X_ij <= 1 - 1/(k-1) for all edges (i,j).\n\
";

char help[]="DSDP Usage: color <graph filename> ";

static int ReadGraph(char*,int *, int *,int**, int **, double **); 
static int ParseGraphline(char*,int*,int*,double*,int*);
static int RandomizedColor(DSDP, SDPCone, int, int[], int[], int);
int MinColoring(int argc,char *argv[]);


int main(int argc,char *argv[]){
  int info;
  info=MinColoring(argc,argv);
  return 0;
}

/*!
\fn int MinColoring(int argc,char *argv[]);
\param argc number of command line arguments
\param argv command line arguments
\ingroup Examples
\brief SDP relaxation of k-coloring problem
 */
int MinColoring(int argc,char *argv[]){

  int i,kk,vari,info;
  int *node1,*node2,nedges,nnodes;
  int *iptr1,*iptr2;
  double *weight,*yy1,*yy2,bb; 
  DSDPTerminationReason reason;
  SDPCone sdpcone;
  BCone bcone;
  DSDP dsdp;

  if (argc<2){ printf("%s\n%s",help2,help); return(1); }

  info = ReadGraph(argv[1],&nnodes,&nedges,&node1,&node2,&weight);
  if (info){ printf("Problem reading file\n"); return 1; }

  info = DSDPCreate(nnodes+nedges,&dsdp); 

  info = DSDPCreateSDPCone(dsdp,1,&sdpcone);
  info = SDPConeSetBlockSize(sdpcone,0,nnodes);
  info = SDPConeSetSparsity(sdpcone,0,nnodes+nedges+1);

  info = DSDPCreateBCone(dsdp,&bcone);

  if (info){ printf("Out of memory\n"); return 1; }


  /* Formulate the problem from the data */  
  /* Create data matrices */

  /* Diagonal elements of X(i,i) must equal 1.0 */
  iptr1=(int*)malloc(nnodes*sizeof(int));
  yy1=(double*)malloc(nnodes*sizeof(double));
  for (i=0;i<nnodes;i++){
    iptr1[i]=(i+2)*(i+1)/2-1; 
    yy1[i]=1.0;
  }
  for (i=0;i<nnodes;i++){
    info=SDPConeSetASparseVecMat(sdpcone,0,i+1,nnodes,1.0,0,iptr1+i,yy1+i,1);
    if (info) printf("ERROR 1: %d \n",i);
    info=DSDPSetDualObjective(dsdp,i+1,1.0);
  }

  /* For each nonzero element (i,j) of the matrix, X(i,j) must be less than 1 - 1/nnodes */
  bb=2-2.0/nnodes;
  iptr2=(int*)malloc(nedges*sizeof(int));
  yy2=(double*)malloc(nedges*sizeof(double));
  for (i=0;i<nedges;i++){
    iptr2[i]=(node1[i])*(node1[i]+1)/2+node2[i]; 
    yy2[i]=1.0;
  }
  info = BConeAllocateBounds(bcone,nedges);
  for (i=0; i<nedges; i++){
    vari=nnodes+i+1;
    info = SDPConeSetSparseVecMat(sdpcone,0,vari,nnodes,0,iptr2+i,yy2+i,1);
    if (info) printf("ERROR 2: %d %d \n",i,vari);
    info = BConeSetPSlackVariable(bcone,vari);
    if (info) printf("ERROR 3: %d %d \n",i,vari);
    info = DSDPSetDualObjective(dsdp,vari,bb);
  }

  
  /* Get read to go */
  info=DSDPSetPotentialParameter(dsdp,5);
  
  for (kk=1; kk<argc-1; kk++){
    if (strncmp(argv[kk],"-dloginfo",8)==0){
      info=DSDPLogInfoAllow(DSDP_TRUE,0);
    } else if (strncmp(argv[kk],"-params",7)==0){
      info=DSDPReadOptions(dsdp,argv[kk+1]);
    } else if (strncmp(argv[kk],"-help",7)==0){
      printf("%s\n",help);
    } 
  }
  info=DSDPSetOptions(dsdp,argv,argc);

  if (info){ printf("Out of memory\n"); return 1; }
  info = DSDPSetStandardMonitor(dsdp,1);

  info = DSDPSetup(dsdp);
  if (info){ printf("Out of memory\n"); return 1; }

  info = DSDPSolve(dsdp);
  if (info){ printf("Numerical error\n"); return 1; }
  info = DSDPStopReason(dsdp,&reason); 

  if (reason!=DSDP_INFEASIBLE_START){ /* Randomized solution strategy */
    info=RandomizedColor(dsdp, sdpcone, nnodes, node1, node2, nedges);
  }

  info = DSDPDestroy(dsdp);
  
  free(node1);free(node2);free(weight);
  free(iptr1);
  free(yy1);
  free(iptr2);
  free(yy2);
  
  return 0;
} 

static int GetXRow(double xmat[],double xrow[],int row,int n){
  int i,i1=row*(row+1)/2;
  for (i=0;i<row;i++){xrow[i]=xmat[i1+i];}
  for (i=row;i<n;i++){xrow[i]=xmat[i1+row];i1+=i+1;}
  return 0;
}

typedef struct {
  int    index;double val;
} orderVec;

static int cut_comp(const void *e1,const void *e2){ /* Used in qsort routine */
   double d1=((orderVec*)e1)->val, d2=((orderVec*)e2)->val;
   if (d1<d2) return (1);
   else if (d1>d2) return (-1);
   return(0);
}

static int RemoveNode(int node, int node1[], int node2[], int *nedges){
  int i,nnedges=*nedges;
  for (i=0;i<nnedges;i++){
    if (node1[i]==node || node2[i]==node){
      node1[i]=node1[nnedges-1];
      node2[i]=node2[nnedges-1];
      nnedges--;
      if (i < nnedges) i--;
    }
  }
  *nedges=nnedges;
  return 0;
}

static int Connected(int n1, int n2, int node1[], int node2[], int nedges){
  int i;
  if (n1==n2) return 1;
  for (i=0;i<nedges;i++){
    if (node1[i]==n1 && node2[i]==n2){ return 1;}
    if (node1[i]==n2 && node2[i]==n1){ return 1;}
  }
  return 0;
}

static int HighDegree(int node1[], int node2[], int nedges, int iwork[], int nnodes){
  int i,nmax=0,maxdegree=-1;
  for (i=0;i<nnodes;i++) iwork[i]=0;
  for (i=0;i<nedges;i++){
    iwork[node1[i]]++; iwork[node2[i]]++;
  }
  for (i=0;i<nnodes;i++){ if (iwork[i]>maxdegree){nmax=i; maxdegree=iwork[i];}  }
  return nmax;
}

static int First(int coloring[], int nnodes){
  int i,nmax=nnodes;
  for (i=0;i<nnodes;i++){
    if (coloring[i]==0){ 
      nmax=i; return nmax;
    }
  }
  return -1;
}

static int RandomizedColor(DSDP dsdp, SDPCone sdpcone, int nnodes, int node1[], int node2[], int nedges){
  int i,j,nodek,nn,info,flag,coloring=0,maxvertex;
  int *degree,*color,*cgroup,ngsize,uncolored=nnodes;
  int tnedges=nedges;
  double *xrow,*xptr;
  orderVec *vorder;

  xrow=(double*)malloc(nnodes*sizeof(double));
  color=(int*)malloc(nnodes*sizeof(int));
  cgroup=(int*)malloc(nnodes*sizeof(int));
  degree=(int*)malloc(nnodes*sizeof(int));
  vorder=(orderVec*)malloc(nnodes*sizeof(orderVec));

  for (i=0;i<nnodes;i++){ color[i]=0;}
  info=DSDPComputeX(dsdp);
  info=SDPConeGetXArray(sdpcone,0,&xptr,&nn);

  while (uncolored>0){

    coloring++;
    
    maxvertex=First(color,nnodes);
    maxvertex=HighDegree(node1,node2,tnedges,degree,nnodes);

    cgroup[0]=maxvertex;ngsize=1;

    info=GetXRow(xptr,xrow,maxvertex,nnodes);

    for (i=0;i<nnodes;i++){vorder[i].index=i; vorder[i].val = xrow[i];}
    qsort( (void*)vorder, nnodes, sizeof(orderVec), cut_comp);

    for (i=0;i<nnodes;i++){
      nodek=vorder[i].index;
      if (color[nodek]==0){
	for (flag=0,j=0;j<ngsize;j++){
	  if (Connected(nodek,cgroup[j],node1,node2,tnedges) ){flag=1;}
	}
	if (flag==0){ cgroup[ngsize]=nodek; ngsize++; }
      }
    }
    for (i=0;i<ngsize;i++){
      color[cgroup[i]]=coloring; uncolored--;
      RemoveNode(cgroup[i],node1,node2,&tnedges);
    }
  }
  printf("\nCOLORS: %d\n",coloring);
  free(xrow);
  free(color);
  free(cgroup);
  free(degree);
  free(vorder);
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
