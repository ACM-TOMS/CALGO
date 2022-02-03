#include "dsdp5.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

/*! \file theta.c
  \brief Read graph complement from file, formulate the Lovasz theta problem, and solve using DSDP.
 */

char help[]="\n\
Compute the Lovasz theta number for a graph.  This number is an upper bound for \n\
the maximum clique of a graph, a lower bound for the mimimal graph coloring, and serves \n\
as a bound for several other combitorial graph problems.  The number is the solution to \n\
a semidfinite program. \n\n\
The file should be the complement of the graph \n\
This file also demonstrates the use of customized data matrices in DSDP.\n\n\
DSDP Usage: theta <graph complement filename> \n";

typedef struct {
  int v1,v2;
} EdgeMat;

#include "../src/sdp/dsdpdatamat_impl.h"
extern int SDPConeAddDataMatrix(SDPCone,int, int, int, char, struct DSDPDataMat_Ops*, void*);


/* These variables and routines are for customized SDP Data Matrices */
static struct  DSDPDataMat_Ops edgematop;
static int EdgeMatOperationsInitialize(struct  DSDPDataMat_Ops*);

static struct  DSDPDataMat_Ops onematops;
static int OneMatOpsInitialize(struct  DSDPDataMat_Ops*);

static struct  DSDPDataMat_Ops identitymatops;
static int IdentitymatOperationsInitialize(struct  DSDPDataMat_Ops*);

static int ReadGraphFromFile(char*,int*, int*, EdgeMat*[]); 
int SetThetaData(DSDP, SDPCone, int, int, EdgeMat[]);
int LovaszTheta(int argc,char *argv[]);

int main(int argc,char *argv[]){
  int info;
  info=LovaszTheta(argc,argv);
  return 0;
}

/*!
\fn int LovaszTheta(int argc,char *argv[]);
\param argc number of command line arguments
\param argv command line arguments
\brief Formulate and solve the Lovasz theta problem.
\ingroup Examples
\sa SetThetaData()
 */
int LovaszTheta(int argc,char *argv[]){

  int info,kk,nedges,nnodes;
  EdgeMat *Edges;
  DSDP dsdp;
  SDPCone  sdpcone;

  if (argc<2){ printf("%s",help); return(1); }

  info = ReadGraphFromFile(argv[1],&nnodes,&nedges,&Edges);
  if (info){ printf("Problem reading file\n"); return 1; }

  info = DSDPCreate(nedges+1,&dsdp); 
  info = DSDPCreateSDPCone(dsdp,1,&sdpcone);
  info = SDPConeSetBlockSize(sdpcone,0,nnodes);
  info = SDPConeUsePackedFormat(sdpcone,0);
  info = SDPConeSetSparsity(sdpcone,0,nedges+1);
  if (info){ printf("Out of memory\n"); return 1; }
  info = SetThetaData(dsdp, sdpcone, nnodes, nedges, Edges);
  if (info){ printf("Out of memory\n"); return 1; }

  info = DSDPSetGapTolerance(dsdp,0.001);
  info = DSDPSetZBar(dsdp,nnodes+1);
  info = DSDPReuseMatrix(dsdp,10);

  for (kk=1; kk<argc-1; kk++){
    if (strncmp(argv[kk],"-dloginfo",8)==0){
      info=DSDPLogInfoAllow(DSDP_TRUE,0);
    } else if (strncmp(argv[kk],"-params",7)==0){
      info=DSDPReadOptions(dsdp,argv[kk+1]);
    } else if (strncmp(argv[kk],"-help",5)==0){
      printf("%s\n",help);
    } 
  }
  info=DSDPSetOptions(dsdp,argv,argc);

  if (info){ printf("Out of memory\n"); return 1; }
  info = DSDPSetStandardMonitor(dsdp,1);
  if (info){ printf("Monitor Problem \n"); return 1; }

  info = DSDPSetup(dsdp);
  if (info){ printf("Out of memory\n"); return 1; }
  if (0==1){info=SDPConeCheckData(sdpcone);}
  
  info = DSDPSolve(dsdp); 
  if (info){ printf("Numberical error\n"); return 1; }
  info=DSDPComputeX(dsdp);DSDPCHKERR(info);
  
  if (0==1){ /* Look at the solution */
    int n; double *xx;
    info=SDPConeGetXArray(sdpcone,0,&xx,&n);
  }

  info = DSDPDestroy(dsdp);
  free(Edges);
  
  return 0;
} /* main */

/*!
\fn int SetThetaData(DSDP dsdp, SDPCone sdpcone, int nodes, int edges, EdgeMat Edge[]);
\param dsdp the solver
\param sdpcone the semidefinite cone
\param nodes number of nodes in graph
\param edges number of edges in graph
\param Edge edges in graph
\brief Given a graph, formulate Lovasz problem and set data.
\ingroup Examples
\sa LovaszTheta
 */
int SetThetaData(DSDP dsdp, SDPCone sdpcone, int nodes, int edges, EdgeMat Edge[]){

  int i,info;

  /* Create data matrices -  these are all custom types */

  /* The c matrix has all elements equal to 1.0 */
  info=OneMatOpsInitialize(&onematops);
  info=SDPConeAddDataMatrix(sdpcone,0,0,nodes,'P',&onematops,0);

  /* For each edge connecting nodes i and j, X(i,j)=X(j,i)=0 */
  info=EdgeMatOperationsInitialize(&edgematop);
  for (i=0; i<edges; i++){
    info = SDPConeAddDataMatrix(sdpcone,0,i+1,nodes,'P',&edgematop,(void*)&Edge[i]);
    info = DSDPSetDualObjective(dsdp,i+1,0.0);
    info = DSDPSetY0(dsdp,i+1,0.0);
    if (info) return info; 
  }

  /* The trace of X must equal 1.0 */
  info = IdentitymatOperationsInitialize(&identitymatops);
  info = SDPConeAddDataMatrix(sdpcone,0,edges+1,nodes,'P',&identitymatops,0);
  info = DSDPSetDualObjective(dsdp,edges+1,1.0);
  info = DSDPSetY0(dsdp,edges+1,-10*nodes-1);

  /* The initial point y is feasible and near the central path */
  info = DSDPSetR0(dsdp,0);
  
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
  for (i=0; i<edges; i++){ E[i].v1=0; E[i].v2=0; }

  while(!feof(fp) && gotone){
    thisline[0]='\0';
    fgets(thisline,BUFFERSIZ,fp); line++;
    info = ParseGraphline(thisline,&row,&col,&value,&gotone);
    if (gotone && k<edges && 
	col < nodes && row < nodes && col >= 0 && row >= 0){
      if (row > col){i=row;row=col;col=i;}
      if (row == col){}
      else { E[k].v1=row; E[k].v2=col; k++;}
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


/* SPecial Matrix where each edge represents a constraint matrix */
static int EdgeMatDestroy(void*);
static int EdgeMatView(void*);
static int EdgeMatVecVec(void*, double[], int, double *);
static int EdgeMatDot(void*, double[], int, int, double *);
static int EdgeMatGetRank(void*, int*, int);
static int EdgeMatFactor(void*);
static int EdgeMatGetEig(void*, int, double*, double[], int, int[], int*);
static int EdgeMatAddRowMultiple(void*, int, double, double[], int);
static int EdgeMatAddMultiple(void*, double, double[], int, int);
static int EdgeMatGetRowNnz(void*, int, int[], int*, int);

static int EdgeMatDestroy(void* AA){
  return 0;
}
static int EdgeMatVecVec(void* A, double x[], int n, double *v){
  EdgeMat*E=(EdgeMat*)A;
  *v=2*x[E->v1]*x[E->v2];
  return 0;
}
static int EdgeMatDot(void* A, double x[], int nn, int n, double *v){
  EdgeMat*E=(EdgeMat*)A;
  int k=E->v2*(E->v2+1)/2 + E->v1;
  *v=2*x[k];
  return 0;
}
static int EdgeMatView(void* A){
  EdgeMat*E=(EdgeMat*)A;
  printf(" Row: %d, Column: %d\n",E->v1,E->v2);
  printf(" Row: %d, Column: %d\n",E->v2,E->v1);
  return 0;
}
static int EdgeMatFactor(void* A){
  return 0;
}
static int EdgeMatGetRank(void *A, int*rank, int n){
  *rank=2;
  return 0;
}
static int EdgeMatGetEig(void*A, int neig, double *eig, double v[], int n,int spind[], int *nind){
  EdgeMat*E=(EdgeMat*)A;
  double tt=1.0/sqrt(2.0);
  memset((void*)v,0,(n)*sizeof(double)); 
  memset((void*)spind,0,(n)*sizeof(int)); 
  if (neig==0){
    v[E->v1]=tt;v[E->v2]=tt;*eig=1;
    spind[0]=E->v1;spind[1]=E->v2; *nind=2;
  } else if (neig==1){
    v[E->v1]=-tt;v[E->v2]=tt;*eig=-1;
    spind[0]=E->v1;spind[1]=E->v2; *nind=2;
  } else { *eig=0;*nind=0;}
  return 0;
}
static int EdgeMatGetRowNnz(void*A, int nrow, int nz[], int *nnzz, int n){
  EdgeMat*E=(EdgeMat*)A;
  if (nrow==E->v1){ nz[E->v2]++; *nnzz=1;} 
  else if (nrow==E->v2){nz[E->v1]++; *nnzz=1;} 
  else {*nnzz=0;}
  return 0;
}
static int EdgeMatAddRowMultiple(void*A, int nrow, double dd, double rrow[], int n){
  EdgeMat*E=(EdgeMat*)A;
  if (nrow==E->v1){ rrow[E->v2]+=dd;} 
  else if (nrow==E->v2){rrow[E->v1]+=dd;} 
  return 0;
}
static int EdgeMatAddMultiple(void*A, double dd, double vv[], int nn, int n){
  EdgeMat*E=(EdgeMat*)A;
  int k=E->v2*(E->v2+1)/2 + E->v1;
  vv[k]+=dd;
  return 0;
}
static int EdgeMatFNorm(void*A, int n, double *fnorm){
  *fnorm=2.0;
  return 0;
}
static int EdgeMatCountNonzeros(void*A, int *nnz, int n){
  *nnz=1;
  return 0;
}
static const char *datamatname="THETA EDGE MATRIX";
static int EdgeMatOperationsInitialize(struct  DSDPDataMat_Ops* edgematoperator){
  int info;
  if (edgematoperator==NULL) return 0;
  info=DSDPDataMatOpsInitialize(edgematoperator); if (info){ return info;}
  edgematoperator->matfactor1=EdgeMatFactor;
  edgematoperator->matgetrank=EdgeMatGetRank;
  edgematoperator->matgeteig=EdgeMatGetEig;
  edgematoperator->matvecvec=EdgeMatVecVec;
  edgematoperator->matrownz=EdgeMatGetRowNnz;
  edgematoperator->matdot=EdgeMatDot;
  edgematoperator->matfnorm2=EdgeMatFNorm;
  edgematoperator->matnnz=EdgeMatCountNonzeros;
  edgematoperator->mataddrowmultiple=EdgeMatAddRowMultiple;
  edgematoperator->mataddallmultiple=EdgeMatAddMultiple;
  edgematoperator->matdestroy=EdgeMatDestroy;
  edgematoperator->matview=EdgeMatView;
  edgematoperator->matname=datamatname;
  edgematoperator->id=25;
  return 0;
}

/* SPecial Matrix where all elements equal negative one */
static int OneMatDestroy(void*);
static int OneMatView(void*);
static int OneMatVecVec(void*, double[], int, double *);
static int OneMatDot(void*, double[], int, int, double *);
static int OneMatGetRank(void*, int*, int);
static int OneMatFactor(void*);
static int OneMatGetEig(void*, int, double*, double[], int, int[], int*);
static int OneMatRowNnz(void*, int, int[], int*, int);
static int OneMatAddRowMultiple(void*, int, double, double[], int);
static int OneMatAddMultiple(void*, double, double[], int,int);


static int OneMatFactor(void*A){return 0;}
static int OneMatGetRank(void *A, int *rank, int n){*rank=1;return 0;}
static int OneMatFNorm2(void*AA, int n, double *fnorm2){*fnorm2=1.0*n*n;return 0;}
static int OneMatCountNonzeros(void*AA, int *nnz, int n){*nnz=n*n;return 0;}
static int OneMatDot(void* A, double x[], int nn, int n, double *v){
  double dtmp=0.0;
  int i,j;
  for (i=0;i<n;i++){
    for (j=0;j<i;j++,x++){dtmp+= (*x);}
    dtmp+= (*x);x++;
  }
  *v=-2*dtmp;
  return 0;
}
static int OneMatVecVec(void* A, double x[], int n, double *v){
  double dtmp=0.0;
  int i;
  for (i=0; i<n; i++){dtmp+=x[i];}
  *v=-dtmp*dtmp;
  return 0;
}
static int OneMatAddMultiple(void*A, double ddd, double vv[], int nn, int n){
  int i,j;
  for (i=0;i<n;i++){
    for (j=0;j<i;j++,vv++){(*vv)+=-ddd;}
    (*vv)+=-ddd; vv++;
  }
  return 0;
}
static int OneMatAddRowMultiple(void*A, int nrow, double ddd, double row[], int n){
  int i;
  for (i=0;i<n;i++){row[i] -= ddd;}
  row[nrow] -= ddd;
  return 0;
}
static int OneMatGetEig(void*A, int neig, double *eig, double v[], int n, int spind[], int *nind){
  int i;
  if (neig==0){ *eig=-1; for (i=0;i<n;i++){ v[i]=1.0; spind[i]=i;} *nind=n; 
  } else {  *eig=0; for (i=0;i<n;i++){ v[i]=0.0; } *nind=0;
  }
  return 0;
}
static int OneMatRowNnz(void*A, int row, int nz[], int *nnz, int n){
  int i;
  for (i=0;i<n;i++){ nz[i]++; }
  *nnz=n;
  return 0;
}
static int OneMatView(void* AA){
  printf("Every element of the matrix is the same: -1\n");
  return 0;
}
static int OneMatDestroy(void* A){return 0;}

static const char *mat1name="THETA ALL ELEMENTS EQUAL -ONE";
static int OneMatOpsInitialize(struct  DSDPDataMat_Ops* mat1ops){
  int info;
  if (mat1ops==NULL) return 0;
  info=DSDPDataMatOpsInitialize(mat1ops); DSDPCHKERR(info);
  mat1ops->matfactor1=OneMatFactor;
  mat1ops->matgetrank=OneMatGetRank;
  mat1ops->matgeteig=OneMatGetEig;
  mat1ops->matvecvec=OneMatVecVec;
  mat1ops->matdot=OneMatDot;
  mat1ops->mataddrowmultiple=OneMatAddRowMultiple;
  mat1ops->mataddallmultiple=OneMatAddMultiple;
  mat1ops->matdestroy=OneMatDestroy;
  mat1ops->matview=OneMatView;
  mat1ops->matrownz=OneMatRowNnz;
  mat1ops->matfnorm2=OneMatFNorm2;
  mat1ops->matnnz=OneMatCountNonzeros;
  mat1ops->id=18;
  mat1ops->matname=mat1name;
  return 0;
}

/* Special Matrix for the Identity */
static int IdentityMatDestroy(void*);
static int IdentityMatView(void*);
static int IdentityMatVecVec(void*, double[], int, double *);
static int IdentityMatDot(void*, double[], int, int, double *);
static int IdentityMatGetRank(void*, int*,int);
static int IdentityMatFactor(void*);
static int IdentityMatGetEig(void*, int, double*, double[], int, int[], int*);
static int IdentityMatAddRowMultiple(void*, int, double, double[], int);
static int IdentityMatAddMultiple(void*, double, double[], int, int);
static int IdentityMatGetRowNnz(void*, int, int[], int*, int);

static int IdentityMatDestroy(void* AA){return 0;}
static int IdentityMatFNorm2(void* AA, int n, double *v){*v=1.0*n;return 0;}
static int IdentityMatGetRank(void *AA, int*rank, int n){*rank=n;return 0;}
static int IdentityMatFactor(void*A){return 0;}
static int IdentityMatCountNonzeros(void*A, int *nnz, int n){*nnz=n;return 0;}
static int IdentityMatVecVec(void* AA, double x[], int n, double *v){
  int i;
  *v=0;
  for (i=0;i<n;i++){ *v+=x[i]*x[i]; }
  return 0;
}
static int IdentityMatDot(void* AA, double x[], int nn, int n, double *v){
  int i;
  double vv=0;
  for (i=0;i<n;i++){ vv+=x[((i+1)*(i+2))/2-1];}
  *v = 2*vv;
  return 0;
}
static int IdentityMatView(void* AA){
  printf("Identity matrix: All Diagonal elements equal 1.0\n");
  return 0;
}
static int IdentityMatGetEig(void*AA, int neig, double *eig, double v[], int n, int spind[], int *nind){
  if (neig<0 || neig>=n){ *eig=0; return 0;} 
  memset((void*)v,0,(n)*sizeof(double)); 
  *eig=1.0;  v[neig]=1.0; spind[0]=neig; *nind=1;
  return 0;
}
static int IdentityMatGetRowNnz(void*A, int nrow, int nz[], int *nnzz, int n){
  if (nrow>=0 && nrow < n){
    *nnzz=1; nz[nrow]++;
  } else { *nnzz=0;    
  }
  return 0;
}
static int IdentityMatAddRowMultiple(void*A, int nrow, double dd, double rrow[], int n){
  rrow[nrow] += dd;return 0;
}
static int IdentityMatAddMultiple(void*A, double dd, double vv[], int nn, int n){
  int i;
  for (i=0;i<n;i++){ vv[(i+1)*(i+2)/2-1] += dd;}
  return 0;
}

static const char *eyematname="THETA IDENTITY MATRIX";
static int IdentitymatOperationsInitialize(struct  DSDPDataMat_Ops* imatops){
  int info;
  if (imatops==NULL) return 0;
  info=DSDPDataMatOpsInitialize(imatops); if (info){return info;}
  imatops->matfactor1=IdentityMatFactor;
  imatops->matgetrank=IdentityMatGetRank;
  imatops->matgeteig=IdentityMatGetEig;
  imatops->matvecvec=IdentityMatVecVec;
  imatops->matrownz=IdentityMatGetRowNnz;
  imatops->matdot=IdentityMatDot;
  imatops->matfnorm2=IdentityMatFNorm2;
  imatops->matnnz=IdentityMatCountNonzeros;
  imatops->mataddrowmultiple=IdentityMatAddRowMultiple;
  imatops->mataddallmultiple=IdentityMatAddMultiple;
  imatops->matdestroy=IdentityMatDestroy;
  imatops->matview=IdentityMatView;
  imatops->id=12;
  imatops->matname=eyematname;
  return 0;
}
