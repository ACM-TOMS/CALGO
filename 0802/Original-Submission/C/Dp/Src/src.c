/* ALC2DLIB-Version 1.0 implemented by WH 22.3.99, all rights reserved */
/*please report your experience  to whoer@statistik.wu-wien.ac.at      */
/*                               or hormannw@boun.edu.tr               */

/* This file is the main program file       alc2d.c
   for explanations see:                    alc2d.h */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "alc2d.h"

extern double UNIFORM(void);

/* The following constants should only be changed, if
the double type is shorter than 64 bit*/

#define EPS 1.e-13 /*tolerance for linecut*/
#define EPS1 1.e-13 /*tolerance for checksameside*/
#define EPS2 1.e-13 /*tolerance for checkvector point*/
#define EPSTW 1.e-6 /*tolerance for a=0 (numerically tested)*/
#define VOLHATMAX 1.e99 /*large positive value at most HUGE_VAL*/



void pr_2dlib_def(void)
/*prints the values of the defined constants*/
{ printf("OUTPUT %d SUNC %d MEMCONTROL %d MEMOUTPUT %d VOLHATMAX %e\n",
	 OUTPUT,SUNC,MEMCONTROL,MEMOUTPUT,VOLHATMAX);
  printf("NC %d CFGT %d EPS %e EPS1 %e EPS2 %e\n",NC,CFGT,EPS,EPS1,EPS2);
}



/*some output-macros for debugging*/
#if (MEMOUTPUT==2)
#define MPR2(a) printf a
#else
#define MPR2(a) 
#endif
#if (MEMOUTPUT>=1)
#define MPR1(a) printf a
#else
#define MPR1(a) 
#endif
#if (OUTPUT==2)
#define OPR2(a) printf a
#else
#define OPR2(a) 
#endif
#if (OUTPUT>=1)
#define OPR1(a) printf a
#else
#define OPR1(a) 
#endif
#if (OUTPUT>=0)
#define OPR0(a) printf a
#else
#define OPR0(a) 
#endif
 
/*some macro-definitions*/
#define OK 1239999
#ifndef NULL
#define NULL (void *) 0
#endif
#define Arrlet2(a,b) do{a[0]=b[0];a[1]=b[1];} while(0)
#define Arrlet2zero(a) do{a[0]=0.;a[1]=0.;} while(0)
#define Arrlet3(a,b) do{a[0]=b[0];a[1]=b[1];a[2]=b[2];} while(0)
#define Arradd2(a,b,c) do{a[0]=b[0]+c[0];a[1]=b[1]+c[1];} while(0)
#define Arrsub2(a,b,c) do{a[0]=b[0]-c[0];a[1]=b[1]-c[1];} while(0)
#define Arrminuslet2(a,c) do{a[0]= -c[0];a[1]= -c[1];} while(0)
#define Arrchangesign2(a) do{a[0]= -a[0];a[1]= -a[1];} while(0)
#define Letstru(a,b) memcpy(&(a),&(b),sizeof(b))
#define Sqr(c) ((c)*(c))
#define Sqrdistance(a,b) (Sqr(b[0]-a[0])+Sqr(b[1]-a[1]))
#define MIN(a,b) ((a)<(b)?a:b)
#define MAX(a,b) ((a)>(b)?a:b)
#define Hatvalue(t,c) (t[1]==0.||t[2]==0.?t[0]+t[1]*c[0]+t[2]*c[1]:t[1]*t[2]*(t[0]/(t[1]*t[2])+c[0]/t[2]+c[1]/t[1]))
#define Directionvector(a,b) do{a[0]= -b[2];a[1]= b[1];} while(0)
/*a[2]..result, b[3]...line, stores the direction vector of b[3] in a[2]*/ 

typedef struct 
{ int nc/*number of corners of polygon+unbounded edges*/,
      nca/*number of corners for which memory was allocated up to now*/,
      update/* 1..update necessary, 0.. no update necessary*/,
      open/*0..closed pol,1..open,2..twoparallel lines,3..one line,4..empty*/,
      index/*index number of that polygon in structure TOTAL*/ ;
  double design[2]/*tangential point*/,
         tang[3]/*tangential space*/,
         (*corner)[2]/*pointer to the array of the vertices in
                                             consecutive order*/,
         sc[2]/*corner with the highest z-value*/,
         s /*value of the tangential plane over sc, i.e. highest value
                                              of the hat over the polygon*/,
         rot[2]/*rotation that rotates the negative gradient of tang
                                                            in x-direction*/,
         a/*slope of the hat in the direction of the gradient*/;
} POLYGON; /*the convex polygon associated with one point of contact */

typedef struct 
{ int open/*kind of region,0(=normal a negative),1(=normal a positive),2 open*/,
      ip/*index of polygon containing this generator region*/;
  double point[2],/*coordinates of the point translated into origin*/
       rot[2],/*rotation to rotate negative gradient into (1,0)*/
       k[2],/*slopes of the two lines from origin*/
       s,/** z-coordinate of point translated into origin*/
       a,/* h(x1,x2)=exp(a*x+const) is the hat over rotated generator region*/
       b,/*for open<=1 the generator region is in the interval (0,b)
           for open=2 b is broadness of the parallel strip*/
       vol,/**volume of this part*/
       vols;/**cumulated volume*/
} GEN; /*All informations about one generator region*/

typedef struct
{ int n,/* current number of touching points*/
      nn,/* maximum number of touching points*/
    gennr,/*current number of generator regions*/
    ngennr,/*current maximum number of generator regions*/
    m,/*current length of guidetable*/
    mm,/*maximum length of guide table*/
    *q/* pointer to the integer array for the guide table*/;
long ok;/*used to check if the correct pointer was passed to sample 
         ok==OK means that the setup was successfull*/
  GEN *gen/*[nngen] pointer to the array of generator regions*/;
  POLYGON *p;/*[nn+1] pointer to the array of POLYGON;
               only necessary as long as new points of contact are added*/
  double (*h)(double x[2],double param[]),
         (*hx)(double x[2],double param[]),(*hy)(double x[2],double param[]),
/*pointers to the functions evaluating the log of the density and
  the partial derivatives of the density*/
    *param,/*pointer to the double array with the current parameter settings*/
    oldhatvol;/*total volume below old hat
                used to check against rounding errors*/
} TOTAL; /*contains all information necessary for set-up and sampling*/


/*******************************************************************/
#if(MEMCONTROL)
/*definition of some auxiliary functions used for memory controll*/
typedef struct{void *adress;unsigned long int size;}HEAP;
static int nrofap=0/*number of allocated pointers*/,
           maxsize=0/*maximal size of allocated memory observed so far*/,
           actsize=0/*current size of allocated memory*/;
static HEAP heap[1000];

void checkheap()
/*computes the current values of the variables maxsize and actsize*/
{ int i;
  actsize=0;
  if(nrofap)  for(i=0;i<nrofap;i++) actsize+=heap[i].size;
  if(maxsize<actsize) maxsize=actsize;
}

void printheap()
/*prints a message containing the current state of the allocated memory*/
{
  printf("there are %d allocated pointers with %d bytes memory,\n", 
                  nrofap,actsize);
  printf("maximal heapsize up to now %d bytes  (adr.of first p. %p)\n",
                  maxsize,heap[0].adress);
}

void *mallocwh(long int size)
/*used instead of malloc() to perform the allocation checks*/
{ void *poi;

  if(size<=0) 
  { if (size<0) printf("ERROR!!mallocwh was called with size <0\n");
    return(NULL);
  }
  poi=malloc(size);
  if(poi!=NULL)
  { heap[nrofap].adress=poi;
    heap[nrofap].size=size;
    nrofap++;
  }
  else
  { printf("Error!!mallocwh memory was not allocated successfully!!\n");
    exit(1);
  }
  checkheap();
#if (MEMCONTROL==2)
  printheap();
#endif
  return(poi);
}

void freewh(void *poi)
/*used instead of free() to perform the allocation checks*/
{ int i,found= -1;
  
  for(i=0;i<nrofap;i++) 
  { if(heap[i].adress==poi) found=i;
    if(found>=0) Letstru(heap[i],heap[i+1]);
  }
  if(found>=0)
    { nrofap--;
      free(poi);
      checkheap();
#if (MEMCONTROL==2)
      printheap();
#endif
    }
  else if(poi!=NULL)
     printf("ERROR!! memtest: a pointer freed which is not in the list\n");
}
void *reallocwh(void *poiold,long int size)
/*used instead of realloc() to perform the allocation checks*/
{ void *poinew;
  int i,found= -1; 


  if(size<=0) printf("ERROR!!reallocwh was called with size <=0\n");
  for(i=0;i<nrofap;i++) 
  { if(heap[i].adress==poiold) found=i;
    if(found>=0) Letstru(heap[i],heap[i+1]);
  }
  if(found>=0) nrofap--;
  else if(poiold!=NULL)
  printf("ERROR!! memtest: a pointer reallocated which is not in the list\n");
  poinew=realloc(poiold,size);
  if(poinew!=NULL)
  { heap[nrofap].adress=poinew;
    heap[nrofap].size=size;
    nrofap++;
  }
  else
  { printf("Error!!reallocwh memory was not allocated successfully!!\n");
    if(size>0) exit(1);
  }
  checkheap();
#if (MEMCONTROL==2)
  printheap();
#endif
  return(poinew);
}
#else
void *mallocwh(int size)
{ return(malloc(size));}
void freewh(void *poi)
{ free(poi);}
void *reallocwh(void *poiold,int size)
{ return(realloc(poiold,size));}
#endif
/**************************************************************************/
/*other utilities for debugging*/

void printgen(GEN *g)
/*prints all information of one generator region*/
{ printf("open %d ip %d point(%4.3e|%4.3e) rot(%4.3e|%4.3e) k(%4.3e|%4.3e)\n",
    g->open,g->ip,g->point[0],g->point[1],g->rot[0],g->rot[1],g->k[0],g->k[1]);
  printf("s %4.3e  a  %4.3e  b  %4.3e  vol %4.3e vols %4.3e\n",
    g->s,g->a,g->b,g->vol,g->vols);
}
void printpoly(POLYGON *poly)
/*prints all information of one polygone*/
{ int i;
  printf("n:%d update:%d open:%d design:(%4.3e|%4.3e) tang:(%4.3e|%4.3e|%4.3e)\n",
         (*poly).nc,(*poly).update,(*poly).open,(*poly).design[0],
         (*poly).design[1],(*poly).tang[0],(*poly).tang[1],(*poly).tang[2]);
  for(i=0;i<=((*poly).nc);i++) 
      printf("%d: (%4.3e|%4.3e)\n",i,(*poly).corner[i][0],(*poly).corner[i][1]);
}
/****************************************************************************/
#if (MEMCONTROL>=1) 
static int ncmax=0;/*maximal number of corners observed for a polygon*/
int returnncmax(void)
/*returns the current value of ncmax*/
{ return(ncmax);}
#else
int returnncmax(void) /*returns -1 as MEMCONTOL==0*/
{ printf("Caution! function returnncmax() does not work correctly as MEMCONTROL not >= 1\n");
  return(-1);}
#endif
/***********************************************************************/
double gg2neu(double ab)
{ double u,h1,h2;
  long int fak=1,i;  

  u=UNIFORM()*(exp(ab)*(ab-1.)+1.);
  h1=ab*ab;
  h2=h1*0.5;
/*printf("h2 %f\n",h2);*/
  if(u<h2) return(sqrt(2.*u));
  for(i=3;i<100;i++)
  {
    h1*=ab;
    fak*=(i-2);
    u-=h2;
    h2=h1/(i*fak);
/*printf("h2 %f i %d fak %d\n",h2,i,fak);*/
    if(u<h2) {/*printf("%d \n",i);*/return(exp((1./i)*log(i*fak*u)));}
  }
  return(ab);
}



#define A (-0.6821555671)
#define EDE 0.3678794411714423216
#define ED2E2 0.06766764161830634595
#define E 2.7182818284590452354
#define ED2E 0.1839397205857211608
#define AGES 1.089775550197435
#define EMULT 1.15747

double gg2inter(double a,double b)
/*generates variates with density proportional to f(x)=x exp(a x)in (0,b)*/
{ double x,ab,tlx,al,col,v,ugr,u,h;

  if(a>0)
    {/*Fall a>0*/
      ab=a*b;
      if(ab<1.) return gg2neu(ab)/a;
      tlx=ab*0.65;
      al=1.+1./tlx;
      col=exp(al*(0.-ab))-1.;
      while(1)
	{
          x=ab+log(UNIFORM()*col+1.)/al;
	  h=x-tlx-al*(x-tlx);
          v=UNIFORM()*tlx;
	  if (v<=x*(1.+h)) return(x/a);
	  if (v<=x*(1.+h*(1.+h*(0.5+h*(1./6.))))) return(x/a);
	  if (v<=x*exp(h)) return(x/a);
	}
    }
  else if(a<0)/*Fall a<0*/
    { ab= -a*b;
      if(ab<EDE) ugr=ab*ab*0.5;
      else if(ab<1.6803) ugr= -ED2E2+EDE*ab;
      else ugr=AGES+(EMULT/A)*exp(A*ab);
      while(1)
	{
          u=UNIFORM()*ugr;
          v=UNIFORM();
	  if (u<ED2E2) 
	    { x=sqrt(2.*u);
	      if (v<=1-x) return(-x/a);
	      if (v<=exp(-x)) return(-x/a);
	    }
	  else if (u<0.5504801833820683) 
	    { x=ED2E+E*u;
	      if (v<=0.25464*E) return(-x/a);
	      if (v<=x*exp(-x+1.)) return(-x/a);
	    }
	  else 
	     { x=log((AGES-u)/(-EMULT/A))*(1./A);
	       h= x*(1.+A);
               if (v*EMULT<=x*(1.-h*(1.-h*(0.5-h*(1./6.))))) return(-x/a);
	       if (v*EMULT<=x*exp(-h)) return(-x/a);
	    }
	}
    }
  else /*Fall a=0  f(x) proportional to x in (0,b)*/
    return(sqrt(UNIFORM())*b);
}

/*******************************************************************/
double marg2(GEN *ge)
/*generates for generator with open=2 random variate of marginal in x-direction*/
{ double u,ak10,a,b;

  ak10=fabs(ge->k[1]-ge->k[0]);
  a=fabs(ge->a);
  b=fabs(ge->b);
  u=UNIFORM()*(ak10/a+b);
  if(u<b)
    {/*parallel strip*/
      u=u/b;/*recycling*/
      return(-log(u)/a);
    }
  else
    {/*"single-vertex"*/
      u=(u-b)/(ak10/a);/*recycling*/
      return(-log(UNIFORM()*UNIFORM())/a);
    }
}
/*******************************************************************/
/*Rotations, described by the matrix    rotx -roty
                                        roty  rotx
are stored in the form (rotx,roty)*/
/*lines are stored according to the equation 0=line[0]+x line[1]+y line[2]*/
/*planes are stored as z=ebene[0]+x ebene[1] + y ebene[2] */
/***********************************************************************/
int projinter(double plane0[3],double plane1[3],double result[3])
/*Computes the (x/y)-projection of the intersection line of 2 planes
  0 is returned if the two planes are nearly parallel 
  result is computed and stored in results and 1 returned otherwise */
{ 
  result[0]=(plane1[0]-plane0[0]);
  result[1]=(plane1[1]-plane0[1]);
  result[2]=(plane1[2]-plane0[2]);
  if(fabs(result[1])<EPS && fabs(result[2])<EPS) return(0);
  else return(1);
}

/***********************************************************************/
int linecut(double g0[3],double g1[3],double result[2])
/*intersection of the 2 lines: g0[3] and g1[3] is computed 
  result is unchanged and 0 returned if the lines are nearly parallel 
  result is computed and stored in results and 1 returned otherwise */
{ double det,help;
  det=g0[1]*g1[2]-g0[2]*g1[1];
/*following if checks against division through 0!*/ 
  if(fabs(g0[1]*g1[2])<1.) help=1.;
  else help=g0[1]*g1[2];
  if ((fabs(det)<EPS)||(fabs(det/help)<EPS)) return(0);
  else 
  { result[0]=(g0[2]*g1[0]-g0[0]*g1[2])/det;
    result[1]=(g0[0]*g1[1]-g1[0]*g0[1])/det;
    return(1);
  }
}
/*****************************************************************/
int checksameside(double g[3],double p1[2],double p2[2])
/*returns
  1...if the points p1 and p2 are off the line g and on the same side of g 
  0...if the points p1 and p2 are off the line g and on different sides of g
  -1... if one of the points is very close to the line g*/
{ int help1,help2;
  double dhelp1,dhelp2;
  dhelp1=g[0]+p1[0]*g[1]+p1[1]*g[2];
  help1=(dhelp1>0.);
  dhelp2=g[0]+p2[0]*g[1]+p2[1]*g[2];
  help2=(dhelp2>0.);
  if (fabs(dhelp1)<EPS1||fabs(dhelp2)<EPS1) {return(-1);}
  if ((help1&&help2)||(!help1&&!help2)) {return(1);}
  return(0);
}
/*****************************************************************/
void changelineeq(double pv[2],double p[2],int flag,double res[3])
/* puts the standard form of a line defined by point-vector or
   point-point into res
   pv[2]...point or direction-vector of the line
   p[2]... point of the line
   flag : 0.. pv is point, 1..pv is direction vector*/
{ 
  if(flag) {res[2]= -pv[0];res[1]=pv[1];} 
  else {res[2]=pv[0]-p[0];res[1]=p[1]-pv[1];} 
  res[0]= -(res[1]*p[0]+res[2]*p[1]);
}
/*****************************************************************/
void intersect(double p1[2],double p2[2],double g[3],double res[2])
/*computes into res[2] the intersection of the line through p1[2] p2[2] 
  and the line g[3]. */ 
{ double helpg[3];
  changelineeq(p1,p2,0,helpg);
  linecut(g,helpg,res);
}
/*****************************************************************/
int checksamesinfty
        (double vec[2],double p[2],double tp[2],double g[3],double res[2])
/*computes for the line (defined by (direction)vec[2] and point p[2])
  and line g[3] the intersection point (and stores it in res[2]); 
  then determines if the infinity point defined by vec p and the 
  tangential point tp[2] are on the same side of g or not.
  1...if the points are on the same side of g 
  0...if the points are on different sides of g
  -1... if the infinity points is very close to the line g
  it is assumed that the tangential point tp is allways off the line*/
{ double helpg[3],helpp[2];

  changelineeq(vec,p,1,helpg);
  if (linecut(g,helpg,res))
  /*the case that the two lines are not parallel*/
    { helpp[0]=res[0]+10.*vec[0]; 
      helpp[1]=res[1]+10.*vec[1]; 
      return( checksameside(g,helpp,tp) );
    }
  else /*case two lines are allmost parallel*/ 
     return( checksameside(g,p,tp) );
}
/*****************************************************************/
 void insertcorner(POLYGON *p,int nr,double point[2])
/*inserts into the array corner of polygon a new point at index nr
  POL.nc is increased by 1*/
{ int i;

  if (p->nc>=p->nca-1)/*da erst unten p->nc erhoeht wird*/ 
  {   
    MPR1(("corner at adress %p with %d points is\n",p->corner,p->nca));
    p->nca+=NC;
    p->corner=reallocwh(p->corner,p->nca*2*sizeof(double));
    if(p->corner==NULL) printf("error realloc no space in heap\n");
    MPR1(("(realloc)is enlargened to %d points (%d bytes) at adress %p\n",p->nca,p->nca*2*sizeof(double),p->corner));
  }
  for(i=p->nca-1;i>nr;i--) Arrlet2(p->corner[i],p->corner[i-1]);
  Arrlet2(p->corner[nr],point);
  p->nc++;
#if (MEMCONTROL>=1)
  if(ncmax<p->nc) ncmax=p->nc;
#endif
}
/*****************************************************************/
#define POL (*poly)

int checkcorner(POLYGON (*poly),double line[3],double helpp[2][2],int a)
/*calls the two versions of check same side depending on the value
 of i, the number of the corner   
  helpp is used to store the intersection points between vector and line,
  which are computed in the function checksamesinfty()*/
{ int n;
  n=POL.nc;
  if (POL.open&&a==0) return(
    checksamesinfty(POL.corner[0],POL.corner[1],POL.design,line,helpp[0]));
  else if(POL.open&&a==POL.nc-1) return(
    checksamesinfty(POL.corner[n-1],POL.corner[n-2],POL.design,line,helpp[1]));
  else return(checksameside(line,POL.corner[a],POL.design));
}
/*****************************************************************/
void intersectpolline(POLYGON (*poly),double line[3],double helpp[2][2],
           int i1,int i2,double intp[2])
/*calls the different versions to compute the intersection with the
  polygon between the neighbouring corners with numbers i1 and i2
  result is stored in intp
  the result of the intersection of vector with line was allready computed in
  the function checksamesinfty() and stored in helpp[2][2]*/ 
{ int k; 
  if (i1>i2) {k=i1;i1=i2;i2=k;}       
  if (POL.open) 
    { if(i1==0&&i2==POL.nc-1) 
	{ intp[0]= -line[2];
	  intp[1]= line[1];
          return;
	}
      else if(i1==0) { Arrlet2(intp,helpp[0]);return;}
      else if(i1==POL.nc-2) { Arrlet2(intp,helpp[1]);return;}
    }
  intersect(POL.corner[i1],POL.corner[i2],line,intp);
}
/*****************************************************************/
void deletecorners(POLYGON (*poly),int corem[])
/*deletes all corners of the polygon with number i and corem[i]<=0*/

{ int i,j=0;

  while(corem[j]>=1&&j<POL.nc) j++;
  for(i=j+1;i<POL.nc;i++) 
    if(corem[i]>=1) 
    { Arrlet2(POL.corner[j],POL.corner[i]);
      j++;
    }
  POL.nc=j;
}

/*****************************************************************/
void checkvectorsigns(POLYGON (*poly))
/* checks the direction of the vectors of open polygons
   and changes them if the touching point is not inside the polygon*/ 
{ double hl[3],hp[2];
  int i;
/*change vector 0*/

  if (POL.nc==3) changelineeq(POL.corner[2],POL.corner[1],1,hl);
/*  else changelineeq(POL.corner[2],POL.corner[1],0,hl);*/
/*stattdessen fuer Spezialproblem falls viele Punkte sehr eng beisammen*/
  else 
  { for(i=2;(i<POL.nc-1)&&(Sqrdistance(POL.corner[i],POL.corner[1])<1.e-10);i++);
    if(i==POL.nc-1) { 
      changelineeq(POL.corner[i],POL.corner[1],1,hl);}
    else changelineeq(POL.corner[i],POL.corner[1],0,hl);
  }
  Arradd2(hp,POL.corner[1],POL.corner[0]);  
  if(!checksameside(hl,hp,POL.design)) Arrchangesign2(POL.corner[0]);  
/*cange vector n-1*/
  if (POL.nc==3) changelineeq(POL.corner[0],POL.corner[1],1,hl);
/*else changelineeq(POL.corner[POL.nc-2],POL.corner[POL.nc-3],0,hl);*/
/*stattdessen fuer Spezialproblem falls viele Punkte sehr eng beisammen*/
  else 
  { for(i=POL.nc-3;(i>0)&&
        (Sqrdistance(POL.corner[POL.nc-2],POL.corner[i])<1.e-10);i--);
    if(i==0) changelineeq(POL.corner[0],POL.corner[POL.nc-2],1,hl);
    else changelineeq(POL.corner[POL.nc-2],POL.corner[i],0,hl);
  }
  Arradd2(hp,POL.corner[POL.nc-2],POL.corner[POL.nc-1]);  
  if(!checksameside(hl,hp,POL.design)) Arrchangesign2(POL.corner[POL.nc-1]);  
}
/*****************************************************************/
void nocornerupdate(POLYGON (*poly),double line[3])
/*a polygon (which has up to now no corners) is updated by the line[]*/
{ double hline1[3],hline2[3],help1,help2,help3;

  POL.update=1;/*polygon is allways considered as updated*/
  switch(POL.open)
  { case 4:/*polygon up to now empty*/
    { POL.open=3;
      Arrlet3(POL.corner[0],line);/*corner is used to store the line*/
      break;
    }
    case 3:/*up to now only one line stored in corner*/
    { 
      Arrlet3(hline1,POL.corner[0]);
      if (linecut(hline1,line,POL.corner[1]))        
      { POL.nc=3;/*case that the two lines have a cutting point*/
	POL.open=1;
        Directionvector(POL.corner[0],line);
        Directionvector(POL.corner[2],hline1);
	checkvectorsigns(poly);
      }
      else/*hline1 and line are allmost parallel
	    thus it is computed if both lines or only one line are stored*/
      { if (fabs(hline1[1])>EPS) help1=hline1[0]*line[1]/hline1[1];
        else help1=hline1[0]*line[2]/hline1[2];
        help2= -line[1]*POL.design[0]-line[2]*POL.design[1];
	if ((help1<help2&&help2<line[0])||(help1>help2&&help2>line[0]))
	{ POL.open=2;/*case: both lines are stored*/
	  /*corner[2] is used to store the 2nd parallel line*/
	  Arrlet3(POL.corner[2],line);
	}
        else if((help2<line[0]&&line[0]<help1)||(help2>line[0]&&line[0]>help1))
	  Arrlet3(POL.corner[0],line);/*case: line is stored, hline1 discarde*/
      /*else case: as hline1 is kept nothing is changed*/
      }
      break;
    }    
    case 2:/*up to now two parallel lines stored in corner*/
    { 
      Arrlet3(hline1,POL.corner[0]);
      Arrlet3(hline2,POL.corner[2]);
      if (linecut(hline1,line,POL.corner[1]))
      /* case: the new line is not parallel to the two old ones*/
      { POL.nc=4;
	POL.open=1;
	linecut(hline2,line,POL.corner[2]);
     /*the direction vector (-line[2],line[1]) is the vector twoards infinity*/
        Directionvector(POL.corner[0],hline1);
	Arrlet2(POL.corner[3],POL.corner[0]);          
	checkvectorsigns(poly);
      }
      else/*case: all 3 lines hline1 hline2 and line are allmost parallel
              thus it is computed if one line has to be replaced*/
      { if (fabs(hline1[1])<EPS) help1=hline1[0]*line[1]/hline1[1];
        else help1=hline1[0]*line[2]/hline1[2];
	if (fabs(hline2[1])<EPS) help2=hline2[0]*line[1]/hline2[1];
	else help2=hline2[0]*line[2]/hline2[2];
        help3= -line[1]*POL.design[0]-line[2]*POL.design[1];
	if ((help3>line[0]&&line[0]>help1)||(help3<line[0]&&line[0]<help1))
          /*case: hline1 is replaced by line*/
        { POL.open=2;
	  Arrlet3(POL.corner[0],line);
        }
        else if((help3>line[0]&&line[0]>help2)||(help3<line[0]&&line[0]<help2))
          Arrlet3(POL.corner[2],line); /*case: hline2 is replaced by line*/  
      /*else case: as line is outside the region nothing is changed*/
      }
    }
  }
}
/**************************************************************/
void polyupdate(POLYGON (*poly),double line[3])
/*the polygon is updated by the line*/
{ double helpp[2][2],ip1[2],ip2[2];
  int i,* corem,nrem=0,ip1nr,ip2nr;

  if(POL.open>1) nocornerupdate(poly,line);/*special case of no corner*/
  else
  { /*cases: POL.open=1 (=open) or 0 (=closed)*/
    corem=mallocwh(sizeof(int)*POL.nca);
    MPR2(("for corem %d byte allocated at adress %p\n",sizeof(int)*POL.nca,corem));

    for(i=0;i<POL.nc;i++)
    { corem[i]=checkcorner(poly,line,helpp,i);
/*corem<=0... corner is removed, corem=1... corner remains*/ 
      if(corem[i]<=0) nrem++;
/*nrem...number of corners that are removed*/
    }
    if(nrem>0)
    { for(i=0;i<POL.nc;i++)
      { if(corem[i]==1 && corem[(i+1)%POL.nc]<=0)
        { ip1nr=(i+1)%POL.nc;/*ip1nr is now the index of the (cyclic) first 0*/
          intersectpolline(poly,line,helpp,i,ip1nr,ip1);          
        }
        else if(corem[i]<=0 && corem[(i+1)%POL.nc]==1)
        { ip2nr=i;/*ip2nr is now the index of the (cyclic) last zero*/
          intersectpolline(poly,line,helpp,i,(i+1)%POL.nc,ip2);          
        }
      }
      if(nrem==1)
      { insertcorner(poly,ip1nr,ip1);
        Arrlet2(POL.corner[ip2nr+1],ip2);        
      }
      else
      { Arrlet2(POL.corner[ip1nr],ip1);        
        Arrlet2(POL.corner[ip2nr],ip2);
        corem[ip1nr]=1;
        corem[ip2nr]=1;
        deletecorners(poly,corem);
      }
      POL.update=1;
      if (ip1nr>ip2nr) POL.open=0;
      if (POL.open) checkvectorsigns(poly);
    }
    MPR1(("corem adress %p is freed\n",corem));
    freewh(corem);
    corem=NULL;
  }
}


/*****************************************************************/
void rotation(double tang[3],double rot[2],double *a)
/* computes a rotation such that the negative gradient of tang[3] points
   in direction of the positiv x-axis; the result is stored in rot[2]
   - the length of the gradient is stored in a*/
{ double t1,t2,st; 
  t1= tang[1];
  t2= tang[2];
  st= sqrt(t1*t1+t2*t2);
  (*a)=-st;
  if (st<=EPSTW) 
  { (*a)=0.;
    rot[0]=1.;
    rot[1]=0.;}
  else 
  { rot[0]= -t1/st;
    rot[1]= t2/st;}
}


/*??perhaps faster when implemented as a macro??*/
void rotate(double rot[2],double vektor[2],double res[2])
/* rotates vektor[2] with rotation rot[2]; result is stored in res[2]*/ 
{ res[0]=rot[0]*vektor[0]-rot[1]*vektor[1];
  res[1]=rot[1]*vektor[0]+rot[0]*vektor[1];
}
void invrotate(double rot[2],double vektor[2],double res[2])
/*rotates vektor[2] with inverse rotation of rot[2];result is stored in res[2]*/
{ res[0]=rot[0]*vektor[0]+rot[1]*vektor[1];
  res[1]= -rot[1]*vektor[0]+rot[0]*vektor[1];
}

/*****************************************************************/
void volumenc(GEN *ge)
/*determines volume of closed generator region (i.e. open<=1)*/
#define GE (*ge)
{ double eda;
  if (fabs(GE.a)>EPSTW)
  { eda=1./GE.a;
 /*   GE.vol= fabs(GE.k[1]-GE.k[0])*exp(GE.s)*eda*(
                                    (GE.b-eda)*exp(GE.a*GE.b)+eda);
*/    GE.vol= (fabs(GE.k[1]-GE.k[0])*(exp(GE.s)*eda*eda+
                                    eda*(GE.b-eda)*exp(GE.a*GE.b+GE.s)));

  }
  else GE.vol=exp(GE.s)*fabs(GE.k[1]-GE.k[0])*GE.b*GE.b*0.5;
               /*case of constant hat*/
}
void volumenu(GEN *ge)
/*determines volume of unbounded generator region (i.e. open==2)*/
#define GE (*ge)
{  
  if (GE.a<-1.e-60)
    GE.vol=(exp(GE.s)/GE.a)*(fabs(GE.k[1]-GE.k[0])/GE.a-fabs(GE.b));
  else { GE.vol= -99999999.;
OPR1(("Probelm: Vol infinit, generator region of open==2 with a>=0\n"));
       }
/*Remark vol of "eineck"=exp(s)*fabs(k1-k0)/(a*a)
         vol of parallel strip= -exp(s)*b/a  where b is the broadness of the strip*/
}
#undef GE

/*****************************************************************/
void  ngenplusplus(int *ngen,int ngennr,double vol)
/*advances ngen, the current number of generator regions only
  if vol>0.    otherwise ngen is not advanced and thus the 
  generator region is overwritten*/
{  if(vol>0.) (*ngen)++;
/*due to the if generator regions with vol<=0. are overwritten*/
#if (MEMCONTROL>=1)
  if(*ngen>ngennr) 
    printf("Error!!!! ngen %d ngennr %d too  regions!!\n",*ngen,ngennr);
#endif
 }
/*****************************************************************/
void comp1stgen(POLYGON *poly,GEN gen[],int *ngen,int ngennr,
		double hpoi0[2],double hpoi1[2],double hpoi2[2],int output)
/*compute everything of the first generator region of an open or closed
  triangular.
 hpoi[0] are the original coordinates of the top, 
 hpoi[1] the translated roated coordinates of the second point, 
 hpoi[2] the coordinates of the rotated vector pointing from the top
 output ... 0 absolute no output*/
{
  gen[*ngen].ip=POL.index;
  gen[*ngen].open=0;/*closed generator region, a<0*/
  gen[*ngen].b=hpoi1[0];
  if (hpoi1[0]>0.) gen[*ngen].k[0]=hpoi1[1]/hpoi1[0];
  else gen[*ngen].k[0]= -99999999.;
  if (hpoi2[0]>0.) gen[*ngen].k[1]=hpoi2[1]/hpoi2[0];
  else {gen[*ngen].k[0]= -99999999.;
      if(output)
	printf("!!!!Error problem!! hpoi2[1]/hpoi2[0] %e/%e\n",hpoi2[1],hpoi2[0]);
  }
/*k[0] and k[1] can be interchanged as well*/
  Arrlet2(gen[*ngen].point,hpoi0);
  Arrlet2(gen[*ngen].rot,POL.rot);
  gen[*ngen].s=Hatvalue(POL.tang,hpoi0);
  gen[*ngen].a=POL.a;
  /*due to round-off errors it is possible that b<0. therefore*/
  if (gen[*ngen].b<EPS) gen[*ngen].vol=0.;
  else volumenc(&(gen[*ngen]));    
}


/*****************************************************************/
void compgen(POLYGON *poly,double point0[2],double point1[2],double opoint1[2],
                GEN gen[],int *ngen,int ngennr,int output)
/*For the triangle the generator regions are computed and stored in gen
  point0[2]...translated rotatd coordinates of second point of triangel
  point1[2]...translated rotatd coordinates of third point,
              (it is assumed that point0[0]<=point1[0])
  opoint1[2]...original coordinates of the third point
  int *ngen number of generator regions stored in gen up to now
  output ... 0 absolute no output*/
{ 

/*the first generator region of the triangel is calculated*/
  comp1stgen(poly,gen,ngen,ngennr,POL.sc,point0,point1,output);
  ngenplusplus(ngen,ngennr,gen[*ngen].vol);

  gen[*ngen].ip=POL.index;
  gen[*ngen].open=1;/*closed generator region a>0*/
  gen[*ngen].s=POL.s+POL.a*point1[0];
  Arrminuslet2(gen[*ngen].rot,POL.rot);
  gen[*ngen].b=point1[0]-point0[0];
  if(gen[*ngen].b>0.)gen[*ngen].k[0]= (point1[1]-point0[1])/gen[*ngen].b;
  else gen[*ngen].k[0]= 99999999.;
/*k[0] and k[1] can be interchanged as well*/
  if(point1[0]>0.) gen[*ngen].k[1]= point1[1]/point1[0];
  else {gen[*ngen].k[1]= 99999999.;
      if(output)
        printf("!!!!Error point1[1]/point1[0] %e/%e\n",point1[1],point1[0]);
  }
  Arrlet2(gen[*ngen].point,opoint1);
  gen[*ngen].a= -POL.a;
  /*due to round-off errors it is possible that b<0. therefore*/
  if (gen[*ngen].b<=0.) gen[*ngen].vol=0.;
  else volumenc(&(gen[*ngen]));    
  ngenplusplus(ngen,ngennr,gen[*ngen].vol);
}

/*****************************************************************/



int compgeno(POLYGON (*poly),GEN gen[],int *ngen,int ngennr,int output)
/*For an infinite biangle the generator regions are computed and stored in gen
  *ngen ...number of generator regions stored in gen up to now
function returns 1 if everything is ok 
                 0 if volume below the hat unbounded occurs
output ...0 absolute no output*/
{ double hpoi[5][2],helpp[2],hb,hs;


/*hpoi[0] are the original coordinates of the top of the biangle, 
 hpoi[1] the translated roated coordinates of the second point, 
 hpoi[2] the coordinates of the rotated vector pointing from the top
 hpoi[3] the vector pointing from the second point, 
 hpoi[4] the original coordinates of the second point
 in the case of a one angle the points hpoi[0] and hpoi[1] are identical!*/

  if(Hatvalue(POL.tang,POL.corner[1])>=Hatvalue(POL.tang,POL.corner[POL.nc-2]))
    { Arrlet2(hpoi[0],POL.corner[1]);
      Arrsub2(helpp,POL.corner[POL.nc-2],hpoi[0]);
      rotate(POL.rot,helpp,hpoi[1]);
      rotate(POL.rot,POL.corner[0],hpoi[2]);
      rotate(POL.rot,POL.corner[POL.nc-1],hpoi[3]);
      Arrlet2(hpoi[4],POL.corner[POL.nc-2]);
    }
    else
    { Arrlet2(hpoi[0],POL.corner[POL.nc-2]);
      Arrsub2(helpp,POL.corner[1],hpoi[0]);
      rotate(POL.rot,helpp,hpoi[1]);
      rotate(POL.rot,POL.corner[POL.nc-1],hpoi[2]);
      rotate(POL.rot,POL.corner[0],hpoi[3]);
      Arrlet2(hpoi[4],POL.corner[1]);
    }   

/*the first (bounded) generator region of the biangel is calculated*/
  comp1stgen(poly,gen,ngen,ngennr,hpoi[0],hpoi[1],hpoi[2],output);

/*the below values of the second generator region are computed using
  values from the first region; therefore *ngen is increased afterwards
  this important as the first region is deleted if vol==0 */
 
  if(POL.open==1 && POL.nc==3) hb=0.;/*case of oneangle*/
  else if(hpoi[1][0]>EPS)
          hb= gen[*ngen].b*(gen[*ngen].k[1]-gen[*ngen].k[0]);
  else hb= -hpoi[1][1];
  hs=gen[*ngen].s+POL.a*hpoi[1][0];
  ngenplusplus(ngen,ngennr,gen[*ngen].vol);

/*the second (unbounded) generator region of the biangel is calculated*/
  gen[*ngen].ip=POL.index;
  gen[*ngen].open=2;/*open generator region*/
  gen[*ngen].a= POL.a;
  gen[*ngen].s=hs;
  gen[*ngen].b=hb;
  Arrlet2(gen[*ngen].rot,POL.rot);

/*hier darf k[1] und k[0] nicht vertauscht werden!! Warum??*/
  if(hpoi[2][0]>1.e-15) gen[*ngen].k[1]=hpoi[2][1]/hpoi[2][0];
  else 
    {
      if(output)
      { OPR0(("Error!!compgeno2 Volume unbounded probably bad starting values\n"));
        OPR0(("hpoi[2][1] %e hpoi[2][0] %e \n",hpoi[2][1],hpoi[2][0]));}
      return(0);}
  if(hpoi[3][0]>1.e-15) gen[*ngen].k[0]=hpoi[3][1]/hpoi[3][0];
  else 
    {
      if(output)
      { OPR0(("Error!!compgeno3 Volume unbounded probably bad starting values\n"));
	OPR0(("hpoi[3][1] %e hpoi[3][0] %e \n",hpoi[3][1],hpoi[3][0]));}
      return(0);
    }
/*k[0] and k[1] can be interchangend as well*/
  Arrlet2(gen[*ngen].point,hpoi[4]);
  volumenu(&(gen[*ngen]));    
  ngenplusplus(ngen,ngennr,gen[*ngen].vol);
  return(1);
}

/*****************************************************************/


#define CP1(i) (i+1)%POL.nc /*cyclic +1  */
#define CM1(i) (POL.nc+i-1)%POL.nc /*cyclic -1  */


int triangulate(POLYGON (*poly),GEN gen[],int *ngen,int ngennr,int output)
/*the polygon (open=1 or 0) is decomposed into triangles and the generator
  regions are computed and stored in gen
  *ngen...number of generator regions stored in gen up to now
function returns 1 if everything is ok 
                 0 if volume below the hat unbounded occurs
output...0 absolute no output*/
{ int itop,i,poln;
  double z,zn,xmin,(*ctrp)[2]/*points to corners of transl. rotat. polygon*/,
         helpp[2],(*corn0p)[2]; 

  if(POL.open>1) return(0);/*as volume is certainly infinity*/

/* this is a trick that allows to make the part of the bounded triangels with
  only one code for bounded and unbounded polygons.
  Only the number of corners POL.nc must be changed and the array containing 
  the corners has to start with corner[1] instead of corner[0]*/
  if(POL.open==1) { poln=POL.nc-2;
		    corn0p=&(POL.corner[1]);}
  else  { poln=POL.nc;
	  corn0p=&(POL.corner[0]);}

/*first the corner itop with the maximal z-value of the hat is searched */
  /* if some z-values are equal, the corner with the minimal x-value is taken*/
  z= -1.e100;
  xmin= 1.e99;
  
  rotation(POL.tang,POL.rot,&(POL.a));
  if(fabs(POL.a)<EPS)
    { z=POL.tang[0];
      for(i=0;i<poln;i++)
       if(corn0p[i][0]<xmin)
          { itop=i; xmin=corn0p[i][0];}
    }
  else for(i=0;i<poln;i++)
    { 
      zn=Hatvalue(POL.tang,corn0p[i]);
      if (fabs(zn-z)<EPS && corn0p[i][0]<xmin) { itop=i; xmin=corn0p[i][0];}
      else if(zn>z) { itop=i; z=zn; xmin=corn0p[i][0];}
    } 
  Arrlet2(POL.sc,corn0p[itop]);
  POL.s=z;

  ctrp=mallocwh(sizeof(double)*2*POL.nca);
MPR1(("for ctrp %d byte allocated at adress %p\n",sizeof(double)*2*POL.nca, ctrp));
/*array ctrp Corners of Translated Rotated Polygon are computed*/
  Arrlet2zero(ctrp[0]);
  for(i=1;i<poln;i++)
  { Arrsub2(helpp,corn0p[(itop+i)%poln],POL.sc);
    rotate(POL.rot,helpp,ctrp[i]);
  }
  for(i=1;i<poln-1;i++) /* the polygon is triangulated and compgen is called*/
    if(ctrp[i][0]<ctrp[i+1][0]) 
       compgen(poly,ctrp[i],ctrp[i+1],corn0p[(itop+i+1)%poln],gen,ngen,ngennr,output);
    else compgen(poly,ctrp[i+1],ctrp[i],corn0p[(itop+i)%poln],gen,ngen,ngennr,output);
  MPR1(("ctrp adress %p is freed\n",ctrp));
  freewh(ctrp);
  ctrp=NULL;

  if(POL.open==1)
    if(!compgeno(poly,gen,ngen,ngennr,output)) return(0);

  return(1);
}
#undef POL

/****************/
#define PO t->p
void initpol(TOTAL (*t),int nr)
/*initialises the  polygon
  if(nr<t->nn) the values of the startpolygon are used as defaults*/
{ double hh,hhx,hhy;

  PO[nr].nca=NC;
  if(PO[nr].corner!=NULL) freewh(PO[nr].corner);
  PO[nr].corner=mallocwh(sizeof(double)*2*PO[nr].nca);
MPR1(("initpol for PO[%d].corner %d points (%d bytes) at adress %p\n",nr,PO[nr].nca,sizeof(double)*2*PO[nr].nca,PO[nr].corner));
  PO[nr].update=1;
  PO[nr].index=nr;
  if(nr==t->nn)
  { PO[nr].open=4;/*empty polygon*/
    PO[nr].nc=0;
  }
  else
  { PO[nr].open=PO[t->nn].open;
    PO[nr].nc=PO[t->nn].nc;
    memcpy(PO[nr].corner[0],PO[t->nn].corner[0],sizeof(double)*2*MAX(PO[t->nn].nc,2));
    /*this maximum for domain a half-plane*/
/* computation of the tangential plane in the point touch*/
    hh=(t->h)(PO[nr].design,t->param);
    hhx=(t->hx)(PO[nr].design,t->param);
    hhy=(t->hy)(PO[nr].design,t->param); 
    PO[nr].tang[0]=hh-PO[nr].design[0]*hhx-PO[nr].design[1]*hhy;
    PO[nr].tang[1]=hhx;
    PO[nr].tang[2]=hhy;
  }
}
/****************/
void removetp(TOTAL (*t),int nr)
/*removes polygon with number nr out of polygon array and closes the gap;
  (*t).n is updated aswell! nothing else is changed!*/
{ int i;
  ((*t).n)--;
MPR1(("Polygon %d deleted therefore corners at adress %p are freed\n",nr,PO[nr].corner));
/*As the polygon nr is deleted the memory allocated for corner must be freed*/
  freewh(PO[nr].corner);
  PO[nr].corner=NULL;
  for(i=nr;i<t->n;i++) Letstru(PO[i],PO[i+1]);
/*PO[*t.n] is marked as empty*/
  PO[(*t).n].corner=NULL;  
  PO[(*t).n].nc=0;
}
/****************/
void initclearpolarr(TOTAL *t,int n,int clear)
/*works on the polygon array
  n.. polygons 0 to n-1 are initialized and possibly cleared
  clear=0... initialization, clear=1.. clearing and initialization*/
{  int i;
   for(i=0;i<n;i++)
   { if (clear&&PO[i].corner!=NULL) 
     { MPR1(("PO[%d] at  adress %p is freed\n",i,PO[i].corner));
       freewh(PO[i].corner);
     }
     PO[i].nc=0;
     PO[i].corner=NULL;
   }
}
/****************/
int compute_tgen(TOTAL *t,int output)
/*computes the generator regions, 
  function returns 1 if everything is ok 
                   0 if volume below the hat unbounded occurs
  output 0, absolute no output*/
{ int genoldnr,i;
/*store old generator regions that can be reused in the beginning of gen*/
  genoldnr=t->gennr;
  t->gennr=0;
  if(genoldnr>0)  
    for(i=0;i<genoldnr;i++) 
      if(PO[t->gen[i].ip].update==0)/*that means this polygon was not updated*/
       { Letstru(t->gen[t->gennr],t->gen[i]);
	(t->gennr)++;
       }
  if(genoldnr==0)
/*last addpoint was interrupted with volume infinity, thus all generator
  regions must be recalculated*/
    for(i=0;i<t->n;i++) PO[i].update=1; 
/*computes the (approximate) total number of generatorregions
  using the formular #gen(1 polygon)=2*(n-2) */
  t->ngennr=t->gennr;
  for(i=0;i<t->n;i++) 
    if(PO[i].update) t->ngennr+=2*(PO[i].nc-2);
  if(t->ngennr<=0) t->ngennr=0;
  else
/*allocate the necessary memory for gen[]*/
  { if(t->gen==NULL) 
    { t->gen=mallocwh(sizeof(GEN)*t->ngennr);
      MPR1(("for t->gen %d byte allocated at adress %p\n",sizeof(GEN)*t->ngennr , t->gen));
    }
    else 
    { t->gen=reallocwh(t->gen,sizeof(GEN)*t->ngennr);
      MPR1(("for t->gen %d byte REallocated at adress %p\n",sizeof(GEN)*t->ngennr , t->gen));
    }
  }
/*triangulate changed polygons(automatically the generator regions
  of these polygons are stored in gen)*/
  for(i=0;i<t->n;i++) 
    if(PO[i].update) 
      if(!triangulate(&(PO[i]),t->gen,&(t->gennr),t->ngennr,output))
      { t->gennr=0;/*shows that the generator regions were not fully computed*/
	return(0); 
      }
  return(1);
}
/************************************/
void compute_cumvol_guidetable(TOTAL *t)
{ double volsh=0.;
  int i,j;

/*cumulated volumes are computed*/

  volsh=0.;
  for(i=0;i<(*t).gennr;i++) 
    { volsh+=(*t).gen[i].vol;
      (*t).gen[i].vols=volsh;
    }
/*we use the cummulated volume to check if there are apparently numerical
  ptoblems when computing the volumes*/
  if(volsh>t->oldhatvol*(1.+1.e-14))/*USE eps instead of constant*/
  { printf("Caution!! new volume below the hat was larger than before1\n");
    printf("Perhaps it is better to use fewer touching points!");
  }
  t->oldhatvol=volsh;
/*computation of guide-table*/
  (*t).q[0]=0;
  i=0;
/*the size of the guide table is proportional to (*t).n
  but during the setup a smaller guide table is used */
  if((*t).n==(*t).nn) (*t).m=(*t).mm;
  else (*t).m=CFGTSU*(*t).mm;
  for(j=1;j<(*t).m;j++)
  { while((*t).gen[i].vols/(*t).gen[(*t).gennr-1].vols<(double)j/(*t).m) i++;
    (*t).q[j]=i;
  }
}
/************************************/
int adddesignpoints(TOTAL *t,int output)
/*computes for one or more newly added points of contact the polygones
  and generator regions, and the necessary updates for the old polygones
  function returns 1 if everything is ok 
                   0 if volume below the hat unbounded occurs
  output 0.. absolut kein output*/
{ int i,j,nalt=0,initflag=0;
  long ok;
  double hger[3];

  for(i=0;i < t->n;i++)
  if(PO[i].corner==NULL) initpol(t,i);
  else  
  { nalt++;
    PO[i].update=0;/*up to now no update necessary*/
  }
/** for the new touching points (it is assumed that the old points are ordered
   at the beginning of the array) all possible lines (projections of the
   cutting lines of the tangentialebene with other tangential planes) are
   computed and possible updates of polygons are made*/
 do
  { initflag=0;
    for(i=t->n-1;i>=nalt;i--)
    for(j=0;j<i;j++)
      { 
	ok=projinter(PO[i].tang,PO[j].tang,hger);
/*ok=0 if the two tangential planes are allmost identical; 
  therefore one point is removed and all polygons are reinitialized*/
        if(ok)
	{ if (!initflag)
          { polyupdate(&(PO[i]),hger);
            polyupdate(&(PO[j]),hger);
	  }
	}
        else/*point i is removed*/
	{ initflag=1;
          OPR0(("\nPoint (%f|%f) was erased \n",PO[i].design[0],PO[i].design[1]));
          OPR0(("because of bad starting values or too many design-points\n"));
          removetp(t,i);
          j=i;
/*the loop continues with i-- to look for other points that should be removed*/
	}
      }
    if (initflag) 
    { nalt=0;
      for(i=0;i<t->n;i++) initpol(t,i);
    }
 }while(initflag);

#if (OUTPUT>=1)
 printf("Polygon-update finished!\n");
 for(i=0;i<(*t).n;i++) printpoly(&((*t).p[i]));
#endif

 ok=compute_tgen(t,output);/*trinagulates all polygons and computes the generator 
 regions, that are stored in t-gen*/

 if(!ok)return(0);
 
#if (OUTPUT>=1)
for(i=0;i<(*t).gennr;i++) {printf("%d:",i);printgen(&((*t).gen[i]));}
#endif


 compute_cumvol_guidetable(t);


#if (OUTPUT>=1)
for(i=0;i<(*t).gennr;i++) {printf("%d:",i);printgen(&((*t).gen[i]));}
#endif

OPR0(("\n points %d, volume below hat %e,(*t).gennr %d \n",(*t).n,(*t).gen[(*t).gennr-1].vols,(*t).gennr));
  return(1);
}/*end of adddesignpoints*/

/********************************************************/
void freepolarr(TOTAL *t)
{ int i;

  for(i=0;i<=(*t).nn;i++)
     if(t->p[i].corner!=NULL) 
     { 
MPR1(("corners von po[%d] at adress %p are freed\n",i,t->p[i].corner));
       freewh(t->p[i].corner);
       t->p[i].corner=NULL;
     }
MPR1(("the polygon array t->p adress %p is freed\n",t->p));
  freewh(t->p);
  t->p=NULL;
}
/****************/
void freesetup(void *t)
/*frees all memory that was assigned during the set-up*/
{ 
  TOTAL *t1;
  if(t!=NULL)
  {
    t1=t;
    if(t1->p!=NULL) freepolarr(t1);
MPR1(("structure of generator regions t->gen  adress %p is freed\n",t1->gen));
    freewh(t1->gen);
    t1->gen=NULL;
MPR1(("structure of guide table t->q  adress %p is freed\n",t1->q));
    freewh(t1->q);
    t1->q=NULL;
MPR1(("array param at  adress %p is freed\n",t1->param));
    freewh(t1->param);
    t1->param=NULL;
MPR1(("structure TOTAL at  adress %p is freed\n",t1));
    freewh(t1);
    t1=NULL;
#if (MEMCONTROL>=1&&OUTPUT>=0)
    printheap();
#endif
  }
}
/************************************/
void sample2d(double res[2],void *t1);
/****************/
void *setup(int nrsp,int nn,double sp[][2],int nruneq,double uneq[][3],
            int nrauxun,
            double (*h)(double x[2],double param[]),
            double (*hx)(double x[2],double param[]),
            double (*hy)(double x[2],double param[]),double param[],int nparam)
/* double sp[][2] .. array of the starting points,
 uneq[][3]..array of the equalities defining the domain of the distribution
equalities are stored in the form   
0 <= uneq[0][0] + uneq[0][1]*x + uneq[0][2]*y   
domain is defined by these equalities toegether with the starting points, 
which define, which of the different regions defined by the inequalities
is the domain (that, where the starting points lie inside.)

!it is assumed that all starting points in sp[] are inside the same region!

    nrsp  number of starting points
    nn    maximal number of points of contact
    nruneq  number of unequalities for domain of distribution

behind these inequalities it can be necessary to store inequalities
for an auxiliary domain (which must be a subset of the original domain)

   nrauxun .. total number of inequalities for auxiliary domain
        we must have    nrauxun >=nruneq
                        nrauxun== nruneq means that no auxiliary 
                        domain is provided

(*h) log of density function of 2-dim. distribution (can have an arbitrary
   number of double parameters)
(*hx) partial derivative of h with respect to the first variable
(*hy) partial derivative of h with respect to the second variable
param[] contains the current setting of the parameters of the distribution
nparam...... number of parameters of the distribution

*/
{ TOTAL *t,*taux;
  int i,ok,hic;
  double help[2];

  if(nn<nrsp||nn<3)
  { printf("ERROR!!setup: there are too few points of contact allowed!\n");
    return(NULL);}
  if(nrsp<1||(((nruneq==0&&nrsp<3)&&nrauxun==nruneq)))
  { printf("ERROR!!the setup was called with too few starting points.\n");
    printf("and without auxiliary domain. No setup was computed\n");
    return(NULL);
  }
  t=mallocwh(sizeof(TOTAL));
MPR1(("for t 1 TOTAL (%d bytes) allocated at adress %p\n",sizeof(TOTAL),t));
/*parameters are saved and necessary memory is allocated*/
  t->nn=nn;
  t->n=nrsp;
  t->mm=CFGT*nn;
  t->ngennr=0;
  t->gennr=0;
  t->h=h;
  t->hx=hx;
  t->hy=hy;
  t->param=mallocwh(sizeof(double)*nparam);
  memcpy(t->param,param,sizeof(double)*nparam);
MPR1(("for t->param %d doubles (%d bytes) at adress %p\n",nparam,sizeof(double)*nparam,t->param));
  t->oldhatvol=VOLHATMAX;
  t->gen=NULL;
  t->q=mallocwh(sizeof(int)*(t->mm));/*for guide table*/
MPR1(("for t->q %d ints (%d bytes) at adress %p\n",t->mm,sizeof(int)*(t->mm),t->q));
  t->p=mallocwh(sizeof(POLYGON)*(nn+1))/*in p[nn] the start polygon is stored*/;
MPR1(("for t->p %d polygons (%d bytes) at adress %p\n",nn+1,sizeof(POLYGON)*(nn+1),t->p));

  initclearpolarr(t,nn+1,0);
  for(i=0;i<nrsp;i++) Arrlet2(PO[i].design,sp[i]);
/*initialisation of Startpolygon*/
  initpol(t,t->nn);
  Arrlet2(PO[nn].design,sp[0]);
  for(i=0;i<nruneq;i++) polyupdate(&(PO[nn]),uneq[i]);

  ok=adddesignpoints(t,0);

/*if the volume below the hat is not bounded a second set-up using
  the auxiliary domain as domain of the distribution is started*/
  if(ok==0&&nrauxun>nruneq)
  { 
OPR0(("as volume unbounded setup is restarted using the auxiliary domain\n"));
    taux=setup(nrsp,nn,sp,nrauxun,uneq,nrauxun,h,hx,hy,t->param,nparam);
    if(taux!=NULL)
    { /*now we sample from the distribution over the auxiliary domain
        to get starting points for the original domain*/
      sample2d(help,taux);
OPR1(("first auxiliary sampling finished taux-n %d t-n %d\n",taux->n,t->n));
      do{
	if(taux->n==nn)/*in this case the polygon structure was dismissed*/
	{ printf("Error!!setup\n");
          printf("too few points allowed for auxiliary domain method\n");
          freesetup(taux);
          freesetup(t);
	  return(NULL);	      
	}
        for(i=0;i<taux->n;i++) Arrlet2(t->p[i].design,taux->p[i].design);
OPR1(("touching points transferred\n"));
/*generator regions and  polygons are cleared to make shure 
  that they are recalculated with the new touching points*/
	t->gennr=0;
        initclearpolarr(t,t->nn,1);
OPR1(("clearing finished  po[nn]\n"));
/*printpoly(&(t->p[t->nn]));*/
	t->n=taux->n;
        ok=adddesignpoints(t,0);
OPR1(("adddesignpoints(t) finished\n"));
        if(!ok) for(hic=0;(taux->n<t->n+1)&&(hic<1000);hic++)
                    { sample2d(help,taux);
OPR1(("auxiliary sampling finished taux-n %d t-n %d\n",taux->n,t->n));
                    }
	if(hic==1000)
	{ printf("Error!!setup, no additional design points are foun in auxiliary domain\n");
          printf("try to use larger auxiliary domain\n");
          freesetup(taux);
          freesetup(t);
	  return(NULL);	      
	}
      }while(!ok);
OPR0(("Starting values found, auxiliary domain points dismissed\n"));
      freesetup(taux);
    }
  }

/*Controll if computation of has worked correctly*/
  if (ok==0||t->gen[t->gennr-1].vols>=VOLHATMAX || t->gen[t->gennr-1].vols<0.) 
  {  
    printf("ErrorSETUP: volume below the hat unbounded!!\n");
    printf("you can try:\n");
    printf("a) to provide an auxilliary domain before starting set-up\n");
    printf("b) use an auxiliary domain and increas nn \n");
    printf("c) use more (or better) starting values\n");
    freesetup(t);
    return(0);
   }
  else 
  {  t->ok = OK;
     return(t);
   }
}


#undef PO
#if(MEMCONTROL>=1)
static unsigned long int widcount=0;
void setwidcount(unsigned long int newval)
/*sets the value of widcount*/
{ widcount=newval;}
unsigned long int returnwidcount(void)
/*returns the current value of widcount*/
{ return(widcount);}
#else /*the functions are defined but return -1 if MEMCONTROLL==1*/
void setwidcount(unsigned long int newval)
{ printf("Caution! function returnwidcount does not work correctly as MEMCONTROL not >= 1\n"); 
}
unsigned long int returnwidcount(void)
{ printf("Caution! function returnwidcount does not work correctly as MEMCONTROL not >= 1\n"); 
return(-1);}
#endif
/****************/
void sample2d(double res[2],void *t1)
/*samples from the 2-dimensional distribution, result is stored in res[2]
  t1 is the pointer that was returned by setup*/
{ TOTAL *t;
  int j;
  double zwert,u,v,helpp[2][2];

  t=t1;
  if(t==NULL)
  { printf("sample2d: es wurde Nullpointer uebergeben!!\n");
    return;}
  else if(t->ok!=OK)  
  { printf("sample2d: es wurde falscher pointer uebergeben!!\n");
    return;}
/*sampling routine*/
  while(1)
  { 
/*first the random index j of the region is generated using the guide table*/ 
    u=UNIFORM();
#if (MEMCONTROL>=1)
widcount++;
#endif
    j=t->q[(int)(t->m*u)];
    u*=t->gen[t->gennr-1].vols;
    while(t->gen[j].vols<u) j++;

/*depending on the "open" the pair "res" of the dominating distribution is 
  generated and the value of the log of the dominating density in that point 
  is evaluated and stored in zwert*/
    if(t->gen[j].open<=1)
    { res[0]=gg2inter(t->gen[j].a,t->gen[j].b);
      res[1]= UNIFORM()*res[0]*(t->gen[j].k[1]-t->gen[j].k[0])+
	res[0]*t->gen[j].k[0];
    }
    else
    { res[0]=marg2(&(t->gen[j]));
      res[1]=UNIFORM()*(t->gen[j].b+res[0]*(t->gen[j].k[1]-t->gen[j].k[0]))+
	res[0]*t->gen[j].k[0];
    }
    zwert=t->gen[j].s+res[0]*t->gen[j].a;
    invrotate(t->gen[j].rot,res,helpp[0]);
    Arradd2(res,helpp[0],t->gen[j].point);
/* the acceptance condition is tested*/
    v=log(UNIFORM());
    if (v<=(t->h)(res,t->param)-zwert) return;
/*in the case of rejection if the number of touching points is not exausted
    the rejected pair is taken as new touching point*/
    else if(t->n<t->nn)
    { t->n++;
      t->p[t->n-1].design[0]=res[0];
      t->p[t->n-1].design[1]=res[1];
/*adddesignpoints is computed only about ??? times*/
      if(t->n<10||(t->n%(t->nn/SUNC))==0||t->n==t->nn) adddesignpoints(t,1);
/*array of polygons is deleted as the number of points of contact was reached*/
      if((t->n==t->nn) &&(t->p!=NULL))
      { freepolarr(t);
	OPR0(("Polygon structure was freed\n"));
      }
    }
  }  
}






