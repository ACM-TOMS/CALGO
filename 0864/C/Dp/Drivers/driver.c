
/* ======================================================================
      3DBPP PROBLEMS, Silvano Martello, David Pisinger, Daniele Vigo
   ====================================================================== */

/* This code generates instances for the three-dimensional bin-packing 
 * problem and solves them using the 3DBPP algorithm by Martello, Pisinger,
 * Vigo. 
 *
 * A description of the test instances is found in the following papers:
 *   S. Martello, D. Pisinger, D. Vigo, E. den Boef, J. Korst (2003)
 *   "Algorithms for General and Robot-packable Variants of the
 *    Three-Dimensional Bin Packing Problem"
 *   submitted TOMS.
 *
 *   S.Martello, D.Pisinger, D.Vigo (2000)
 *   "The three-dimensional bin packing problem"
 *   Operations Research, 48, 256-267
 *
 * The algorithm can either read an instance from a file or randomly
 * generate 10 instances.
 *
 * If the instance is read from a file, five arguments should be given:
 *   filename  A filename in which the test instance is found. The format 
 *             of the file is:
 *                n W H D
 *                w_1 h_1 d_1
 *                :
 *                w_n h_n d_n
 *             where 
 *                n is the number of items, 
 *                W,H,D is the size of the bin, and 
 *                w_j,h_j,d_j is the size of box j.
 *   nodelimit maximum number of decision nodes to be explored in the
 *             main branching tree. If set to zero, the algorithm will
 *             run until an optimal solution is found (or timelimit or
 *             iterlimit is reached). Measured in thousands (see IUNIT).
 *   iterlimit maximum number of iterations in the ONEBIN algorithm
 *             which packs a single bin. If set to zero, the algorithm will
 *             run until an optimal solution is found (or timelimit or
 *             nodelimit is reached). Measured in thousands (see IUNIT).
 *   timelimit Time limit for solving the problem expressed in seconds.
 *             If set to zero, the algorithm will run until an optimal
 *             solution is found; otherwise it terminates after timelimit
 *             seconds with a heuristic solution. 
 *   nodeused  returns the number of branch-and-bound nodes investigated,
 *             measured in thousands (see constant IUNIT in 3dbpp.c).
 *   iterused  returns the number of iterations in ONEBIN algorithm,
 *             measured in thousands (see constant IUNIT in 3dbpp.c).
 *   timeused  returns the time used in miliseconds
 *   packingtype 
 *             Desired packing type. If set to zero, the algorithm will
 *             search for an optimal general packing; if set to one, it
 *             will search for a robot packing.
 *             will search for a robot packing.
 * If the code should randomly generate 10 instances, seven arguments 
 * should be given:
 *   n         The size of the instance, i.e., number of boxes.
 *   bindim    The size of the bin, typically 40-100.
 *   type      An integer saying which randomly generated instance should
 *             be generated. A value between 1-9 selects one of the instance 
 *             types described in the above papers. Value 10-11 generates
 *             1D and 2D instances.
 *   nodelimit as above
 *   iterlimit as above
 *   timelimit as above
 *   packingtype
 *             as above
 *
 * Results are written to standard output
 * 
 * (c) Copyright 1998, 2003, 2005
 *
 *   David Pisinger                        Silvano Martello, Daniele Vigo
 *   DIKU, University of Copenhagen        DEIS, University of Bologna
 *   Universitetsparken 1                  Viale Risorgimento 2
 *   Copenhagen, Denmark                   Bologna, Italy
 * 
 * This code can be used free of charge for research and academic purposes 
 * only. 
 */

#define RANDOMTESTS    10     /* Number of test to run for each type */
#define MAXBOXES      201     /* Max number of boxes plus one */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>


/* ======================================================================
				   macros
   ====================================================================== */

#define srand(x)     srand48x(x)
#define randm(x)     (lrand48x() % (long) (x))
#define rint(a,b)    (randm((b)-(a)+1) + (a))

#define TRUE  1           /* logical variables */
#define FALSE 0

#define VOL(i)                 ((i)->w * (stype) (i)->h * (i)->d)
#define DIF(i,j)               ((int) ((j) - (i) + 1))


/* ======================================================================
				 type declarations
   ====================================================================== */

typedef int           boolean; /* logical variable      */
typedef int           ntype;   /* number of states,bins */
typedef short         itype;   /* can hold up to W,H,D  */
typedef int           stype;   /* can hold up to W*H*D  */

typedef int (*funcptr) (const void *, const void *);

/* box record */
typedef struct irec {
  ntype    no;           /* box number */
  itype    w;            /* box x-size */
  itype    h;            /* box y-size */
  itype    d;            /* box z-size */
  itype    x;            /* optimal x-position */
  itype    y;            /* optimal y-position */
  itype    z;            /* optimal z-position */
  ntype    bno;          /* bin number */
  stype    vol;          /* volume of box */
} box;


/* ======================================================================
				global variables
   ====================================================================== */

int TESTS;


/* ======================================================================
                                binpack3d
   ====================================================================== */

void binpack3d(int n, int W, int H, int D,
               int *w, int *h, int *d,
               int *x, int *y, int *z, int *bno,
               int *lb, int *ub,
               int nodelimit, int iterlimit, int timelimit,
               int *nodeused, int *iterused, int *timeused,
               int packingtype);


/* =======================================================================
                                random
   ======================================================================= */

/* to ensure that the same instances are generated on all machines */
/* here follows C-versions of SRAND48, and LRAND48.  */

unsigned int _h48, _l48;

void srand48x(long s)
{
  _h48 = s;
  _l48 = 0x330E;
}

int lrand48x(void)
{
  _h48 = (_h48 * 0xDEECE66D) + (_l48 * 0x5DEEC);
  _l48 = _l48 * 0xE66D + 0xB;
  _h48 = _h48 + (_l48 >> 16);
  _l48 = _l48 & 0xFFFF;
  return (_h48 >> 1);
}


/* ======================================================================
				specialbin
   ====================================================================== */

void specialbin(box *f, box *l, int W, int H, int D)
{
  box *m, *i, *j, *k;
  int w, h, d;

  if (f == l) { f->w = W; f->h = H; f->d = D; return; }
  if (DIF(f,l) == 5) {
    w = W/3; h = H/3; i = f+1; j = f+2; k = f+3; 
    i->w =   W-w; i->h =     h; i->d = D;
    j->w =     w; j->h =   H-h; j->d = D;
    k->w =   W-w; k->h =     h; k->d = D; 
    l->w =     w; l->h =   H-h; l->d = D;
    f->w = W-2*w; f->h = H-2*h; f->d = D;
    return;
  }

  m = f + (l-f) / 2;
  for (;;) {
    switch (rint(1,3)) {
      case 1: if (W < 2) break;
              w = rint(1,W-1); 
  	      specialbin(f, m, w, H, D);
	      specialbin(m+1, l, W-w, H, D);
	      return;
      case 2: if (H < 2) break;
	      h = rint(1,H-1); 
	      specialbin(f, m, W, h, D); 
	      specialbin(m+1, l, W, H-h, D);
	      return;
      case 3: if (D < 2) break;
	      d = rint(1,D-1);
	      specialbin(f, m, W, H, d); 
	      specialbin(m+1, l, W, H, D-d);
	      return;
    }
  }
}


/* ======================================================================
				randomtype
   ====================================================================== */

void randomtype(box *i, int W, int H, int D, int type)
{
  itype w, h, d, t;

  if (type <= 5) { /* Martello, Vigo */
    t = rint(1,10); if (t <= 5) type = t;
  }

  w = h = d = 0;
  switch (type) {
    /* Martello, Vigo */
    case  1: w = rint(1,W/2);   h = rint(2*H/3,H); d = rint(2*D/3,D); break;
    case  2: w = rint(2*W/3,W); h = rint(1,H/2);   d = rint(2*D/3,D); break;
    case  3: w = rint(2*W/3,W); h = rint(2*H/3,H); d = rint(1,D/2);   break;
    case  4: w = rint(W/2,H);   h = rint(H/2,H);   d = rint(D/2,D);   break;
    case  5: w = rint(1,W/2);   h = rint(1,H/2);   d = rint(1,D/2);   break;

    /* Berkey, Wang */
    case 6: w = rint(1,10);     h = rint(1,10);    d = rint(1,10);    break;
    case 7: w = rint(1,35);     h = rint(1,35);    d = rint(1,35);    break;
    case 8: w = rint(1,100);    h = rint(1,100);   d = rint(1,100);   break;

    /* 2D and 1D cases */
    case 10: w = rint(1,W);     h = rint(1,H);     d = D;             break;
    case 11: w = rint(1,W/2);   h = H;             d = D;             break;
  }
  i->w = w; i->h = h; i->d = d;
  i->x = 0; i->y = 0; i->z = 0;
}


/* ======================================================================
				allgood
   ====================================================================== */

boolean allgood(stype totvol, box *f, box *l)
{
  box *j, *m;
  stype vol;

  for (vol = 0, j = f, m = l+1; j != m; j++) {
    if ((j->w < 1) || (j->h < 1) || (j->d < 1)) return FALSE;
    vol += VOL(j);
  }
  return (vol == totvol);
}


/* ======================================================================
				maketest
   ====================================================================== */

void maketest(box *f, box *l, itype *W, itype *H, itype *D,
	      stype bdim, int type)
{
  register box *i, *k, *m;
  int no;

  /* set bin dimensions */
  *W = bdim; *H = bdim; *D = bdim;

  /* make maxtypes box types */
  for (i = f, m = l+1, no = 1; i != m; i++, no++) {
    randomtype(i, *W, *H, *D, type);
    i->no = no;
  }

  /* guillotine cut three bins */
  if (type == 9) {
    no = DIF(f,l)/3;
    k = f + no; m = k + no;
    for (;;) {
      specialbin(f  , k, *W, *H, *D); 
      specialbin(k+1, m, *W, *H, *D); 
      specialbin(m+1, l, *W, *H, *D); 
      if (allgood(3*bdim*bdim*bdim, f, l)) break;
    }
  }
}


/* ======================================================================
				readtest
   ====================================================================== */

int readtest(box *tab, itype *W, itype *H, itype *D, char *file)
{
  FILE *in;
  box *i;
  int n, w, h, d;

  in = fopen(file, "r");
  if (in == NULL) { printf("wrong filename"); exit(-1); }

  fscanf(in,"%d %d %d %d", &n, &w, &h, &d);
  *W = w; *H = h; *D = d;
  for (i = tab; i < tab+n; i++) {
    fscanf(in,"%d %d %d", &w, &h, &d);
    i->w = w; i->h = h; i->d = d;
  }
  fclose(in);
  return n;
}


/* ======================================================================
                                printboxes
   ====================================================================== */

void printboxes(int n, int W, int H, int D, int *w, int *h, int *d, 
                int *x, int *y, int *z, int *bno)
{
  int i;

  printf("%d (%d,%d,%d)\n", n, W, H, D);
  for (i = 0; i < n; i++) {
    printf("%2d (%2d %2d %2d) : Bin %2d (%2d, %2d, %2d)\n",
           i,w[i],h[i],d[i],bno[i],x[i],y[i],z[i]);
  }
}


/* ======================================================================
                                prepareboxes
   ====================================================================== */

void prepareboxes(box *f, box *l, int *w, int *h, int *d)
{
  box *i;
  int k;

  for (i = f, k = 0; i != l+1; i++, k++) {
    w[k] = i->w; h[k] = i->h; d[k] = i->d;
  }
}


/* ======================================================================
				main
   ====================================================================== */

int main(int argc, char *argv[])
{
  int v, n;
  itype W, H, D;
  int bdim, type, packingtype;
  int nodelimit, iterlimit, timelimit, nodeused, iterused, timeused;
  box tab[MAXBOXES];
  int w[MAXBOXES], h[MAXBOXES], d[MAXBOXES];
  int x[MAXBOXES], y[MAXBOXES], z[MAXBOXES], bno[MAXBOXES];
  int ub, lb, solved, gap, sumnode, sumiter;
  double time, sumtime, deviation, sumub, sumlb, sumdev;
  char file[1000];

  if (argc == 6) {
    strcpy(file, argv[1]);
    nodelimit   = atoi(argv[2]);
    iterlimit   = atoi(argv[3]);
    timelimit   = atoi(argv[4]);
    packingtype = atoi(argv[5]);
    printf("3DBPP PROBLEM %s %d %d %d %d\n", 
           file, nodelimit, iterlimit, timelimit, packingtype);
    n = readtest(tab, &W, &H, &D, file);
    bdim = 0;
    type = 0;
    TESTS = 1;
  } 
  if (argc == 8) {
    n           = atoi(argv[1]);
    bdim        = atoi(argv[2]);
    type        = atoi(argv[3]);
    nodelimit   = atoi(argv[4]);
    iterlimit   = atoi(argv[5]);
    timelimit   = atoi(argv[6]);
    packingtype = atoi(argv[7]);
    TESTS = RANDOMTESTS;
  }
  if ((argc != 6) && (argc != 8)) {
    printf("3DBPP PROBLEM\n");
    printf("n = ");
    scanf("%d", &n);
    printf("bindim = ");
    scanf("%d", &bdim);
    printf("type = ");
    scanf("%d", &type);
    printf("nodelimit = ");
    scanf("%d", &nodelimit);
    printf("iterlimit = ");
    scanf("%d", &iterlimit);
    printf("timelimit = ");
    scanf("%d", &timelimit);
    printf("packingtype = ");
    scanf("%d", &packingtype);
    TESTS = RANDOMTESTS;
  }

  printf("3DBPP PROBLEM %d %d %d %d %d %d %d\n", 
          n, bdim, type, nodelimit, iterlimit, timelimit, packingtype);
  sumnode = sumiter = sumtime = sumdev = sumub = sumlb = gap = solved = 0;
  for (v = 1; v <= TESTS; v++) {
    srand(v+n); /* initialize random generator */
    if (type != 0) maketest(tab, tab+n-1, &W, &H, &D, bdim, type);
    prepareboxes(tab, tab+n-1, w, h, d);
    binpack3d(n, W, H, D, w, h, d, x, y, z, bno, &lb, &ub, 
              nodelimit, iterlimit, timelimit, 
              &nodeused, &iterused, &timeused, 
              packingtype);
    time = timeused * 0.001;
    printf("%2d : lb %2d z %2d node %9d iter %9d time %6.2f\n", 
            v, lb, ub, nodeused, iterused, time); 
    sumnode += nodeused;
    sumiter += iterused;
    sumtime += time; 
    gap += ub - lb;
    deviation = (ub - lb) / (double) lb;
    sumdev += deviation;
    sumub += ub;
    sumlb += lb;
    if (lb == ub) solved++;
    if (type == 0) printboxes(n, W, H, D, w, h, d, x, y, z, bno);
  }
  printf("n           = %d\n", n);
  printf("bdim        = %d\n", bdim);
  printf("type        = %d\n", type);
  printf("nodelimit   = %d\n", nodelimit);
  printf("iterlimit   = %d\n", iterlimit);
  printf("timelimit   = %d\n", timelimit);
  printf("packingtype = %d\n", packingtype);
  printf("solved      = %d\n", solved);
  printf("ub          = %.1f\n",   sumub / (double) TESTS);
  printf("lb          = %.1f\n",   sumlb / (double) TESTS);
  printf("gap         = %.1f\n",     gap / (double) TESTS);
  printf("dev         = %.2f\n",  sumdev / (double) TESTS);
  printf("nodes       = %.0f\n", sumnode / (double) TESTS);
  printf("iterations  = %.0f\n", sumiter / (double) TESTS);
  printf("time        = %.2f\n", sumtime / (double) TESTS);

  return 0; /* correct termination */
}

