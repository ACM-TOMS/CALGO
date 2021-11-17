/*
   ----------------------------------------------------------------
   File hov_reverse.c of ADOL-C version 1.6 as of January 1,   1995
   ----------------------------------------------------------------
   Contains the routine hov_reverse (higher-order-vector reverse 
   mode).

*/

#ifdef __cplusplus
extern "C" {
#endif


/* 
  There are four basic versions of the procedure `reverse', which
  are optimized for the cases of scalar or vector reverse sweeps
  with first or higher derivatives, respectively. In the calling
  sequence this distinction is apparent from the type of the
  parameters `lagrange' and `results'. The former may be left out
  and the integer parameters `depen', `indep', `degre', and `nrows'
  must be set or default according to the following matrix of
  calling cases. 

           no lagrange         double* lagrange     double** lagrange

double*   gradient of scalar   weight vector times    infeasible 
results   valued function      Jacobian product       combination

          ( depen = 1 ,         ( depen > 0 ,         
	    degre = 0 ,           degre = 0 ,              ------
	    nrows = 1 )           nrows = 1 )

double**  Jacobian of vector   weight vector times     weight matrix
results   valued function      Taylor-Jacobians        times Jacobian
           
	  ( 0 < depen           ( depen > 0 ,          ( depen > 0 ,
	      = nrows ,           degre > 0 ,            degre = 0 ,
	    degre = 0 )           nrows = 1 )            nrows > 0 )

double*** full family of         ------------          weigth matrix x
results   Taylor-Jacobians       ------------          Taylor Jacobians

*/ 

#include "dvlparms.h" /* Developers Parameters */

/* Necessary Includes */

#include "usrparms.h"
#include "oplate.h"
#include "taputil1.h"
#include "taputil2.h"
#include "taputil3.h"
#include "tayutil.h"

/* Static Locals */
static short tag;


/* External memory management routines from driversc.c */

double** myalloc2(int, int);
double*** myalloc3(int, int, int);

#define maxinc(a,b) if ((a) < (b)) (a) = (b)

/***************************************************************************/
/* Higher Order - Vector Reverse.                                          */
/***************************************************************************/
void hov_reverse(short tnum,/* tape id */
		 int depen,          /* consistency chk on # of dependents */
		 int indep,          /* consistency chk on # of independents */
		 int degre,          /* highest derivative degre */
		 int nrows,          /* # of Jacobian rows being calculated */
		 double **lagrange,  /* domain weight vector */
		 double ***results,  /* matrix of coefficient vectors */
		 short ** nonzero )  /* structural sparsity  pattern  */
{   
  unsigned char operation;
  int tape_stats[11];        /* tape stats */

  locint result=0;
  locint arg=0;
  locint res=0;
  locint loc1=0;
  locint loc2=0;
  double stored_val=0.0;
  double r0b,divs;
  locint result_v=0;
  locint loc1_v=0;
  locint loc2_v=0;
  double *d=0;
  locint size = 0;

  int i,l,j,k,k1,p,pd,pdk;
  static double** A;
  static revreal** Tr;
  static int kax, rax, pax;

  double x,y,r_0;
  double r0;
  int indexi;
  int indexd;
  int rev_location_cnt;
  int dep_cnt;
  int indep_cnt;
  int buffer;
  int taycheck;
  int numdep,numind;
  static revreal *A1, *A2, *Ares, *Tr1, *Tr2, *Tres;
  static revreal *Atemp, *Atemp2, *Trtemp;

  tag = tnum;   /* tag is global which indicates which tape to look at */
  pd = nrows;
  k = degre+1;
  k1 = k+1;
  pdk = pd*k1;

  tapestats(tag,tape_stats);

  indep_cnt = tape_stats[0];
  dep_cnt = tape_stats[1];
  rev_location_cnt = tape_stats[2];
  buffer =  tape_stats[4];

  set_buf_size(buffer);

  if ((depen != dep_cnt)||(indep != indep_cnt))
  {
    printf("ADOL-C error: reverse sweep on tape %d aborted!\n",tag);
    printf("Number of dependent and/or independent variables passed");
    printf(" to reverse is\ninconsistant with number ");
    printf("recorded on tape %d \n",tag);
    exit (-1);
  }


  indexi = indep_cnt - 1;
  indexd = dep_cnt - 1;

  if (k compsize kax || rev_location_cnt compsize  rax || pd compsize  pax)
    {
      if(pax)
      {
	/* delete Atemp;
	   delete Atemp2;
	   delete Trtemp; */
	free((char*) Atemp);
	free((char*) Atemp2);
	free((char*) Trtemp);
        free((char*) *Tr);
        free((char*) Tr);
        free((char*) *A);
	free((char*) A);
      }
      Atemp = (revreal *)malloc(sizeof(revreal)*k1);
      Atemp2 = (revreal *)malloc(sizeof(revreal)*k1);
      Trtemp = (revreal *)malloc(sizeof(revreal)*k);

    /* Atemp = new revreal[k1];
       Atemp2 = new revreal[k1];
       Trtemp = new revreal[k]; */
    Tr = myalloc2(rev_location_cnt,k);
    A = myalloc2(rev_location_cnt,pdk);
    kax = k;
    pax = pd;
    rax = rev_location_cnt;
  }
  taylor_back2(Tr,&numdep,&numind,&taycheck);

  if(taycheck != degre)   
    {
      printf("\n ADOL-C error: reverse fails because it was not preceeded \n");
      printf("by a forward sweep with degree of %d degree was %d !!!!!!\n",taycheck,degre);
      exit(-2);
    };

  if((numdep != depen)||(numind != indep))
    {
      printf("\n ADOL-C error: reverse fails on tape %d because the number of\n"
	     ,tag);
      printf("independent and/or dependent variables given to reverse are\n");
      printf("inconsistant with that of the internal taylor array.\n");
      exit(-2);
    }

  
  /* Set up the tape */

  /* my_init_rev_sweep(tag);*/
  init_rev_sweep(tag);

  /* Get the last operation */

  operation=get_op_r();
  
  while (operation != start_of_tape) 
    {

      /* Switch Statement to execute operands (in reverse) */

      switch (operation) {

	/* Markers */

      case end_of_int:
        get_loc_block_r(); /* Get the next int block */
        break;
      case end_of_val:
        get_val_block_r(); /* Get the next val block */
        break;
      case end_of_op:
        get_op_block_r();
        operation = get_op_r();/* Skip next operation, it's another end_of_op */
        break;
      case start_of_tape:
      case end_of_tape:
	break;

	/* Scalar Operations */

      case int_adb_a:
        loc1 = get_locint_r();
        loc2 = get_locint_r();
	/* get_taylors(loc1,k); */
	A1 = A[loc2];
	Ares = A[loc1];
	for (p=0;p<pd;p++)
	  {
	    if  (0==*Ares)
	      {
		Ares += k1;
		A1 += k1;
	      }
	    else 
	      {
		maxinc(*A1,*Ares);
		A1++; *Ares++ = 0;
		for (i=0;i<k;i++)
		  {
		    *A1++ += *Ares;
		    *Ares++ = 0 ;
                  }
	      }
	  }
	break;
      case assign_a:
        loc1 = get_locint_r();
        loc2 = get_locint_r();
	get_taylors(loc1,k);
	A1 = A[loc2];
	Ares = A[loc1];
	for (p=0;p<pd;p++)
	  {
	    if  (0==*Ares)
	      {
		Ares += k1;
		A1 += k1;
	      }
	    else 
	      {
		maxinc(*A1,*Ares);
		A1++; *Ares++ = 0;
		for (i=0;i<k;i++)
		  {
		    *A1++ += *Ares;
		    *Ares++ = 0 ;
                  }
	      }
	  }
	break;
	
      case assign_d:
        loc1 = get_locint_r();
        stored_val = get_val_r();
	Ares = A[loc1];
	for (p=0;p<pdk;p++)
	  *Ares++ = 0.0;
	get_taylors(loc1,k);
	break;
      case int_adb_d:
        loc1 = get_locint_r();
        stored_val = get_val_r();
	Ares = A[loc1];
	for (p=0;p<pdk;p++)
	  *Ares++ = 0.0;
	/* get_taylors(loc1,k); */
	break;
	
      case assign_ind:
        loc1 = get_locint_r();
	get_taylors(loc1,k);
	Ares = A[loc1];
	for (p=0;p<pd;p++)
	  {
	    if(nonzero)
	      nonzero[p][indexi] = (int)*Ares;
	    Ares++;
	    for (i=0;i<k;i++)
	      results[p][indexi][k-1-i]= *Ares++; 
	  }
	indexi--;
	break;
	
      case assign_dep:
        loc1 = get_locint_r();
	Ares = A[loc1];
	for (p=0;p<pd;p++)
	  {
	    Ares[k] = lagrange[p][indexd];
	    if(Ares[k]) *Ares = 1.0;
	    Ares += k1;
	  }
	indexd--;
	break;
	
      case eq_plus_d:
        result=get_locint_r();
        stored_val = get_val_r();
	get_taylors(result,k);
	break;
	
      case eq_plus_a: 
        result=get_locint_r();
        loc1=get_locint_r();
	get_taylors(result,k);
	Ares=A[result];
	A1=A[loc1];
	for (p=0;p<pd;p++)
	  if  (0 == *Ares)
	    {
	      Ares += k1;
	      A1 += k1;
	    }
	  else
	    {
	      maxinc(*A1,*Ares);
	      A1++; Ares++;
	      for (i=0;i<k;i++)
		*A1++ += *Ares++;
	    }
	break;
	
      case eq_min_d:
        result=get_locint_r();
        stored_val = get_val_r();
	get_taylors(result,k);
	break;
	
      case eq_min_a: 
        result=get_locint_r();
        loc1=get_locint_r();
	get_taylors(result,k);
	Ares=A[result];
	A1=A[loc1];
	for (p=0;p<pd;p++)
	  if  (0==*Ares)
	    {
	      Ares += k1;
	      A1 += k1;
	    }
	  else 
	    {
	      maxinc(*A1,*Ares);
	      A1++; Ares++;
	      for (i=0;i<k;i++)
		*A1++ -= *Ares++;
	    }
	break;
	
      case eq_mult_d:
        result=get_locint_r();
        stored_val = get_val_r();
	get_taylors(result,k);
	Ares = A[result];
	for (p=0;p<pd;p++)
	  if  (0==*Ares++)
	    Ares += k;
	  else for (i=0;i<k;i++)
	    *Ares++  *= stored_val;
	break;
	
      case eq_mult_a:
        result=get_locint_r();
        loc1=get_locint_r();
	get_taylors(result,k);
	Tr1 = Tr[loc1];
	Tres = Tr[result];
	A1 = A[loc1];
	Ares = A[result];
	for (p=0;p<pd;p++)
	  {
	    if (0==*Ares) 
	      {
		A1 += k1;
		Ares += k1;
	      }
	    else
	      {
		maxinc(*Ares,2.0);
		maxinc(*A1,*Ares);
		A1++; Ares++;
		for (i=k-1;i>=0;i--)
		  {
		    x = 0.0;
		    y = 0.0;
		    Atemp[i] = Ares[i];
		    Ares[i] = 0;
		    for (j=i;j<k;j++)
		      {
			x+=Atemp[j]*Tres[j-i];
			y+=Atemp[j]*Tr1[j-i];
		      }
		    A1[i] +=x;
		    Ares[i] +=y;
		  }
		A1 += k;
		Ares += k;
              }
	  }
	break;
	
      case plus_a_a:
        result =get_locint_r();
        loc2=get_locint_r();
        loc1=get_locint_r();
	Ares = A[result];
	A1 = A[loc1];
	A2 = A[loc2];
	for (p=0;p<pd;p++)
	  if  (0==*Ares)
	    {
	      Ares+=k1;
	      A1 += k1;
	      A2 += k1;
	    }
	  else
	    {
	      maxinc(*A1,*Ares);
	      maxinc(*A2,*Ares);
	      A2++; A1++; Ares++;
	      for (i=0;i<k;i++)
		{
		  *A1++ += *Ares;
		  *A2++ += *Ares++;
		}
	    }
	break;
	
      case plus_d_a:
        result = get_locint_r();
        loc1=get_locint_r();
        stored_val = get_val_r();
	Ares = A[result];
	A1 = A[loc1];
	for (p=0;p<pd;p++)
	  if  (0==*Ares)
	    {
	      Ares+=k1;
	      A1 += k1;
	    }
	  else
	    {
	      maxinc(*A1,*Ares);
	      A1++; Ares++;
	      for (i=0;i<k;i++)
	        *A1++ += *Ares++;
	    }
	break;
	
      case min_a_a:
        result =get_locint_r();
        loc2=get_locint_r();
        loc1=get_locint_r();
	Ares = A[result];
	A1 = A[loc1];
	A2 = A[loc2];
	for (p=0;p<pd;p++)
	  if  (0==*Ares)
	    {
	      Ares+=k1;
	      A1 += k1;
	      A2 += k1;
	    }
	  else
	    {
	      maxinc(*A1,*Ares);
	      maxinc(*A2,*Ares);
	      A2++; A1++; Ares++;
	      for (i=0;i<k;i++)
	        {
		  *A1++ += *Ares;
		  *A2++ -= *Ares++;
                }
	    }
	break;
	
      case min_d_a:
        result = get_locint_r();
        loc1=get_locint_r();
        stored_val = get_val_r();
	Ares = A[result];
	A1 = A[loc1];
	for (p=0;p<pd;p++)
	  if (0==*Ares)
	    {
	      Ares+=k1;
	      A1 += k1;
	    }
	  else
	    {
	      maxinc(*A1,*Ares);
	      A1++; Ares++;
	      for (i=0;i<k;i++)
		*A1++ -= *Ares++;
	    }
	break;
	
      case mult_a_a:
        result =get_locint_r();
        loc2=get_locint_r();
        loc1=get_locint_r();
	Tr1 = Tr[loc1];
	Tr2 = Tr[loc2];
	Ares = A[result];
	A2 = A[loc2];
	A1 = A[loc1];
	for (p=0;p<pd;p++)
	  if (0==*Ares)
	    { 
	      A1 += k1;
	      A2 += k1;
	      Ares += k1;
	    }
	  else
	    {
	      revreal comp = (*Ares > 2.0) ? *Ares : 2.0 ;
	      maxinc(*A1,comp);
	      maxinc(*A2,comp);
	      A1++; A2++; Ares++;
              for (i=0;i<k;i++)
                {
		  x = 0.0;
		  y = 0.0;
		  for (j=i;j<k;j++)
		    {
		      x+=Ares[j]*Tr1[j-i];
		      y+=Ares[j]*Tr2[j-i];
		    }
		  *A2++ +=x;
		  *A1++ +=y;
                }
              Ares += k;
	    } 
	break;
	
      case mult_d_a:
        result = get_locint_r();
        loc1=get_locint_r();
        stored_val = get_val_r();
	Ares = A[result];
	A1= A[loc1];
	for (p=0;p<pd;p++)
	  if  (0==*Ares)
	    {
	      Ares += k1;
	      A1 += k1;
	    }
	  else
	    {
              maxinc(*A1,*Ares);
	      A1++; Ares++;
	      for (i=0;i<k;i++)
		*A1++ += stored_val* (*Ares++);
	    }
	break;
	
      case div_a_a:
        result =get_locint_r();
        loc2=get_locint_r();
        loc1=get_locint_r();
	Tr2 = Tr[loc2];
	Tres = Tr[result];
	r0=1/(*Tr2);
	Ares = A[result];
	A2 = A[loc2];
	A1 = A[loc1];
	for (p=0;p<pd;p++)
	  {
	    if (0==*Ares)
	      {
		Ares += k1;
		A1 += k1;
		A2 += k1;
	      }
	    else
	      {
		maxinc(*A1,3.0);
		maxinc(*A1,*Ares);
		maxinc(*A2,3.0);
		maxinc(*A2,*Ares);
		A1++; A2++; Ares++;
		for (i=k-2;i>=0;i--)
		  {
		    x = 0.0;
		    for (j=i+1;j<k;j++)
		      x+=Ares[j]*Tr2[j-i];
		    Ares[i]-=x*r0;
		  }
		r_0 = 0.0;
		for (i=0;i<k;i++)
		  r_0+=Ares[i]*Tres[i];
		*A2++ -= r_0*r0;
		for (i=1;i<k;i++)
		  {
		    x = 0.0;
		    for (j=i;j<k;j++)
		      x+=Ares[j]*Tres[j-i];
		    *A2++ -=x*r0;
		  }
		for (i=0;i<k;i++)
		  *A1++ +=r0*Ares[i];
		Ares += k;
              }
	  }
	break;
	
      case div_d_a:
        result = get_locint_r();
        loc1=get_locint_r();
        stored_val = get_val_r();
	Tr1 = Tr[loc1];
	r0=1/(*Tr1);
	Tres = Tr[result];
	Ares = A[result];
	A1 = A[loc1];
            for (p=0;p<pd;p++)
              { 
	      if (0==*Ares)
		{
		Ares += k1;
		A1 += k1;
		}
             else
	      {
	      maxinc(*A1,*Ares);
	      maxinc(*A1,3.0);
	      A1++; Ares++;
              for (i=k-2;i>=0;i--)
                {
                x = 0.0;
                for (j=i+1;j<k;j++)
                   x += Ares[j]*Tr1[j-i];
                Ares[i] -= r0*x;
                }
              x = 0.0;
              for (i=0;i<k;i++)
                x+=Ares[i]*Tres[i];
              r_0=Tr1[0]*x;
              *A1++ -=r_0*r0*r0;
              for (i=1;i<k;i++)
                {
                x = 0.0;
                for (j=i;j<k;j++)
                    x+=Ares[j]*Tres[j-i];
                *A1++ -=r0*x;
                }
              Ares += k;
      	      }
              }
            break;

        case pow_op:
        result = get_locint_r();
        loc1=get_locint_r();
        stored_val = get_val_r();
	    Tr1 = Tr[loc1];
	    Tres = Tr[result];
	    Ares = A[result];
	    A1 = A[loc1];
            if (Tr1[0]==0.0)
              r0=0.0;
            else
	      r0 = 1/Tr1[0];
            for (p=0;p<pd;p++)
              {
              if (0==*Ares)
                {
                Ares += k1;
                A1 += k1;
                }
              else
                {
                if (Tr1[0]==0.0){
                *A1++ =5.0; Ares++;
                }
                else {
                maxinc(*A1,*Ares);
                maxinc(*A1,3.0);
                A1++; Ares++;
                }
                r0b = 0.0;
                for (i=k-1;i>0;i--)
                  {
                  x = Ares[i]*r0/i;
                  r0b += Ares[i]*Tres[i];
                  divs = -i;
                  for (j=i-1;j>=0;j--)
                    { divs += stored_val+1;
                       Ares[j] += x*divs*Tr1[i-j];
                       A1[i-j] += x*divs*Tres[j];
                    }
                  }
                *A1 += (*Ares*stored_val*Tres[0]-r0b)*r0;
                Ares += k;
                A1 += k;
               }
              }
          break; 

        case death_not:
        loc2 = get_locint_r();
        loc1 = get_locint_r();
            for (i=loc1;i<=loc2;i++)
              {
	       A1 = A[i];
	       for (p=0;p<pdk;p++)
                  *A1++ = 0.0;
               get_taylors(i,k);
              }
            break;

        case exp_op:
        result = get_locint_r();
        loc1=get_locint_r();
            Tr1 = Tr[loc1];
            Tres = Tr[result];
            Ares = A[result];
            A1 = A[loc1];
            for(i=1;i<k;i++)
               Trtemp[i] = i * Tr1[i];
            for (p=0;p<pd;p++)
              {
	      if (0==*Ares)
	      {
	       A1 += k1;
	       Ares += k1;
	       }
              else
	      {
               maxinc(*A1,*Ares);
               maxinc(*A1,4.0);
               A1++; Ares++;
              for (i=k-2;i>=0;i--)
                {
                Atemp[i+1] = Ares[i+1] / (i+1);
                x = 0.0;
                for (j=i+1;j<k;j++)
                    x+=Atemp[j]*Trtemp[j-i];
                Ares[i]+=x;
                }
              *A1++ += Ares[0]*Tres[0];
              for (i=1;i<k;i++)
                {
                x = 0.0;
                for (j=i;j<k;j++)
                   x+=Atemp[j]*Tres[j-i];
                *A1++ +=i*x;
                }
              Ares += k;
     	      }
              }
            break;
            
        case sin_op:
        result = get_locint_r();
        loc2 = get_locint_r();
        loc1 = get_locint_r();
        stored_val = get_val_r();
            Tr1 = Tr[loc1];
            Tr2 = Tr[loc2];
            Tres = Tr[result];
            Ares = A[result];
            A2 = A[loc2];
            A1 = A[loc1];
            for(i=1;i<k;i++)
                Trtemp[i] = i * Tr1[i];
            for (p=0;p<pd;p++)
              {
	      if (0==*Ares)
	       {
		A1 += k1;
		A2 += k1;
		Ares += k1;
	       }
              else
	      {
              maxinc(*A1,*Ares);
              maxinc(*A1,4.0);
              maxinc(*A2,*Ares);
              maxinc(*A2,4.0);
              A1++; A2++; Ares++;
              for (i=k-2;i>=0;i--)
                {
                x = 0.0;
                y = 0.0;
                Atemp[i+1] = Ares[i+1] / (i+1);
                Atemp2[i+1] = A2[i+1] / (i+1);
                for (j=i+1;j<k;j++)
                  {
                  x+=Atemp[j]*Trtemp[j-i];
                  y+=Atemp2[j]*Trtemp[j-i];
                  }
                A2[i]+=x;
                Ares[i]-=y;
                }
              *A1++ += Ares[0]*Tr2[0]-A2[0]*Tres[0];
              for (i=1;i<k;i++)
                {
                x = 0.0;
                for (j=i;j<k;j++)
                   x += Tr2[j-i]*Atemp[j]-Atemp2[j]*Tres[j-i] ;
                *A1++ +=i*x;
                }
              A2 += k;
	      Ares += k;
              }
              }
            break;
            
        case cos_op:
        result = get_locint_r();
        loc2 = get_locint_r();
        loc1 = get_locint_r();
        stored_val = get_val_r();
            Tr1 = Tr[loc1];
            Tres = Tr[result];
            Tr2 = Tr[loc2];
            Ares = A[result];
            A2 = A[loc2];
            A1 = A[loc1];
            for(i=1;i<k;i++)
               Trtemp[i] = i * Tr1[i];
            for (p=0;p<pd;p++)
              {
	      if (0==*Ares)
	      {
	      A1 += k1;
	      A2 += k1;
	      Ares += k1;
	      }
	      else
	      {
              maxinc(*A1,*Ares);
              maxinc(*A1,4.0);
              maxinc(*A2,*Ares);
              maxinc(*A2,4.0);
              A1++; A2++; Ares++;
              for (i=k-2;i>=0;i--)
                {
                x = 0.0;
                y = 0.0;
                Atemp[i+1] = Ares[i+1] / (i+1);
                Atemp2[i+1] = A2[i+1] / (i+1);
                for (j=i+1;j<k;j++)
                  {
                  x+=Atemp2[j]*Trtemp[j-i];
                  y+=Atemp[j]*Trtemp[j-i];
                  }
                Ares[i]+=x;
                A2[i]-=y;
                }
              *A1++ += A2[0]*Tres[0]-Ares[0]*Tr2[0];
              for (i=1;i<k;i++)
              {
                x = 0.0;
                for (j=i;j<k;j++)
                   x += Tres[j-i]*Atemp2[j]-Atemp[j]*Tr2[j-i] ;
                *A1++ +=i*x;
              }
	      A2 += k;
	      Ares += k;
	      }
	    } 
            break;
            
        case sqrt_op:
        result = get_locint_r();
        loc1=get_locint_r();
            Ares = A[result];
            A1 = A[loc1];
            Tres = Tr[result];
            if ((*Tres)==0.0)
              r0=0.0;
            else
              r0=1.0/(*Tres);
            for (p=0;p<pd;p++)
              {
		if (0==*Ares)
		  {
		    A1 += k1;
		    Ares += k1;
		  }
		else
		  {
                    if ((*Tres)==0.0) {
                    *A1++ = 5.0; Ares++;
                    } 
                    else {
                    maxinc(*A1,*Ares);
                    maxinc(*A1,4.0);
                    A1++; Ares++;
                    }
		    for (i=k-2;i>0;i--)
		      {
			x = 0.0;
			for (j=i+1;j<k;j++)
			  x += Ares[j]*Tres[j-i];
			Ares[i]-=r0*x;
		      }
		    x = 0.0;
		    r_0 = r0/2;
		    for (i=1;i<k;i++)
		      {
			x+=Ares[i]*Tres[i];
			A1[i] += r_0*Ares[i];
		      }
		    Ares[0] -= x*r0;
		    A1[0] += r_0*Ares[0];
		    Ares +=k;
		    A1 += k;
		  }
              }
            break;
        case abs_val:
        result = get_locint_r();
        loc1=get_locint_r();
            Ares = A[result];
            A1 = A[loc1];
            Tr1 = Tr[loc1];
            x=0.0;
            j=0;
            for (i=0;i<k;i++) {
              if((x==0.0) && (Tr1[i]!=0.0)) {
                j=i;
                if(Tr1[i]<0.0)
                  x = -1.0;
                else
                  x= 1.0;
              } /* end if */
            } /* end for */
            for (p=0;p<pd;p++) {
              if (0==*Ares) {
                A1 += k1;
                Ares += k1;
              } /* endif */
              else {
                if(Tr1[0]==0.0) {
                  *A1++=5.0; Ares++;
                } /* end if */
                else {
                maxinc(*A1,*Ares);
                A1++; Ares++;
                } /* end else */
                for (i=j;i<k;i++)
                  A1[i] += x*Ares[i];
                Ares +=k;
                A1 += k; 
              } /* end else */
            } /* endfor */
            break;


	    
        case asin_op:
        case acos_op:
        case atan_op:
        case gen_quad:
        result = get_locint_r();
        loc2 = get_locint_r();
        loc1 = get_locint_r();
        stored_val = get_val_r();
            Ares = A[result];
            A2 = A[loc2];
            A1 = A[loc1];
            Tr1 = Tr[loc1];
            Tr2 = Tr[loc2];
            for(i=1;i<k;i++)
	      Trtemp[i] = i*Tr1[i];
            for (p=0;p<pd;p++)
	      {
		if(0==*Ares)
		  {
		    A1 += k1;
		    A2 += k1;
		    Ares += k1;
		  }
		else
		  {
                    maxinc(*A1,*Ares);
                    maxinc(*A1,4.0);
                    maxinc(*A2,*Ares);
                    maxinc(*A2,4.0);
                    A1++; A2++; Ares++;
		    for (i=k-1;i>0;i--)
		      {
			Atemp[i] = Ares[i] / i;
			x = 0.0;
			for (j=i;j<k;j++)
			  x+=Atemp[j]*Trtemp[j-i+1];
			A2[i-1]+=x;
		      }
		    *A1++ +=Tr2[0]*Ares[0];
		    for (i=1;i<k;i++)
		      {
			x = 0.0;
			for (j=i;j<k;j++)
			  x+=Atemp[j]*Tr2[j-i];
			*A1++ +=i*x;
		      }
		    A2 += k; 
		    Ares += k;
		  }
	      }
            break;      
            
        case log_op:
        result = get_locint_r();
        loc1=get_locint_r();
            Ares = A[result];
            A1 = A[loc1];
            Tr1 = Tr[loc1];
            Tres = Tr[result];
            r0=1/Tr1[0];
            for(i=1;i<k;i++)
	      Trtemp[i] = i*Tres[i];
            for (p=0;p<pd;p++)
              {
		if (0==*Ares)
		  {
		    A1 += k1;
		    Ares += k1;
		  }
		else
		  { 
                    maxinc(*A1,*Ares);
                    maxinc(*A1,4.0);
                    A1++; Ares++;
		    for (i=k-2;i>0;i--)
		      {
			x = 0.0;
			Atemp[i+1] = Ares[i+1] / (i+1);
			for (j=i+1;j<k;j++)
			  x+=Atemp[j]*Tr1[j-i];
			Ares[i]-=x*i*r0;
		      }
		    x = 0.0;
		    for (i=1;i<k;i++)
		      x += Ares[i]*Tres[i];
		    r_0 = x/r0;
		    *A1++ +=r0*Ares[0]-r0*r0*r_0;
		    for (i=1;i<k;i++)
		      {
			x=Ares[i];
			for (j=i+1;j<k;j++)
			  x-=Atemp[j]*Trtemp[j-i];
			*A1++ +=r0*x;
		      }
		    Ares += k;
		  }
              }
            break;

	/* Vector Operations */

      case int_av_av:
        result_v= get_locint_r();
        size = get_locint_r();
        loc1_v = get_locint_r();
	for (l=0;l<size;l++)
	  {
	    loc1 = result_v + l;   /* Location of left-hand-side  */
	    loc2 = loc1_v + l;   /* Location of right-hand-side */
	    /* code for int_adb_a */
	    A1 = A[loc2];
	    Ares = A[loc1];
	    for (p=0;p<pd;p++)
	      {
		if  (0==*Ares)
		  {
		    Ares += k1;
		    A1 += k1;
		  }
		else 
		  {
		    maxinc(*A1,*Ares);
		    A1++; *Ares++ = 0;
		    for (i=0;i<k;i++)
		      {
			*A1++ += *Ares;
			*Ares++ = 0 ;
		      }
		  }
	      }
	    
	  }
	break;
      case assign_dv:
        loc1_v = get_locint_r();
        size = get_locint_r();
        d = get_val_v_r(size);
	for (l=size-1;l>=0;l--)
	  {
	    loc1 = loc1_v + l;           /* Location of left-hand-side */
	    stored_val = d[l];      /* Value of right-hand-side   */     
	    /* code for assign_d */
	    Ares = A[loc1];
	    for (p=0;p<pdk;p++)
	      *Ares++ = 0.0;
	    get_taylors(loc1,k);
	  }
	reset_val_r();
	break;
      case assign_av:
        loc1_v = get_locint_r();
        size = get_locint_r();
        loc2_v = get_locint_r();
	for (l=size-1;l>=0;l--)
	  {
	    loc1 = loc1_v + l;             /* Location of left-hand-side  */
	    loc2 = loc2_v + l;             /* Location of right-hand-side */
	    /* code for assign_a */
	    get_taylors(loc1,k);
	    A1 = A[loc2];
	    Ares = A[loc1];
	    for (p=0;p<pd;p++)
	      {
		if  (0==*Ares)
		  {
		    Ares += k1;
		    A1 += k1;
		  }
		else 
		  {
		    maxinc(*A1,*Ares);
		    A1++; *Ares++ = 0;
		    for (i=0;i<k;i++)
		      {
			*A1++ += *Ares;
			*Ares++ = 0 ;
		      }
		  }
	      }
	  }
	break;
      case assign_indvec:
        loc1_v = get_locint_r();
        size = get_locint_r();
	for (l=size-1;l>=0;l--)
	  {
	    loc1 = loc1_v + l;             /* Location of the left-hand-side */
	    /* code for assign_ind */
	    get_taylors(loc1,k);
	    Ares = A[loc1];
	    for (p=0;p<pd;p++)
	      {
		if(nonzero)
		  nonzero[p][indexi] = (int)*Ares;
		Ares++;
		for (i=0;i<k;i++)
		  results[p][indexi][k-1-i]= *Ares++; 
	      }
	    indexi--;
	  }
	reset_val_r();
	break;
      case assign_depvec:
        loc1_v = get_locint_r();
        size = get_locint_r();
	for (l=size-1;l>=0;l--)
	  {
	    loc1 = loc1_v + l;             /* Location of the left-hand-side */
	    /* code for assign_dep */
	    Ares = A[loc1];
	    for (p=0;p<pd;p++)
	      {
		Ares[k] = lagrange[p][indexd];
		if(Ares[k]) *Ares = 1.0;
		Ares += k1;
	      }
	    indexd--;
	  }
	break;
      case eq_min_av:
        result_v = get_locint_r();
        size = get_locint_r();
        loc1_v = get_locint_r();
	for (l=size-1;l>=0;l--)
	  {
	    result = result_v + l;            /* Location of left-hand-side  */
	    loc1   = loc1_v + l;              /* Location on right-hand-side */
	    /* code for eq_min_a */ 
	    get_taylors(result,k);
	    Ares=A[result];
	    A1=A[loc1];
	    for (p=0;p<pd;p++)
	      if  (0==*Ares)
		{
		  Ares += k1;
		  A1 += k1;
		}
	      else 
		{
		  maxinc(*A1,*Ares);
		  A1++; Ares++;
		  for (i=0;i<k;i++)
		    *A1++ -= *Ares++;
		}
	  }
	break;
      case eq_plus_av:
        result_v = get_locint_r();
        size = get_locint_r();
        loc1_v = get_locint_r();
	for (l=size-1;l>=0;l--)
	  {
	    result = result_v + l;            /* Location of left-hand-side  */
	    loc1   = loc1_v   + l;            /* Location on right-hand-side */
	    /* code for eq_plus_a */
	    get_taylors(result,k);
	    Ares=A[result];
	    A1=A[loc1];
	    for (p=0;p<pd;p++)
	      if  (0 == *Ares)
		{
		  Ares += k1;
		  A1 += k1;
		}
	      else
		{
		  maxinc(*A1,*Ares);
		  A1++; Ares++;
		  for (i=0;i<k;i++)
		    *A1++ += *Ares++;
		}
	    
	  }
	break;
      case eq_mult_av_a:
        result_v = get_locint_r();
        size = get_locint_r();
        loc1 = get_locint_r();
	for (l=size-1;l>=0;l--)
	  {
	    result = result_v + l;      /* Location of the left-hand-side  */
	    /* loc1   = fixed;             Location on the right-hand-side */
	    
	    /* code for eq_mult_a*/
	    get_taylors(result,k);
	    Tr1 = Tr[loc1];
	    Tres = Tr[result];
	    A1 = A[loc1];
	    Ares = A[result];
	    for (p=0;p<pd;p++)
	      {
		if (0==*Ares) 
		  {
		    A1 += k1;
		    Ares += k1;
		  }
		else
		  {
		    maxinc(*Ares,2.0);
		    maxinc(*A1,*Ares);
		    A1++; Ares++;
		    for (i=k-1;i>=0;i--)
		      {
			x = 0.0;
			y = 0.0;
			Atemp[i] = Ares[i];
			Ares[i] = 0;
			for (j=i;j<k;j++)
			  {
			    x+=Atemp[j]*Tres[j-i];
			    y+=Atemp[j]*Tr1[j-i];
			  }
			A1[i] +=x;
			Ares[i] +=y;
		      }
		    A1 += k;
		    Ares += k;
		  }
	      }
	  }
	break;
      case eq_mult_av_d:
        result_v = get_locint_r();
        size = get_locint_r();
        stored_val = get_val_r();
	for (l=size-1;l>=0;l--)
	  {
	    result = result_v + l;      /* Location of the left-hand-side  */
	    /* stored_val = fixed;            value on the right-hand-side */
	    /* code for eq_mult_d*/
	    get_taylors(result,k);
	    Ares = A[result];
	    for (p=0;p<pd;p++)
	      if  (0==*Ares++)
		Ares += k;
	      else for (i=0;i<k;i++)
		*Ares++  *= stored_val;
	  }
	break;
      case sub_av_av:
        result_v = get_locint_r();
        size     = get_locint_r();
        loc2_v   = get_locint_r();
        loc1_v   = get_locint_r();
	for (l=0;l<size;l++)
	  {
	    loc2   = loc2_v    + l;       /* Location of var 2  */
	    loc1   = loc1_v    + l;       /* Location of var 1  */
	    result = result_v  + l;       /* Location of result */
	    /* code for min_a_a */
	    Ares = A[result];
	    A1 = A[loc1];
	    A2 = A[loc2];
	    for (p=0;p<pd;p++)
	      if  (0==*Ares)
		{
		  Ares+=k1;
		  A1 += k1;
		  A2 += k1;
		}
	      else
		{
		  maxinc(*A1,*Ares);
		  maxinc(*A2,*Ares);
		  A2++; A1++; Ares++;
		  for (i=0;i<k;i++)
		    {
		      *A1++ += *Ares;
		      *A2++ -= *Ares++;
		    }
		}
	    
	  }
	break;
      case plus_av_av:
        result_v = get_locint_r();
        size     = get_locint_r();
        loc2_v   = get_locint_r();
        loc1_v   = get_locint_r();
	for (l=0;l<size;l++)
	  {
	    loc2   = loc2_v    + l;       /* Location of var 2  */
	    loc1   = loc1_v    + l;       /* Location of var 1  */
	    result = result_v  + l;       /* Location of result */
	    /* code for plus_a_a */
	    Ares = A[result];
	    A1 = A[loc1];
	    A2 = A[loc2];
	    for (p=0;p<pd;p++)
	      if  (0==*Ares)
		{
		  Ares+=k1;
		  A1 += k1;
		  A2 += k1;
		}
	      else
		{
		  maxinc(*A1,*Ares);
		  maxinc(*A2,*Ares);
		  A2++; A1++; Ares++;
		  for (i=0;i<k;i++)
		    {
		      *A1++ += *Ares;
		      *A2++ += *Ares++;
		    }
		}
	    
	  }
	break;
      case dot_av_av:
        result_v = get_locint_r();
        size     = get_locint_r();
        loc2_v   = get_locint_r();
        loc1_v   = get_locint_r();
	result = result_v;
	for (l=0;l<size;l++)
	  {
	    loc2 = loc2_v + l;
	    loc1 = loc1_v + l;
	    /* code for mult_a_a  */
	    Tr1 = Tr[loc1];
	    Tr2 = Tr[loc2];
	    Ares = A[result];
	    A2 = A[loc2];
	    A1 = A[loc1];
	    for (p=0;p<pd;p++)
	      if (0==*Ares)
		{ 
		  A1 += k1;
		  A2 += k1;
		  Ares += k1;
		}
	      else
		{
		  revreal comp = (*Ares > 2.0) ? *Ares : 2.0 ;
		  maxinc(*A1,comp);
		  maxinc(*A2,comp);
		  A1++; A2++; Ares++;
		  for (i=0;i<k;i++)
		    {
		      x = 0.0;
		      y = 0.0;
		      for (j=i;j<k;j++)
			{
			  x+=Ares[j]*Tr1[j-i];
			  y+=Ares[j]*Tr2[j-i];
			}
		      *A2++ +=x;
		      *A1++ +=y;
		    }
		  Ares += k;
		} 
	  }
	Ares = A[result];
	for (p=0;p<pdk;p++)
	  *Ares++ = 0.0;
	
	break;
      case mult_d_av:
        result_v = get_locint_r();
        size = get_locint_r();
        loc1_v = get_locint_r();
        stored_val =get_val_r();
	for (l=0;l<size;l++)
	  {
	    loc1   = loc1_v   + l;   /* Location on the right-hand-side */
	    result = result_v + l;   /* location of the result */
	    /* stored_val = Fixed double value */
	    
	    /* code for mult_d_a */
	    Ares = A[result];
	    A1= A[loc1];
	    for (p=0;p<pd;p++)
	      if  (0==*Ares)
		{
		  Ares += k1;
		  A1 += k1;
		}
	      else
		{
		  maxinc(*A1,*Ares);
		  A1++; Ares++;
		  for (i=0;i<k;i++)
		    *A1++ += stored_val* (*Ares++);
		}
	  }
	break;
      case div_av_a:
        result_v = get_locint_r();
        size     = get_locint_r();
        loc2_v   = get_locint_r();
        loc1_v   = get_locint_r();
	loc2= loc2_v;     /* Location of adouble */
	for (l=0;l<size;l++)
	  {
	    loc1   = loc1_v   + l;  /* Location of right-hand-side vector[l] */
	    result = result_v + l;  /* Location of the result */
	    
	    /* code for div_a_a */
	    Tr2 = Tr[loc2];
	    Tres = Tr[result];
	    r0=1/(*Tr2);
	    Ares = A[result];
	    A2 = A[loc2];
	    A1 = A[loc1];
	    for (p=0;p<pd;p++)
	      {
		if (0==*Ares)
		  {
		    Ares += k1;
		    A1 += k1;
		    A2 += k1;
		  }
		else
		  {
		    maxinc(*A1,3.0);
		    maxinc(*A1,*Ares);
		    maxinc(*A2,3.0);
		    maxinc(*A2,*Ares);
		    A1++; A2++; Ares++;
		    for (i=k-2;i>=0;i--)
		      {
			x = 0.0;
			for (j=i+1;j<k;j++)
			  x+=Ares[j]*Tr2[j-i];
			Ares[i]-=x*r0;
		      }
		    r_0 = 0.0;
		    for (i=0;i<k;i++)
		      r_0+=Ares[i]*Tres[i];
		    *A2++ -= r_0*r0;
		    for (i=1;i<k;i++)
		      {
			x = 0.0;
			for (j=i;j<k;j++)
			  x+=Ares[j]*Tres[j-i];
			*A2++ -=x*r0;
		      }
		    for (i=0;i<k;i++)
		      *A1++ +=r0*Ares[i];
		    Ares += k;
		  }
	      }
	  }
	break;
      case mult_av_a:
        result_v = get_locint_r();
        size     = get_locint_r();
        loc2_v   = get_locint_r();
        loc1_v   = get_locint_r();
	loc2= loc2_v;   /* Location of adouble */
	for (l=0;l<size;l++)
	  {
	    loc1   = loc1_v   + l; /* Location of right-hand-side vector[l] */
	    result = result_v + l; /* Location of the result */
	    
	    /* code for mult_a_a */
	    Tr1 = Tr[loc1];
	    Tr2 = Tr[loc2];
	    Ares = A[result];
	    A2 = A[loc2];
	    A1 = A[loc1];
	    for (p=0;p<pd;p++)
	      if (0==*Ares)
		{ 
		  A1 += k1;
		  A2 += k1;
		  Ares += k1;
		}
	      else
		{
		  revreal comp = (*Ares > 2.0) ? *Ares : 2.0 ;
		  maxinc(*A1,comp);
		  maxinc(*A2,comp);
		  A1++; A2++; Ares++;
		  for (i=0;i<k;i++)
		    {
		      x = 0.0;
		      y = 0.0;
		      for (j=i;j<k;j++)
			{
			  x+=Ares[j]*Tr1[j-i];
			  y+=Ares[j]*Tr2[j-i];
			}
		      *A2++ +=x;
		      *A1++ +=y;
		    }
		  Ares += k;
		} 
	    
	  }
	break;
      case mult_a_av:
        result_v = get_locint_r();
        size     = get_locint_r();
        loc2_v   = get_locint_r();
        loc1_v   = get_locint_r();
	loc1= loc2_v;      /* Location of the adouble */
	for (l=0;l<size;l++)
	  {
	    loc2= loc1_v+l;        /* Location of right hand side vector[l]  */
	    result = result_v + l; /* Location of the result */
	    
	    /* code for mult_a_a */
	    Tr1 = Tr[loc1];
	    Tr2 = Tr[loc2];
	    Ares = A[result];
	    A2 = A[loc2];
	    A1 = A[loc1];
	    for (p=0;p<pd;p++)
	      if (0==*Ares)
		{ 
		  A1 += k1;
		  A2 += k1;
		  Ares += k1;
		}
	      else
		{
		  revreal comp = (*Ares > 2.0) ? *Ares : 2.0 ;
		  maxinc(*A1,comp);
		  maxinc(*A2,comp);
		  A1++; A2++; Ares++;
		  for (i=0;i<k;i++)
		    {
		      x = 0.0;
		      y = 0.0;
		      for (j=i;j<k;j++)
			{
			  x+=Ares[j]*Tr1[j-i];
			  y+=Ares[j]*Tr2[j-i];
			}
		      *A2++ +=x;
		      *A1++ +=y;
		    }
		  Ares += k;
		} 
	  }
	break;
#ifdef conditional
      case cond_assign:
        result = get_locint_r();
        loc2 = get_locint_r();
        loc1 = get_locint_r();
        loc1_v = get_locint_r(); 
	get_taylors(result,k);
	A1 = A[loc1];
	Ares = A[result];
        A2 = A[loc2];
        Tr1 = Tr[loc1_v];
        if (*Tr1>0)
        {
	  for (p=0;p<pd;p++)
	  {
	    if  (0==*Ares)
	    {
              Ares += k1;
	      A1 += k1;
	    } /* endif */
	    else 
	    {
              maxinc(*A1,*Ares);
              A1++; *Ares++ = 0;
              for (i=0;i<k;i++)
              {
                 *A1++ += *Ares;
                 *Ares++ = 0 ;
              } /* endfor */
            } /* endelse */
	  } /* endfor */
        } /* endif */
        else
        {
          for (p=0;p<pd;p++)
          {
            if  (0==*Ares)
            {
              Ares += k1;
              A2 += k1;
            } /* endif */
            else
            {
              maxinc(*A2,*Ares);
              A2++; *Ares++ = 0;
              for (i=0;i<k;i++)
              {
                 *A2++ += *Ares;
                 *Ares++ = 0 ;
              } /* endfor */
            } /* endelse */
          } /* endfor */
        } /* endelse */
	break;
      case cond_assign_s:
        result = get_locint_r();
        loc1 = get_locint_r();
        loc1_v = get_locint_r(); 
	get_taylors(result,k);
	A1 = A[loc1];
	Ares = A[result];
        Tr1 = Tr[loc1_v];
        if (*Tr1>0)
        {
          for (p=0;p<pd;p++)
          {
            if  (0==*Ares)
            {
              Ares += k1;
              A1 += k1;
            } /* endif */
            else
            {
              maxinc(*A1,*Ares);
              A1++; *Ares++ = 0;
              for (i=0;i<k;i++)
              {
                 *A1++ += *Ares;
                 *Ares++ = 0 ;
              } /* endfor */
            } /* endelse */
          } /* endfor */
        } /* endif */
	break;
      case subscript:
        result = get_locint_r();
        loc1 = get_locint_r();
        loc1_v = get_locint_r(); 
	get_taylors(result,k);
	A1 = A[loc1_v+(int)(Tr[loc1][0])];
	Ares = A[result];
	for (p=0;p<pd;p++)
	  {
	    if  (0==*Ares)
	      {
		Ares += k1;
		A1 += k1;
	      }
	    else 
	      {
		maxinc(*A1,*Ares);
		A1++; *Ares++ = 0;
		for (i=0;i<k;i++)
		  {
		    *A1++ += *Ares;
		    *Ares++ = 0 ;
                  }
	      }
	  }
	break;
      case subscript_l:
        result = get_locint_r();
        loc1 = get_locint_r();
        loc1_v = get_locint_r(); 
	get_taylors(loc1_v+(int)(Tr[loc1][0]),k);
	A1 = A[result];
	Ares = A[loc1_v+(int)(Tr[loc1][0])];
	for (p=0;p<pd;p++)
	  {
	    if  (0==*Ares)
	      {
		Ares += k1;
		A1 += k1;
	      }
	    else 
	      {
		maxinc(*A1,*Ares);
		A1++; *Ares++ = 0;
		for (i=0;i<k;i++)
		  {
		    *A1++ += *Ares;
		    *Ares++ = 0 ;
                  }
	      }
	  }
	break;
      case subscript_ld:
        loc1 = get_locint_r();
        loc1_v = get_locint_r(); 
        stored_val = get_val_r();
	Ares = A[loc1_v+(int)(Tr[loc1][0])];
	for (p=0;p<pdk;p++)
	  *Ares++ = 0.0;
	break;
      case m_subscript:
        result = get_locint_r();
        size = get_locint_r();
        loc1 = get_locint_r();
        loc1_v = get_locint_r(); 
        for (l=size-1;l>=0;l--)
        {
          res=result+l;
          arg=loc1_v+(int)(Tr[loc1][0])*size+l;
	  A1 = A[arg];
          Ares = A[res];
	  for (p=0;p<pd;p++)
	  {
	    if  (0==*Ares)
	      {
		Ares += k1;
		A1 += k1;
	      }
	    else 
	      {
		maxinc(*A1,*Ares);
		A1++; *Ares++ = 0;
		for (i=0;i<k;i++)
		  {
		    *A1++ += *Ares;
		    *Ares++ = 0 ;
                  }
	      }
	  }
        } /* endfor */
	break;
      case m_subscript_l:
        result = get_locint_r();
        size = get_locint_r();
        loc1 = get_locint_r();
        loc1_v = get_locint_r(); 
        for (l=size-1;l>=0;l--)
        {
          res=result+l;
          arg= loc1_v+(int)(Tr[loc1][0])*size+l;
	  get_taylors(arg,k);
	  A1 = A[res];
	  Ares = A[arg];
	  for (p=0;p<pd;p++)
	  {
	    if  (0==*Ares)
	      {
		Ares += k1;
		A1 += k1;
	      }
	    else 
	      {
		maxinc(*A1,*Ares);
		A1++; *Ares++ = 0;
		for (i=0;i<k;i++)
		  {
		    *A1++ += *Ares;
		    *Ares++ = 0 ;
                  }
	      }
	  }
        } /* endfor */
	break;
      case m_subscript_ld:
        size = get_locint_r();
        loc2 = get_locint_r();
        loc1 = get_locint_r();
        loc1_v = get_locint_r(); 
        d = get_val_v_r(size);
        for (l=size-1;l>=0;l--)
        {
          arg=loc1_v+(int)(Tr[loc1][0])*size+l+loc2;
	  Ares = A[arg];
	  for (p=0;p<pdk;p++)
	    *Ares++ = 0.0;
        } /* endfor */
	break;
#endif	 
   
      default:
	printf("ADOL-C error: Fatal error in hov_reverse on operation %d\n",
	       operation);
	exit(-1);
	break;
      }

      /* Get the next operation */

      operation=get_op_r();
    }
  end_sweep();
}

#ifdef __cplusplus
}
#endif

