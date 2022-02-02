/*
   ----------------------------------------------------------------
   File fos_reverse.c of ADOL-C version 1.6 as of January 1,   1995
   ----------------------------------------------------------------
   Contains the routine fos_reverse (first-order-scalar reverse 
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

static short tag;


/****************************************************************************/
/* First-order scalar reverse.                                              */
/****************************************************************************/
void fos_reverse(short tnum,         /* tape id */
		 int depen,          /* consistency chk on # of dependents */
		 int indep,          /* consistency chk on # of independents */
		 double* lagrange,
		 double* results)    /*  coefficient vectors */
{
  unsigned char operation;
  int tape_stats[11];  /* tape stats */

  locint result=0;
  locint arg=0;
  locint res=0;
  locint loc1=0;
  locint loc2=0;
  double stored_val=0.0;
  locint result_v=0;
  locint loc1_v=0;
  locint loc2_v=0;
  double *d = 0;
  locint size = 0;

  static int rax;
  static double* As;
  static revreal* Trs;
  int i,l,kl;
  double r0;
  double* Ad;

  int indexi;
  int indexd;
  int rev_location_cnt;
  int dep_cnt;
  int indep_cnt;

  int buffer;

  int taycheck;
  int numdep,numind;

  tag = tnum;   /*tag is global which indicates which tape to look at */
    

    
  tapestats(tag,tape_stats);

  indep_cnt = tape_stats[0];
  dep_cnt = tape_stats[1];
  rev_location_cnt = tape_stats[2];
  buffer =  tape_stats[4];


  
  set_buf_size(buffer);

  if ((depen != dep_cnt)||(indep != indep_cnt))
    {
      printf("ADOL-C error: Reverse sweep on tape %d  aborted!\n",tag);
      printf("Number of dependent and/or independent variables ");
      printf("passed to reverse is\ninconsistant with number ");
      printf("recorded on tape %d \n",tag);
      exit (-1);
    }
  
  indexi = indep_cnt - 1;
  indexd = dep_cnt - 1;
/*
  for (kl =0; kl < indep_cnt; kl++) results[kl]=NaN;
*/
  
  if(rev_location_cnt compsize rax) 
    {
      if(rax)
	{ 
	  free((char*) Trs);
	  free((char*) As);
	}
      Trs = (revreal*) malloc(rev_location_cnt*sizeof(revreal));
      As = (revreal*) malloc(rev_location_cnt*sizeof(revreal));
      rax = rev_location_cnt;
    }
  Ad = As;
  for(i=0;i<rev_location_cnt;i++) *Ad++ = 0;
  taylor_back(Trs,&numdep,&numind,&taycheck);
  
  if(taycheck < 0)   
    {
      printf("\n ADOL-C error: reverse fails because it was not preceeded \n");
      printf("by a forward sweep with degree 0 \n");
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

  init_rev_sweep(tag);


  /* Read the last operation */

  operation=get_op_r();
  
  while (operation != start_of_tape) 
    {

      /* Switch Statemnet to execute operations in reverse */
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
        operation = get_op_r(); /* Skip next operation,it's another end_of_op */
        break;
      case start_of_tape:
      case end_of_tape:
	break;

	/* Scalar Operations */
      case int_adb_a:
        loc1 = get_locint_r();
        loc2 = get_locint_r();
	As[loc2] += As[loc1];
	As[loc1] = 0.0;
	break;
      case assign_a:
        loc1 = get_locint_r();
        loc2 = get_locint_r();
	get_taylor(loc1);
	As[loc2] += As[loc1];
	As[loc1] = 0.0;
	break;
      case int_adb_d:
        loc1 = get_locint_r();
        stored_val = get_val_r();
	As[loc1] = 0;
	break;
      case assign_d:
        loc1 = get_locint_r();
        stored_val = get_val_r();
	As[loc1] = 0;
	get_taylor(loc1);
	break;
      case assign_ind:
        loc1 = get_locint_r();
	get_taylor(loc1);
	results[indexi] =As[loc1];
	indexi--;
	break;
      case assign_dep:
        loc1 = get_locint_r();
	As[loc1] =lagrange[indexd];
	indexd--;
	break;
      case eq_plus_d:
        result=get_locint_r();
        stored_val = get_val_r();
	get_taylor(result);
	break;
      case eq_plus_a: 
        result=get_locint_r();
        loc1=get_locint_r();
	get_taylor(result);
	As[loc1] +=As[result];
	break;
      case eq_min_d:
        result=get_locint_r();
        stored_val = get_val_r();
	get_taylor(result);
	break;
      case eq_min_a: 
        result=get_locint_r();
        loc1=get_locint_r();
	get_taylor(result);
	As[loc1] -= As[result];
	break;
      case eq_mult_d:
        result=get_locint_r();
        stored_val = get_val_r();
	get_taylor(result);
	As[result] *=stored_val;
	break;
      case eq_mult_a:
        result=get_locint_r();
        loc1=get_locint_r();
	get_taylor(result);
	r0 = As[result];
	As[result] = 0;
	As[loc1] += r0*Trs[result];
	As[result] += r0*Trs[loc1] ; 
	break;
      case plus_a_a:
        result =get_locint_r();
        loc2=get_locint_r();
        loc1=get_locint_r();
	As[loc1]+=As[result];
	As[loc2]+=As[result];
	break;
      case plus_d_a:
        result = get_locint_r();
        loc1=get_locint_r();
        stored_val = get_val_r();
	As[loc1] += As[result];
	break;
      case min_a_a:
        result =get_locint_r();
        loc2=get_locint_r();
        loc1=get_locint_r();
	As[loc1]+=As[result];
	As[loc2]-=As[result];
	break;
      case min_d_a:
        result = get_locint_r();
        loc1=get_locint_r();
        stored_val = get_val_r();
	As[loc1] -= As[result];
	break;
      case mult_a_a:
        result =get_locint_r();
        loc2=get_locint_r();
        loc1=get_locint_r();
	As[loc2] += 
	  As[result]*Trs[loc1];
	As[loc1] += 
	  As[result]*Trs[loc2];
	break;
      case mult_d_a:
        result = get_locint_r();
        loc1=get_locint_r();
        stored_val = get_val_r();
	As[loc1] +=
	  stored_val*As[result];
	break;
      case div_a_a:
        result =get_locint_r();
        loc2=get_locint_r();
        loc1=get_locint_r();
	r0= As[result]/Trs[loc2];
	As[loc1] += r0;
	As[loc2] -=
	  r0*Trs[result];
	break;
      case div_d_a:
        result = get_locint_r();
        loc1=get_locint_r();
        stored_val = get_val_r();
	As[loc1] -= As[result]*
	  Trs[result]/Trs[loc1];
	break;
      case pow_op:
        result = get_locint_r();
        loc1=get_locint_r();
        stored_val = get_val_r();
        if (Trs[loc1]==0.0)
          r0=0.0;
        else
          r0=1.0/Trs[loc1];
	As[loc1] += As[result]*
           stored_val* Trs[result]*r0;
	break;
      case death_not:
        loc2 = get_locint_r();
        loc1 = get_locint_r();
	for (i=loc1;i<=loc2;i++)
	  {
	    As[i] = 0.0;
	    get_taylor(i);
	  }
	break;
      case exp_op:
        result = get_locint_r();
        loc1=get_locint_r();
	As[loc1] +=
	  As[result]*Trs[result];
	break;
      case sin_op:
        result = get_locint_r();
        loc2 = get_locint_r();
        loc1 = get_locint_r();
        stored_val = get_val_r();
	As[loc1] +=
	  As[result]*Trs[loc2];
	break;
      case cos_op:
        result = get_locint_r();
        loc2 = get_locint_r();
        loc1 = get_locint_r();
        stored_val = get_val_r();
	As[loc1] -=
	  As[result]*Trs[loc2];
	break;
      case sqrt_op:
        result = get_locint_r();
        loc1=get_locint_r();
        if (Trs[result]==0.0)
          r0=0.0;
        else
          r0=0.5/Trs[result];
	As[loc1] += r0*As[result];
	break;
      case abs_val:
        result = get_locint_r();
        loc1=get_locint_r();
        if(Trs[loc1]<0.0)
          As[loc1]-=As[result];
        else if(Trs[loc1]>0.0)
          As[loc1]+=As[result];
        break;

      case asin_op:
      case acos_op:
      case atan_op:
      case gen_quad:
        result = get_locint_r();
        loc2 = get_locint_r();
        loc1 = get_locint_r();
        stored_val = get_val_r();
	As[loc1] += 
	  As[result]*Trs[loc2];
	break;      
      case log_op:
        result = get_locint_r();
        loc1=get_locint_r();
	As[loc1] += 
	  As[result]/Trs[loc1];
	break;
 
      /* Vector operations */

      case int_av_av:
        result_v= get_locint_r();
        size = get_locint_r();
        loc1_v = get_locint_r();
	for (l=0;l<size;l++)
	  {
	    loc1 = result_v + l;   /* Location of left-hand-side  */
	    loc2 = loc1_v + l;   /* Location of right-hand-side */
	    /* code for int_adb_a */
	    As[loc2] += As[loc1];
	    As[loc1] = 0.0;
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
	    As[loc1] = 0;
	    get_taylor(loc1);
	  }
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
	    get_taylor(loc1);
	    As[loc2] += As[loc1];
	    As[loc1] = 0.0;
	  }
	break;
      case assign_indvec:
        loc1_v = get_locint_r();
        size = get_locint_r();
	for (l=size-1;l>=0;l--)
	  {
	    loc1 = loc1_v + l;             /* Location of left-hand-side */
	    /* code for assign_ind */
	    get_taylor(loc1);
	    results[indexi] =As[loc1];
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
	    As[loc1] =lagrange[indexd];
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
	    get_taylor(result);
	    As[loc1] -= As[result];
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
	    get_taylor(result);
	    As[loc1] +=As[result];
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
	    get_taylor(result);
	    r0 = As[result];
	    As[result] = 0;
	    As[loc1] += r0*Trs[result];
	    As[result] += r0*Trs[loc1] ; 
	    
	  }
	break;
      case eq_mult_av_d:
        result_v = get_locint_r();
        size = get_locint_r();
        stored_val = get_val_r();
	for (l=size-1;l>=0;l--)
	  {
	    result = result_v + l;      /* Location of the left-hand-side  */
	    /* code for eq_mult_d*/
	    get_taylor(result);
	    As[result] *=stored_val;
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
	    As[loc1]+=As[result];
	    As[loc2]+=As[result];
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
	    As[loc1]+=As[result];
	    As[loc2]-=As[result];
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
	    As[loc2] += As[result]*Trs[loc1];
	    As[loc1] += As[result]*Trs[loc2];
	  }
	As[result] = 0;

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
	    As[loc1] += stored_val*As[result];
	    /* Be Careful here also see line 753 */
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
	    r0= As[result]/Trs[loc2];
	    As[loc1] += r0;
	    As[loc2] -= r0*Trs[result];
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
	    As[loc2] += As[result]*Trs[loc1];
	    As[loc1] += As[result]*Trs[loc2];
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
	    loc2= loc1_v+l;      /* Location of right hand side vectore[l]  */
	    result = result_v + l; /* Location of the result */
	    
	    /* code for mult_a_a */
	    As[loc2] += As[result]*Trs[loc1];
	    As[loc1] += As[result]*Trs[loc2];
	  }
	break;
	
#ifdef conditional
      case cond_assign:
        result = get_locint_r();
        loc2 = get_locint_r();
        loc1 = get_locint_r();
        loc1_v = get_locint_r();
	get_taylor(result);
        if(Trs[loc1_v]>0)
        {
	  As[loc1] += As[result];
	  As[result] = 0.0;
        } /* endif */
        else
        {
          As[loc2] += As[result];
          As[result] = 0.0;
        } /* endelse */
	break;
      case cond_assign_s:
        result = get_locint_r();
        loc1 = get_locint_r();
        loc1_v = get_locint_r();
	get_taylor(result);
        if(Trs[loc1_v]>0)
        {
	  As[loc1] += As[result];
	  As[result] = 0.0;
        } /* endif */
	break;
      case subscript:
        result = get_locint_r();
        loc1 = get_locint_r();
        loc1_v = get_locint_r(); 
	get_taylor(result);
	As[loc1_v+(int)(Trs[loc1])] += As[result];
	As[result] = 0.0;
	break;
      case subscript_l:
        result = get_locint_r();
        loc1 = get_locint_r();
        loc1_v = get_locint_r(); 
	get_taylor(loc1_v+(int)(Trs[loc1]));
	As[result] += As[loc1_v+(int)(Trs[loc1])];
	As[loc1_v+(int)(Trs[loc1])] = 0.0;
	break;
      case subscript_ld:
        loc1 = get_locint_r();
        loc1_v = get_locint_r(); 
        stored_val = get_val_r();
	As[loc1_v+(int)(Trs[loc1])] = 0.0;
	break;
      case m_subscript:
        result = get_locint_r();
        size = get_locint_r();
        loc1 = get_locint_r();
        loc1_v = get_locint_r(); 
        for (l=size-1;l>=0;l--)
        {
          res=result+l;
          arg=loc1_v+(int)(Trs[loc1])*size+l;
	  As[arg] += As[res];
	  As[res] = 0.0;
        } /* endfor */
	break;
      case m_subscript_l:
        result = get_locint_r();
        size = get_locint_r();
        loc1 = get_locint_r();
        loc1_v = get_locint_r(); 
        for (l=size-1;l>=0;l--)
        {
          arg=loc1_v+(int)(Trs[loc1])*size+l;
          res=result+l;
	  get_taylor(arg);
	  As[res] += As[arg];
	  As[arg] = 0.0;
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
          arg=loc1_v+(int)(Trs[loc1])*size+l+loc2;
	  As[arg] = 0.0;
        } /* endfor */
	break;
#endif	


      default:
	/* Die here, we screwed up */
	printf("ADOL-C error: Fatal error in fos_reverse on operation %d\n",
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

