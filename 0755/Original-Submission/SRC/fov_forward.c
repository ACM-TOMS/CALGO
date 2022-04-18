/*
   ----------------------------------------------------------------
   File fov_forward.c of ADOL-C version 1.6 as of January 1,   1995
   ----------------------------------------------------------------
   Contains the routine fov_forward (first-order-vector forward
   mode).

*/

#ifdef __cplusplus
extern "C" {
#endif


#include "dvlparms.h" /* Developers Parameters */
#include "usrparms.h"
#include "oplate.h"
#include "taputil1.h"
#include "taputil2.h"
#include "taputil3.h"
#include "tayutil.h"
#include "adutilsc.h"

static short tag;

static int for_location_cnt;
static int dep_cnt;
static int ind_cnt;

/****************************************************************************/
/* First Order Vector version of the forward mode.                         */
/****************************************************************************/
void fov_forward(short tnum,         /* tape id */
		 int depcheck,       /* consistency chk on # of dependents */
		 int indcheck,       /* consistency chk on # of independents */
		 int p,              /* # of taylor series */
                 double *basepoint,  /* independent variable values */
		 double **argument,  /* Taylor coefficients (input) */
                 double *valuepoint, /* Taylor coefficients (output) */
		 double **taylors)   /* matrix of coifficient vectors */
/* the order of the indices in argument and taylors is [var][taylor] */
{
  unsigned char operation;
  int tape_stats[11];  /* tape stats */

  locint size = 0;
  locint result=0;
  locint result_v=0;
  locint loc1=0;
  locint loc1_v=0;
  locint loc2=0;
  locint loc2_v=0;
  double stored_val=0.0;
  double *d = 0;

  int  l, pl;
  double r0;
  locint arg,arg1,arg2,res;
  int indexi =0; 
  int indexd =0;
  double  *Targ,*Tres,*Tqo,*Targ1,*Targ2,*Tloc1,*Tloc2; 
  double coval,divs;
  int buffer; 
  static int fax, kax;
  static double *Tdum; 
  double* T0;
  double** T;

#ifdef inf_num
double i_num=inf_num;
double i_den=inf_den;
double n_num=non_num;
double n_den=non_den;
double InfVal;
double NoNum;
#endif

  tag = tnum;         /*tag is global which specifies which tape to look at */
  
  tapestats(tag,tape_stats);

  ind_cnt = tape_stats[0];
  dep_cnt = tape_stats[1];
  for_location_cnt = tape_stats[2];
  buffer =  tape_stats[4];

  set_buf_size(buffer);
  T0 = (double*)malloc(for_location_cnt*sizeof(double));
  T = (double**)myalloc2(for_location_cnt,p); 
  if ((depcheck != dep_cnt)||(indcheck != ind_cnt))
    {
      printf("ADOL-C error: forward sweep on tape %d  aborted!\n",tag);
      printf("Number of dependent and/or independent variables passed to forward is\ninconsistant with number recorded on tape %d \n",tag);
      exit (-1);
    } /* endif */
  

  /* Initialize the Forward Sweep */
  init_for_sweep(tag);


  operation=get_op_f();
  while (operation !=end_of_tape)
    {
      switch (operation){
      case end_of_op:
        get_op_block_f();
        operation=get_op_f(); 
        /* Skip next operation, it's another end_of_op */
        break;
      case end_of_int:
        get_loc_block_f();
        break;
      case end_of_val:
        get_val_block_f();
        break;
      case start_of_tape:
      case end_of_tape:
	break;
      case int_adb_a:  /* initialize an adouble */	
        arg = get_locint_f();
        res = get_locint_f();
	Targ = T[arg];
        Tres = T[res]; 
        T0[res]= T0[arg];
        for (l=0;l<p;l++)
	  Tres[l]=Targ[l];
	break;
      case assign_a:  /* assign an adouble variable an adouble value. (=) */
        arg = get_locint_f();
        res = get_locint_f();
	Targ = T[arg];
        Tres = T[res];
        T0[res]= T0[arg];
        for (l=0;l<p;l++)
          Tres[l]=Targ[l];
        break;
      case int_adb_d: /* initialize an double */
        res = get_locint_f();
	Tres = T[res];
	T0[res]=get_val_f();
        for (l=0;l<p;l++)
	  Tres[l]=0;
	break;
      case assign_d: /* assign an adouble variable a float value. (=) */
        res = get_locint_f();
	Tres = T[res];
	T0[res]=get_val_f();
        for (l=0;l<p;l++)
          Tres[l]=0;
	break;
      case assign_ind: /* assign an adouble variable an independent 
                          float value (<<=) */
        res = get_locint_f();
	Tres = T[res];
        T0[res]=basepoint[indexi];
        for (l=0;l<p;l++)
	  Tres[l]=argument[indexi][l];
	++indexi;
	break;
      case assign_dep: /* assign a float variable a 
                          dependent adouble value. (>>=) */
	res = get_locint_f();
	Tres = T[res];
        valuepoint[indexd]=T0[res];
        if (taylors != 0 )  
        {
          for (l=0;l<p;l++)
	    taylors[indexd][l] = Tres[l];
        } /* endif */
	indexd++;
	break;
      case eq_plus_d: /* Add a floating point to an adouble. (+=) */
        T0[get_locint_f()]+= get_val_f();
	break;
      case eq_plus_a: /* Add an adouble to another adouble. (+=) */
        arg=get_locint_f();
        res=get_locint_f();
	Tres = T[res];
	Targ = T[arg];
        T0[res]+= T0[arg];
        for (l=0;l<p;l++)
	  Tres[l]+=Targ[l];
	break;
      case eq_min_d: /* Subtract a floating point from an adouble. (-=) */
	T0[get_locint_f()]-= get_val_f();
	break;
      case eq_min_a: /* Subtract an adouble from another adouble. (-=) */
        arg=get_locint_f();
        res=get_locint_f();
	Tres = T[res];
	Targ = T[arg];   
        T0[res]-=T0[arg];
        for (l=0;l<p;l++)
	  Tres[l]-=Targ[l];
	break;
      case eq_mult_d:  /* Multiply an adouble by a float. (*=) */
        res=get_locint_f();
        coval = get_val_f();
	Tres = T[res];
        T0[res]*=coval;
        for (l=0;l<p;l++)
	  Tres[l]*=coval;
	break;
      case eq_mult_a: /* Multiply one adouble by another. (*=) */
        arg=get_locint_f();
        res=get_locint_f();
	Tres = T[res];
	Targ = T[arg];
        for( l = 0; l < p; l++)
          Tres[l]=T0[res]*Targ[l]+Tres[l]*T0[arg];
        T0[res]*=T0[arg];
	break;
      case plus_a_a: /* Add two adoubles. (+) */
        arg1=get_locint_f();
        arg2=get_locint_f();
        res =get_locint_f();
	Tres = T[res];
	Targ1 = T[arg1];
	Targ2 = T[arg2];
        T0[res]=T0[arg1]+T0[arg2];
        for( l = 0; l < p; l++)
	  Tres[l]=Targ1[l]+Targ2[l];
	break;
      case plus_d_a: /* Add an adouble and a double. (+) */
        arg=get_locint_f();
        res = get_locint_f();
	Tres = T[res];
	Targ = T[arg];
	T0[res]=T0[arg]+get_val_f();
        for( l = 0; l < p; l++)
	  Tres[l]=Targ[l];
	break;
      case min_a_a: /* Subtraction of two adoubles. (-) */
        arg1=get_locint_f();
        arg2=get_locint_f();
        res =get_locint_f();
	Tres = T[res];
	Targ1 = T[arg1];
	Targ2 = T[arg2];
        T0[res]=T0[arg1]-T0[arg2];
        for( l = 0; l < p; l++)
          Tres[l]=Targ1[l]-Targ2[l];
	break;
      case min_d_a: /* Subtract an adouble from a double. (-) */
        arg=get_locint_f();
        res = get_locint_f();
        Tres = T[res];
        Targ = T[arg];
	T0[res]=get_val_f()-T0[arg];
        for( l = 0; l < p; l++)
	  Tres[l] = -Targ[l];
	break;
      case mult_a_a: /* Multiply two adoubles. (*) */
        arg1=get_locint_f();
        arg2=get_locint_f();
        res =get_locint_f();
        Tres = T[res];
        Targ1 = T[arg1];
        Targ2 = T[arg2];
        T0[res]=T0[arg1]*T0[arg2];
        for( l = 0; l < p; l++)
          Tres[l]=T0[arg1]*Targ2[l]+Targ1[l]*T0[arg2];
	break;
      case mult_d_a: /* Multiply an adouble by a double. (*) */
        arg=get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        Tres = T[res];
        Targ = T[arg];
        T0[res]=T0[arg]*coval;
        for( l = 0; l < p; l++)
	  Tres[l]=Targ[l]*coval;
	break;
      case div_a_a: /* Divide an adouble by an adouble. (/) */
        arg1=get_locint_f();
        arg2=get_locint_f();
        res =get_locint_f();
        Tres = T[res];
        Targ1 = T[arg1];
        Targ2 = T[arg2];
	divs   = 1/T0[arg2];
        T0[res]=T0[arg1]*divs;
        for( l = 0; l < p; l++)
	  Tres[l]=Targ1[l]*divs+T0[res]*(-Targ2[l]*divs);
	break;
      case div_d_a: /* Division double - adouble. (/) */
        arg=get_locint_f();
        res = get_locint_f();
        Tres = T[res];
        Targ = T[arg];
	divs = 1/T0[arg];
	T0[res]=get_val_f()*divs;
        for( l = 0; l < p; l++)
	  Tres[l]=T0[res]*(-Targ[l]*divs);
	break;
      case sqrt_op: /* Compute sqrt of adouble. */
        arg=get_locint_f();
        res= get_locint_f();
        Targ = T[arg];
	Tres = T[res];
	T0[res]=sqrt(T0[arg]);
        for( l = 0; l < p; l++){
          if (T0[arg]==0.0){
#ifdef inf_num
  InfVal=i_num/i_den;
  NoNum=n_num/n_den;
#endif
            if (Targ[l]>0.0)
              r0=InfVal;
            else if (Targ[l]<0.0) 
              r0=NoNum;
            else
              r0=0.0;
          } /* end if */
          else {
            r0 = 0.5/T0[res];
          } /* end else */
	  Tres[l]=r0*Targ[l];
        } /* end for */
	break;
      case abs_val: /* Compute fabs of adouble. */
        arg=get_locint_f();
        res= get_locint_f();
        Targ = T[arg];
        Tres = T[res];
        T0[res]=fabs(T0[arg]);
        for( l = 0; l < p; l++)
          if (T0[arg]<0.0)
            Tres[l]= -Targ[l];
          else if(T0[arg]>0.0)
            Tres[l]=Targ[l];
          else
            Tres[l]= fabs(Targ[l]);
        break;

      case exp_op: /* exponent operation */
        arg=get_locint_f();
        res= get_locint_f();
	Tres = T[res];
	Targ = T[arg];
	T0[res]=exp(T0[arg]);
        for( l = 0; l < p; l++)
	  Tres[l]=T0[res]*Targ[l];
	break;
      case sin_op: /* sine operation */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
	Targ = T[arg2];
	Tres = T[res]; 
	Tqo = T[arg1];
	T0[arg2]=cos(T0[arg1]);
	T0[res] = sin(T0[arg1]);
        for( l = 0; l < p; l++)
	{
	    Targ[l] = -(T0[res]*Tqo[l]);
	    Tres[l]=T0[arg2]*Tqo[l];
        } /* endfor */
	break;
      case cos_op:  /* cosine operation */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
	Targ = T[arg2];
	Tres = T[res];
	Tqo = T[arg1];
        T0[arg2]=sin(T0[arg1]);
        T0[res] = cos(T0[arg1]);
        for( l = 0; l < p; l++)
        {
	  Targ[l] = T0[res]*Tqo[l];
	  Tres[l] = -(T0[arg2]*Tqo[l]);
        } /* endfor */
	break;
      case asin_op: 
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        Targ = T[arg2];
        Tres = T[res];
        Tqo = T[arg1];
	T0[res]=asin(T0[arg1]);
        for( l = 0; l < p; l++)
	  Tres[l]=T0[arg2]*Tqo[l];
	break;
      case acos_op:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        Targ = T[arg2];
        Tres = T[res];
        Tqo = T[arg1];
	T0[res]=acos(T0[arg1]);
        for( l = 0; l < p; l++)
	  Tres[l]=T0[arg2]*Tqo[l];
	break;
      case atan_op:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        Targ = T[arg2];
        Tres = T[res];
        Tqo = T[arg1];
	T0[res]=atan(T0[arg1]);
        for( l = 0; l < p; l++)
	  Tres[l]=T0[arg2]*Tqo[l];
	break;
      case gen_quad:  
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        Targ = T[arg2];
        Tres = T[res];
        Tqo = T[arg1];
	T0[res]=coval;
        for( l = 0; l < p; l++)
	  Tres[l]=T0[arg2]*Tqo[l];
	break;
      case log_op:
        arg=get_locint_f();
        res= get_locint_f();
	Tres = T[res];
	Targ = T[arg];
	divs = 1.0/T0[arg];
	T0[res]=log(T0[arg]);
        for( l = 0; l < p; l++)
	  Tres[l]=Targ[l]*divs;
	break;
      case pow_op:
        arg=get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
	Tres = T[res];
	Targ = T[arg];
	T0[res] = pow(T0[arg],coval);
        for( l = 0; l < p; l++){
          if (T0[arg]==0.0){
#ifdef inf_num
  InfVal=i_num/i_den;
  NoNum=n_num/n_den;
#endif
            if (Targ[l]>0.0)
              r0=InfVal;
            else if (Targ[l]<0.0)
              r0=NoNum;
            else
              r0=0.0;
          } /* end if */
          else {
            r0 = 1.0/T0[arg];
          } /* end else */
          Tres[l] = T0[res]*Targ[l]*coval*r0;
        } /* end for */
	break;

      
      case int_av_av:
        loc1_v = get_locint_f();
        size = get_locint_f();
        result_v= get_locint_f();
        for (l=0;l<size;l++)
        {
          res = result_v + l;   /* Location of left-hand-side  */
          arg = loc1_v + l;   /* Location of right-hand-side */
          Targ = T[arg];
          Tres = T[res];
          T0[res]= T0[arg];
          for (pl=0;pl<p;pl++)
            Tres[pl]=Targ[pl];
        } /* endfor */
        break;
      case assign_dv:
        size = get_locint_f();
        loc1_v = get_locint_f();
        d = get_val_v_f(size);
        for (l=0;l<size;l++)
        {
            res = loc1_v + l;           /* Location of left-hand-side */
            Tres = T[res];
            T0[res]=d[l];   
            for (pl=0;pl<p;pl++)
	      Tres[pl]=0;
        } /* endfor */
        break;
      case assign_av:
        loc2_v = get_locint_f();
        size = get_locint_f();
        loc1_v = get_locint_f();
        for (l=0;l<size;l++)
        {
          res = loc1_v + l;             /* Location of the left-hand-side  */
          arg = loc2_v + l;             /* Location of the right-hand-side */
          Targ = T[arg];
          Tres = T[res];
          T0[res]=T0[arg];
          for (pl=0;pl<p;pl++)
            Tres[pl]=Targ[pl];
        } /* endfor */
        break;
      case assign_indvec:
        size = get_locint_f();
        loc1_v = get_locint_f();
        for (l=0;l<size;l++)
        {
          res = loc1_v + l;             /* Location of the left-hand-side */
          Tres = T[res];
          T0[res]=basepoint[indexi];
          for (pl=0;pl<p;pl++)
            Tres[pl]=argument[indexi][pl];
          ++indexi;
        } /* endfor */
        break;
      case assign_depvec:
        size = get_locint_f();
        loc1_v = get_locint_f();
        for (l=0;l<size;l++)
        {
          res = loc1_v + l;             /* Location of the left-hand-side */
          Tres = T[res];
          valuepoint[indexd]=T0[res];
          if (taylors != 0 )  
          {
            for (pl=0;pl<p;pl++)
	      taylors[indexd][pl] = Tres[pl];
          } /* endif */
          indexd++;
        } /* endfor */
        break;
      case eq_min_av:
        loc1_v = get_locint_f();
        size = get_locint_f();
        result_v = get_locint_f();
        for (l=0;l<size;l++)
        {
          res = result_v + l;            /* Location of the left-hand-side  */
          arg = loc1_v + l;            /* Location on the right-hand-side */
          Tres = T[res];
          Targ = T[arg];
          T0[res]-=T0[arg];
          for (pl=0;pl<p;pl++)
            Tres[pl] -= Targ[pl];
        } /* endfor */
        break;
      case eq_plus_av:
        loc1_v = get_locint_f();
        size = get_locint_f();
        result_v = get_locint_f();
        for (l=0;l<size;l++)
        {
          res = result_v + l;            /* Location of the left-hand-side  */
          arg = loc1_v   + l;            /* Location on the right-hand-side */
          Tres = T[res];
          Targ = T[arg];
          T0[res]+= T0[arg];
          for (pl=0;pl<p;pl++)
            Tres[pl]+=Targ[pl];
        } /* endfor */
        break;
      case eq_mult_av_a:
        arg = get_locint_f();
        size = get_locint_f();
        result_v = get_locint_f();
        for (l=0;l<size;l++)
        {
          res = result_v + l;      /* Location of the left-hand-side  */
          Tres = T[res];
          Targ = T[arg];
          for (pl=0;pl<p;pl++)
	    Tres[pl]=T0[res]*Targ[pl]+Tres[pl]*T0[arg];
          T0[res]*=T0[arg]; 
        } /* endfor */
        break;
      case eq_mult_av_d:
        size = get_locint_f();
        result_v = get_locint_f();
        coval = get_val_f();
        for (l=0;l<size;l++)
        {
          res = result_v + l;      /* Location of the left-hand-side  */
          Tres = T[res];
          T0[res]*=coval;
          for (pl=0;pl<p;pl++)
            Tres[pl]*=coval;
        } /* endfor */
        break;
      case plus_av_av:
        loc1_v   = get_locint_f();
        loc2_v   = get_locint_f();
        size     = get_locint_f();
        result_v = get_locint_f();
        for (l=0;l<size;l++)
        {
          arg2   = loc2_v    + l;       /* Location of var 2  */
          arg1   = loc1_v    + l;       /* Location of var 1  */
          res = result_v  + l;       /* Location of result */
          Tres = T[res];
          Targ1 = T[arg1];
          Targ2 = T[arg2];
          T0[res]=T0[arg1]+T0[arg2];
          for (pl=0;pl<p;pl++)
            Tres[pl]=Targ1[pl]+Targ2[pl];
        } /* endfor */
        break;
      case sub_av_av:
        loc1_v   = get_locint_f();
        loc2_v   = get_locint_f();
        size     = get_locint_f();
        result_v = get_locint_f();
        for (l=0;l<size;l++)
        {
          arg2   = loc2_v    + l;       /* Location of var 2  */
          arg1   = loc1_v    + l;       /* Location of var 1  */
          res = result_v  + l;       /* Location of result */
          Tres = T[res];
          Targ1 = T[arg1];
          Targ2 = T[arg2];
          T0[res]=T0[arg1]-T0[arg2];
          for (pl=0;pl<p;pl++)
	    Tres[pl]=Targ1[pl]-Targ2[pl];
        } /* endfor */
        break;
      case dot_av_av:
        loc1_v   = get_locint_f();
        loc2_v   = get_locint_f();
        size     = get_locint_f();
        res = get_locint_f();
        Tres = T[res];
        T0[res]=0;
        for (pl=0;pl<p;pl++)
          Tres[pl]=0.0;
        for (l=0;l<size;l++)
        {
          arg2 = loc2_v + l;
          arg1 = loc1_v + l;
          Targ1 = T[arg1];
          Targ2 = T[arg2];
          Tres = T[res];
          T0[res]+=T0[arg1]*T0[arg2];
          for (pl=0;pl<p;pl++)
	    Tres[pl]+=T0[arg1]*Targ2[pl]+Targ1[pl]*T0[arg2];
        } /* endfor */
        break;
      case mult_d_av:
        loc1_v = get_locint_f();
        size = get_locint_f();
        result_v = get_locint_f();
        coval =get_val_f();
        for (l=0;l<size;l++)
        {
          arg = loc1_v   + l;   /* Location on the right-hand-side */
          res = result_v + l;   /* location of the result */
          Tres = T[res];
          Targ = T[arg];
          T0[res]=T0[arg]*coval;
          for (pl=0;pl<p;pl++)
            Tres[pl]=Targ[pl]*coval;
        } /* endfor */
        break;
      case div_av_a:
        loc1_v   = get_locint_f();
        arg2     = get_locint_f();
        size     = get_locint_f();
        result_v = get_locint_f();
        for (l=0;l<size;l++)
        {
          arg1 = loc1_v   + l;  /* Location of the right-hand-side vector */
          res = result_v + l;  /* Location of the result */
          /* code for div_a_a */
          Tres = T[res];
          Targ1 = T[arg1];
          Targ2 = T[arg2];
          divs   = 1/T0[arg2];
          T0[res]=T0[arg1]*divs;
          for (pl=0;pl<p;pl++)
	    Tres[pl]=(Targ1[pl]-T0[res]*Targ2[pl])*divs;
        } /* endfor */
        break;
      case mult_av_a:
        loc1_v   = get_locint_f();
        arg2     = get_locint_f();
        size     = get_locint_f();
        result_v = get_locint_f();
        for (l=0;l<size;l++)
         {
          arg1 = loc1_v   + l; /* Location of the right-hand-side vector */
          res = result_v + l; /* Location of the result */
          Tres = T[res];
          Targ1 = T[arg1];
          Targ2 = T[arg2];
          T0[res]=T0[arg1]*T0[arg2];
          for (pl=0;pl<p;pl++)
	    Tres[pl]=T0[arg1]*Targ2[pl]+Targ1[pl]*T0[arg2];
        } /* endfor */
        break;
      case mult_a_av:
        loc1_v   = get_locint_f();
        arg1     = get_locint_f();
        size     = get_locint_f();
        result_v = get_locint_f();
        for (l=0;l<size;l++)
        {
          arg2= loc1_v+l;   /* Location of right hand side vector */
          res = result_v + l; /* Location of the result */
          Tres = T[res];
          Targ1 = T[arg1];
          Targ2 = T[arg2];
          T0[res]=T0[arg1]*T0[arg2];
          for (pl=0;pl<p;pl++)
            Tres[pl]=T0[arg1]*Targ2[pl]+Targ1[pl]*T0[arg2];
        } /* endfor */
        break;

#ifdef conditional
      case cond_assign:
        arg = get_locint_f();
        loc1 = get_locint_f();
        loc2 = get_locint_f();
        res = get_locint_f();
        Tres = T[res];
        Targ = T[arg];
        Tloc1 = T[loc1];
        Tloc2 = T[loc2];
        if (T0[arg]>0)
        {
          T0[res]=T0[loc1];
	  for( l = 0; l < p; l++)
	    Tres[l]=Tloc1[l];
        } /* endif */
        else
        {
          T0[res]=T0[loc2];
	  for( l = 0; l < p; l++)
	    Tres[l]=Tloc2[l];
        } /* endelse */
        break;
      case cond_assign_s:
        arg = get_locint_f();
        loc1 = get_locint_f();
        res = get_locint_f();
        Tres = T[res];
        Targ = T[arg];
        Tloc1 = T[loc1];
        if (T0[arg]>0)
	{
          T0[res]=T0[loc1];
	  for(l = 0; l < p; l++)
	    Tres[l]=Tloc1[l];
        } /* endif */
        break;
      case subscript:
        loc1_v=get_locint_f();  /* Base */
        loc1=get_locint_f();/* pointer to the variable containing the offset */
        res=get_locint_f();
        Tres=T[res];
        arg=loc1_v+(int)(T0[loc1]);
        Targ = T[arg];
        T0[res]=T0[arg];
        for( l = 0; l < p; l++)
          Tres[l]=Targ[l];
        break;
      case subscript_l:
        loc1_v=get_locint_f();  /* Base */
	loc1=get_locint_f();
        arg=loc1_v+(int)(T0[loc1]);
        res=get_locint_f();
        Tres = T[res];
        Targ = T[arg];
        T0[arg]=T0[res];
        for( l = 0; l < p; l++)
          Targ[l]=Tres[l];
        break;
      case subscript_ld:
        loc1_v=get_locint_f();  /* Base */
        loc1=get_locint_f();  /* pointer to the variable containing the offset */
        res = loc1_v+(int)(T0[loc1]);
        Tres = T[res];
        T0[res]=get_val_f();
        for( l = 0; l < p; l++)
          Tres[l]=0;
        break;
      case m_subscript:
        loc1_v=get_locint_f();  /* Base */
        loc1=get_locint_f();/* pointer to the variable containing the offset */
        size=get_locint_f();
        result=get_locint_f();
        for (l=0;l<size;l++)
        {
          res=result+l;
          Tres=T[res];
          arg=loc1_v+(int)(T0[loc1])*size+l;
          Targ = T[arg];
          T0[res]=T0[arg];
          for(pl=0;pl<p;pl++)
            Tres[pl]=Targ[pl];
        } /* endfor */
        break;
      case m_subscript_l:
        loc1_v=get_locint_f();  /* Base */
	loc1=get_locint_f();
        size=get_locint_f();
        arg=get_locint_f(); /* RHS */
        for (l=0;l<size;l++)
        {
          res=loc1_v+(int)(T0[loc1])*size+l;
          Tres = T[res];
          Targ = T[arg+l];
          T0[res]=T0[arg+l];
          for(pl=0;pl<p;pl++)
            Tres[pl]=Targ[pl];
        } /* endfor */
        break;
      case m_subscript_ld:
        loc1_v=get_locint_f();  /* Base */
        loc1=get_locint_f();  /* pointer to the variable containing the offset */
        loc2=get_locint_f(); /* offset in the vector itself */
        size=get_locint_f();
        d = get_val_v_f(size);
        for (l=0;l<size;l++)
        {
          res = loc1_v+(int)(T0[loc1])*size+l+loc2;
          Tres = T[res];
          T0[res]=d[l];
          for(pl=0;pl<p;pl++)
            Tres[pl]=0;
        } /* endfor */
        break;
#endif

      case death_not:
        loc1=get_locint_f();
        loc1=get_locint_f();
  	break;
      default:
	/* Die here, we screwed up */
	printf("ADOL-C error: Fatal error in fov_forward for op %d\n",operation);
	break;
	
      } /* endswitch */
      operation=get_op_f();
    }  /* endwhile */
  end_sweep();
} /* end fov_forward */



#ifdef __cplusplus
}
#endif

