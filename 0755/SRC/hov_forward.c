/*
   ----------------------------------------------------------------
   File hov_forward.c of ADOL-C version 1.6 as of January 1,   1995
   ----------------------------------------------------------------
   Contains the routine hov_forward (higher-order-vector forward
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
static int degree;

/****************************************************************************/
/* Higher Order Vector version of the forward mode.                         */
/****************************************************************************/
void hov_forward(short tnum,         /* tape id */
		 int depcheck,       /* consistency chk on # of dependents */
		 int indcheck,       /* consistency chk on # of independents */
		 int gdegree,        /* highest derivative degree */
		 int p,              /* # of taylor series */
                 double *basepoint,  /* independent variable values */
		 double ***argument, /* Taylor coefficients (input) */
                 double *valuepoint, /* Taylor coefficients (output) */
		 double ***taylors)  /* matrix of coifficient vectors */
/* the order of the indices in argument and taylors is [var][taylor][deriv] */
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

  int k, i, j ,l, rloc,pl;
  double r0;
  locint arg,arg1,arg2,res;
  int indexi =0; 
  int indexd =0;
  double  *Targ,*Tres,*Tqo,*Targ1,*Targ2,*Tloc1,*Tloc2; 
  double coval,divs;
  int buffer; 
  double x,y;
  static int fax, kax;
  static double *Tdum; 
  double *z;
  int even;
  double* T0;
  double*** T;

#ifdef inf_num
double i_num=inf_num;
double i_den=inf_den;
double n_num=non_num;
double n_den=non_den;
double InfVal;
double NoNum;
#endif

  degree = gdegree;
  tag = tnum;         /*tag is global which specifies which tape to look at */
  k = degree + 1;
  
  tapestats(tag,tape_stats);

  ind_cnt = tape_stats[0];
  dep_cnt = tape_stats[1];
  for_location_cnt = tape_stats[2];
  buffer =  tape_stats[4];

  set_buf_size(buffer);
  T0 = (double*)malloc(for_location_cnt*sizeof(double));
  T = (double***)myalloc3(for_location_cnt,p,degree); 
  z = (double*)malloc(k*sizeof(double));
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
	Targ = *T[arg];
        Tres = *T[res]; 
        T0[res]= T0[arg];
        for (l=0;l<p;l++)
        {
          for (i=0;i<k-1;i++)
	    Tres[i]=Targ[i];
          Targ+=degree;
          Tres+=degree;
        } /* endfor*/
	break;
      case assign_a:  /* assign an adouble variable an adouble value. (=) */
        arg = get_locint_f();
        res = get_locint_f();
	Targ = *T[arg];
        Tres = *T[res];
        T0[res]= T0[arg];
        for (l=0;l<p;l++)
        {
          for (i=0;i<k-1;i++)
            Tres[i]=Targ[i];
          Targ+=degree;
          Tres+=degree;
        } /* endfor*/
        break;
      case int_adb_d: /* initialize an double */
        res = get_locint_f();
        coval = get_val_f();
	Tres = *T[res];
	T0[res]=coval;
        for (l=0;l<p;l++)
        {
	  for (i=0;i<k-1;i++)
	    Tres[i]=0;
          Tres+=degree;
        } /* endfor*/
	break;
      case assign_d: /* assign an adouble variable a float value. (=) */
        res = get_locint_f();
        coval = get_val_f();
	Tres = *T[res];
	T0[res]=coval;
        for (l=0;l<p;l++)
        {
          for (i=0;i<k-1;i++)
            Tres[i]=0;
          Tres+=degree;
        } /* endfor*/
	break;
      case assign_ind: /* assign an adouble variable an independent 
                          float value (<<=) */
        res = get_locint_f();
	Tres = *T[res];
        T0[res]=basepoint[indexi];
        for (l=0;l<p;l++)
        {
	  for (i=0;i<k-1;i++)
	    Tres[i]=argument[indexi][l][i];
          Tres+=degree;
        } /* endfor*/
	++indexi;
	break;
      case assign_dep: /* assign a float variable a 
                          dependent adouble value. (>>=) */
	res = get_locint_f();
	Tres = *T[res];
        valuepoint[indexd]=T0[res];
        if (taylors != 0 )  
        {
          for (l=0;l<p;l++)
          {
	    for (i=0;i<k-1;i++)
	      taylors[indexd][l][i] = Tres[i];
            Tres+=degree;
          } /* endfor*/
        } /* endif */
	indexd++;
	break;
      case eq_plus_d: /* Add a floating point to an adouble. (+=) */
        res=get_locint_f();
        coval = get_val_f();
        T0[res]+= coval;
	break;
      case eq_plus_a: /* Add an adouble to another adouble. (+=) */
        arg=get_locint_f();
        res=get_locint_f();
	Tres = *T[res];
	Targ = *T[arg];
        T0[res]+= T0[arg];
        for (l=0;l<p;l++)
        {
          for (i=0;i<k-1;i++)
	    Tres[i]+=Targ[i];
          Tres+= degree;
          Targ+=degree;
        } /* endfor */
	break;
      case eq_min_d: /* Subtract a floating point from an adouble. (-=) */
        res=get_locint_f();
        coval = get_val_f();
	T0[res]-=coval;
	break;
      case eq_min_a: /* Subtract an adouble from another adouble. (-=) */
        arg=get_locint_f();
        res=get_locint_f();
	Tres = *T[res];
	Targ = *T[arg];   
        T0[res]-=T0[arg];
        for (l=0;l<p;l++)
        {
	  for (i=0;i<k-1;i++)
	    Tres[i]-=Targ[i];
          Tres+= degree;
          Targ+=degree;
        } /* endfor */
	break;
      case eq_mult_d:  /* Multiply an adouble by a float. (*=) */
        res=get_locint_f();
        coval = get_val_f();
	Tres = *T[res];
        T0[res]*=coval;
        for (l=0;l<p;l++)
        {
	  for (i=0;i<k-1;i++)
	    Tres[i]*=coval;
          Tres+= degree;
        } /* endfor */
	break;
      case eq_mult_a: /* Multiply one adouble by another. (*=) */
        arg=get_locint_f();
        res=get_locint_f();
	Tres = *T[res];
	Targ = *T[arg];
        for( l = 0; l < p; l++)
        {
          for (i=k-2;i>=0;i--)
          {
  	    x=T0[res]*Targ[i];
   	    for (j=0;j<i;j++)
              x+=Tres[j]*Targ[i-j-1];
            x+=Tres[i]*T0[arg]; 
            Tres[i]=x;
          } /* endfor */
          Tres += degree;
          Targ += degree;
        } /* endfor */
        T0[res]*=T0[arg];
	break;
      case plus_a_a: /* Add two adoubles. (+) */
        arg1=get_locint_f();
        arg2=get_locint_f();
        res =get_locint_f();
	Tres = *T[res];
	Targ1 = *T[arg1];
	Targ2 = *T[arg2];
        T0[res]=T0[arg1]+T0[arg2];
        for( l = 0; l < p; l++)
        {
	  for (i=0;i<k-1;i++)
	    Tres[i]=Targ1[i]+Targ2[i];
          Tres += degree;
          Targ1+= degree;
          Targ2+= degree;
        } /* endfor */
	break;
      case plus_d_a: /* Add an adouble and a double. (+) */
        arg=get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
	Tres = *T[res];
	Targ = *T[arg];
	T0[res]=T0[arg]+coval;
        for( l = 0; l < p; l++)
        {
	  for (i=0;i<k-1;i++)
	    Tres[i]=Targ[i];
          Tres += degree;
          Targ += degree;
        } /* endfor */
	break;
      case min_a_a: /* Subtraction of two adoubles. (-) */
        arg1=get_locint_f();
        arg2=get_locint_f();
        res =get_locint_f();
	Tres = *T[res];
	Targ1 = *T[arg1];
	Targ2 = *T[arg2];
        T0[res]=T0[arg1]-T0[arg2];
        for( l = 0; l < p; l++)
        {
          for (i=0;i<k-1;i++)
            Tres[i]=Targ1[i]-Targ2[i];
          Tres += degree;
          Targ1+= degree;
          Targ2+= degree;
        } /* endfor */
	break;
      case min_d_a: /* Subtract an adouble from a double. (-) */
        arg=get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        Tres = *T[res];
        Targ = *T[arg];
	T0[res]=coval-T0[arg];
        for( l = 0; l < p; l++)
        {
	  for (i=0;i<k-1;i++)
	    Tres[i] = -Targ[i];
          Tres += degree;
          Targ += degree;
        } /* endfor */
	break;
      case mult_a_a: /* Multiply two adoubles. (*) */
        arg1=get_locint_f();
        arg2=get_locint_f();
        res =get_locint_f();
        Tres = *T[res];
        Targ1 = *T[arg1];
        Targ2 = *T[arg2];
        T0[res]=T0[arg1]*T0[arg2];
        for( l = 0; l < p; l++)
        {
	  for (i=0;i<k-1;i++)
	  {
            x=T0[arg1]*Targ2[i];
            for (j=0;j<i;j++)
              x+=Targ1[j]*Targ2[i-j-1];
            x+=Targ1[i]*T0[arg2];
            Tres[i]=x;
	  } /* endfor*/
          Tres += degree;
          Targ1 += degree;
 	  Targ2 += degree;
        } /* endfor */
	break;
      case mult_d_a: /* Multiply an adouble by a double. (*) */
        arg=get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        Tres = *T[res];
        Targ = *T[arg];
        T0[res]=T0[arg]*coval;
        for( l = 0; l < p; l++)
        {
	  for (i=0;i<k-1;i++)
	    Tres[i]=Targ[i]*coval;
          Tres += degree;
          Targ += degree;
        } /* endfor */
	break;
      case div_a_a: /* Divide an adouble by an adouble. (/) */
        arg1=get_locint_f();
        arg2=get_locint_f();
        res =get_locint_f();
        Tres = *T[res];
        Targ1 = *T[arg1];
        Targ2 = *T[arg2];
	divs   = 1/T0[arg2];
        T0[res]=T0[arg1]*divs;
        for( l = 0; l < p; l++)
        {
	  for (i=0;i<k-1;i++)
	  {
	    x=Targ1[i]*divs;
            x+=T0[res]*(-Targ2[i]*divs);
	    z[i]= -Targ2[i]*divs;
	    for (j=0;j<i;j++)
	      x+=Tres[j]*z[i-j-1];
	    Tres[i]=x;
	  } /* endfor */
          Tres += degree;
          Targ1 += degree;
          Targ2 += degree;
        } /* endfor */
	break;
      case div_d_a: /* Division double - adouble. (/) */
        arg=get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        Tres = *T[res];
        Targ = *T[arg];
	divs = 1/T0[arg];
	T0[res]=coval*divs;
        for( l = 0; l < p; l++)
        {
	  for (i=0;i<k-1;i++)
	  {
	    x=0;
	    z[i]=-Targ[i]*divs;
            x+=T0[res]*z[i];
	    for (j=0;j<i;j++)
	      x+=Tres[j]*z[i-j-1];
	    Tres[i] = x;
	  } /* endfor */
          Tres += degree;
          Targ += degree;
        } /* endfor */
	break;
      case sqrt_op: /* Compute sqrt of adouble. */
        arg=get_locint_f();
        res= get_locint_f();
        Targ = *T[arg];
	Tres = *T[res];
	T0[res]=sqrt(T0[arg]);
        for( l = 0; l < p; l++)
        {
          if (Targ[0]==0.0){
#ifdef inf_num
  InfVal=i_num/i_den;
  NoNum=n_num/n_den;
#endif
            r0=0.0;
            for (i=1;i<k;i++){
              if (Targ[i]>0.0){
                r0=InfVal;
                i=k;
              } /* end if */
              if (Targ[i]<0.0){
                r0=NoNum;
                i=k;
              } /* end if */
            } /* end for */
          } /* end if */
          else {
            r0 = 0.5/T0[res];
          } /* end else */
	  even =0;
	  for (i=1;i<k;i++)
	  {
	    x=0;
	    for (j=1;2*j<i;j++)
	      x+=Tres[j-1]*Tres[i-j-1];
	    x *= 2;
	    if(even) 
              x += Tres[i/2-1]*Tres[i/2-1];
	    even = !even;
	    Tres[i-1]=r0*(Targ[i-1]-x);
	  } /* endfor */
          Tres += degree;
          Targ += degree;
        } /* endfor */
	break;
      case abs_val: /* Compute fabs of adouble. */
        arg=get_locint_f();
        res= get_locint_f();
        Targ = *T[arg];
        Tres = *T[res];
        T0[res]=fabs(T0[arg]);
        y=0.0;
        if (Targ[0]!=0.0)
        {
          if (Targ[0]<0.0)
            y = -1.0;
          else
            y=1.0;
        }
        for( l = 0; l < p; l++)
        {
          x=y;
          for (i=0;i<k-1;i++)
          {
            if((x==0.0) && (Targ[i]!=0.0))
            {
              if (Targ[i]<0.0)
                x=-1.0;
              else
                x=1.0;
            }
            Tres[i] = x*Targ[i];
          } /* endfor */
          Tres += degree;
          Targ += degree;
        } /* endfor */
        break;

      case exp_op: /* exponent operation */
        arg=get_locint_f();
        res= get_locint_f();
	Tres = *T[res];
	Targ = *T[arg];
	T0[res]=exp(T0[arg]);
        for( l = 0; l < p; l++)
        {
	  for (i=0;i<k-1;i++)
	  { 
            double x=0;
	    z[i] = (i+1)*Targ[i];
            x+= T0[res]*z[i];
	    for (j=0;j<i;j++)
	      x+= Tres[j]*z[i-j-1];
	    Tres[i]=x/(i+1);
	  } /* endfor */
          Tres += degree;
          Targ += degree;
        } /* endfor */
	break;
      case sin_op: /* sine operation */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
	Targ = *T[arg2];
	Tres = *T[res]; 
	Tqo = *T[arg1];
	T0[arg2]=cos(T0[arg1]);
	T0[res] = sin(T0[arg1]);
        for( l = 0; l < p; l++)
        {
          for (i=0;i<k-1;i++)
	  {
	    x=0;y=0;
	    z[i] = (i+1)*Tqo[i];
            x+= T0[arg2]*z[i];
            y+= T0[res]*z[i];
	    for (j=0;j<i;j++)
	      {
		x+= Targ[j]*z[i-j-1]; 
		y+= Tres[j]*z[i-j-1];
	      } /* endfor */
	    Targ[i] = -y/(i+1);
	    Tres[i]=x/(i+1);
	  } /* endfor */
          Tres += degree;
          Targ += degree;
          Tqo += degree;
        } /* endfor */
	break;
      case cos_op:  /* cosine operation */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
	Targ = *T[arg2];
	Tres = *T[res];
	Tqo = *T[arg1];
        T0[arg2]=sin(T0[arg1]);
        T0[res] = cos(T0[arg1]);
        for( l = 0; l < p; l++)
        {
	  for (i=0;i<k-1;i++)
	  {
	    x=0;y=0;
	    z[i] = (i+1)*Tqo[i];
            x+= T0[arg2]*z[i];
            y+= T0[res]*z[i];
	    for (j=0;j<i;j++)
	      {
		x+= Targ[j]*z[i-j-1]; 
		y+= Tres[j]*z[i-j-1];
	      } /* endfor */
	    Targ[i] = y/(i+1);
	    Tres[i] = -x/(i+1);
	  } /* endfor */
          Tres += degree;
          Targ += degree;
          Tqo += degree;
        } /* endfor */
	break;
      case asin_op: 
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        Targ = *T[arg2];
        Tres = *T[res];
        Tqo = *T[arg1];
	T0[res]=asin(T0[arg1]);
        for( l = 0; l < p; l++)
        {
	  for (i=0;i<k-1;i++)
	  {
	    x=0;
	    z[i] = (i+1)*Tqo[i];
            x+= T0[arg2]*z[i];
	    for (j=0;j<i;j++)
	      x+= Targ[j]*z[i-j-1]; 
	    Tres[i]=x/(i+1);
	  } /* endfor */
          Tres += degree;
          Targ += degree;
          Tqo += degree;
        } /* endfor */
	break;
      case acos_op:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        Targ = *T[arg2];
        Tres = *T[res];
        Tqo = *T[arg1];
	T0[res]=acos(T0[arg1]);
        for( l = 0; l < p; l++)
        {
	  for (i=0;i<k-1;i++)
	  {
	    x=0;
	    z[i] = (i+1)*Tqo[i];
            x+= T0[arg2]*z[i];
	    for (j=0;j<i;j++)
	      x+= Targ[j]*z[i-j-1];
	    Tres[i]=x/(i+1);
	  } /* endfor */
          Tres += degree;
          Targ += degree;
          Tqo += degree;
        } /* endfor */
	break;
      case atan_op:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        Targ = *T[arg2];
        Tres = *T[res];
        Tqo = *T[arg1];
	T0[res]=atan(T0[arg1]);
        for( l = 0; l < p; l++)
        {
	  for (i=0;i<k-1;i++)
	  {
	    x=0;
	    z[i] = (i+1)*Tqo[i];
            x+= T0[arg2]*z[i];
	    for (j=0;j<i;j++)
	      x+= Targ[j]*z[i-j-1]; 
	    Tres[i]=x/(i+1);
	  } /* endfor */
          Tres += degree;
          Targ += degree;
          Tqo += degree;
        } /* endfor */
	break;
      case gen_quad:  
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        Targ = *T[arg2];
        Tres = *T[res];
        Tqo = *T[arg1];
	T0[res]=coval;
        for( l = 0; l < p; l++)
        {
	  for (i=0;i<k-1;i++)
	  {
	    x=0;
	    z[i] = (i+1)*Tqo[i];
            x+= T0[arg2]*z[i];
	    for (j=0;j<i;j++)
	      x+= Targ[j]*z[i-j-1];
	    Tres[i]=x/(i+1);
	  } /* endfor */
          Tres += degree;
          Targ += degree;
          Tqo += degree;
        } /* endfor */
	break;
      case log_op:
        arg=get_locint_f();
        res= get_locint_f();
	Tres = *T[res];
	Targ = *T[arg];
	divs = 1.0/T0[arg];
	T0[res]=log(T0[arg]);
        for( l = 0; l < p; l++)
        {
	  for (i=0;i<k-1;i++)
	  {
	    x=0;
	    for (j=1;j<i+1;j++)
	      x+= Targ[i-j]*z[j-1]; 
	    Tres[i]=(Targ[i]-x/(i+1))*divs;
	    z[i] = (i+1)*Tres[i];
	  } /* endfor */
          Tres += degree;
          Targ += degree;
        } /* endfor */
	break;
      case pow_op:
        arg=get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
	Tres = *T[res];
	Targ = *T[arg];
	T0[res] = pow(T0[arg],coval);
	r0 = 1/T0[arg];
        for( l = 0; l < p; l++)
        {
          if (Targ[0]==0.0){
#ifdef inf_num
  InfVal=i_num/i_den;
  NoNum=n_num/n_den;
#endif
            r0=0.0;
            for (i=1;i<k;i++){
              if (Targ[i]>0.0){
                r0=InfVal;
                i=k;
              } /* end if */
              if (Targ[i]<0.0){
                r0=NoNum;
                i=k;
              } /* end if */
            } /* end for */
          } /* end if */
          else {
            r0 = 1.0/T0[arg];
          } /* end else */
	  for (i=1;i<k;i++)
	  {
	    x = 0 ;
	    y = coval*i;
            x += T0[res]*Targ[i-1]*y;
            y -= coval+1;
	    for (j=1;j<i;j++)
	      {
		x += Tres[j-1]*Targ[i-j-1]*y;
		y -= coval+1;
	      } /* endfor */
	    Tres[i-1] = x*r0/i;
	  } /* endfor */
          Tres += degree;
          Targ += degree;
        } /* endfor */
	break;

      
      case int_av_av:
        loc1_v = get_locint_f();
        size = get_locint_f();
        result_v= get_locint_f();
        for (l=0;l<size;l++)
        {
          res = result_v + l;   /* Location of left-hand-side  */
          arg = loc1_v + l;   /* Location of right-hand-side */
          /* code for int_adb_a */
          Targ = *T[arg];
          Tres = *T[res];
          T0[res]= T0[arg];
          for (pl=0;pl<p;pl++)
          {
            for (i=0;i<k-1;i++)
              Tres[i]=Targ[i];
            Tres += degree;
            Targ += degree;
          } /* endfor */
        } /* endfor */
        break;
      case assign_dv:
        size = get_locint_f();
        loc1_v = get_locint_f();
        d = get_val_v_f(size);
        for (l=0;l<size;l++)
        {
            res = loc1_v + l;           /* Location of left-hand-side */
            coval = d[l];      /* Value of right-hand-side   */     
            /* code for assign_d */
            Tres = *T[res];
            T0[res]=coval;
            for (pl=0;pl<p;pl++)
            {
              for (i=0;i<k-1;i++)
	        Tres[i]=0;
              Tres += degree;
            } /* endfor */
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
          /* code for assign_a */
          rloc = loc1;
          Targ = *T[arg];
          Tres = *T[res];
          T0[res]=T0[arg];
          for (pl=0;pl<p;pl++)
          {
            for (i=0;i<k-1;i++)
              Tres[i]=Targ[i];
            Tres += degree;
            Targ += degree;
          } /* endfor */
        } /* endfor */
        break;
      case assign_indvec:
        size = get_locint_f();
        loc1_v = get_locint_f();
        for (l=0;l<size;l++)
        {
          res = loc1_v + l;             /* Location of the left-hand-side */
          /* code for assign_ind */
          Tres = *T[res];
          T0[res]=basepoint[indexi];
          for (pl=0;pl<p;pl++)
          {
            for (i=0;i<k-1;i++){
              Tres[i]=argument[indexi][pl][i];}
            Tres += degree;
          } /* endfor */
          ++indexi;
        } /* endfor */
        break;
      case assign_depvec:
        size = get_locint_f();
        loc1_v = get_locint_f();
        for (l=0;l<size;l++)
        {
          res = loc1_v + l;             /* Location of the left-hand-side */
          /* code for assign_dep */
          Tres = *T[res];
          valuepoint[indexd]=T0[res];
          if (taylors != 0 )  
          {
            for (pl=0;pl<p;pl++)
            {
              for (i=0;i<k-1;i++)
	        taylors[indexd][pl][i] = Tres[i];
              Tres += degree;
            } /* endfor */
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
          /* code for eq_min_a */ 
          Tres = *T[res];
          Targ = *T[arg];
          T0[res]-=T0[arg];
          for (pl=0;pl<p;pl++)
          {
            for (i=0;i<k-1;i++)
              Tres[i] -= Targ[i];
            Tres+= degree;
            Targ+=degree;
          } /* endfor */
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
          /* code for eq_plus_a */
          Tres = *T[res];
          Targ = *T[arg];
          T0[res]+= T0[arg];
          for (pl=0;pl<p;pl++)
          {
            for (i=0;i<k-1;i++)
              Tres[i]+=Targ[i];
            Tres+= degree;
            Targ+= degree;
          } /* endfor */
        } /* endfor */
        break;
      case eq_mult_av_a:
        arg = get_locint_f();
        size = get_locint_f();
        result_v = get_locint_f();
        for (l=0;l<size;l++)
        {
          res = result_v + l;      /* Location of the left-hand-side  */
          /* code for eq_mult_a*/
          Tres = *T[res];
          Targ = *T[arg];
          for (pl=0;pl<p;pl++)
          {
            for (i=k-2;i>=0;i--)
            {
	      x=T0[res]*Targ[i];
	      for (j=0;j<i;j++)
	        x+=Tres[j]*Targ[i-j-1];
              x+=Tres[i]*T0[arg];
	      Tres[i]=x;
            } /* endfor */
            Tres+= degree;
            Targ+= degree;
          } /* endfor */
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
          /* code for eq_mult_d*/
          Tres = *T[res];
          T0[res]*=coval;
          for (pl=0;pl<p;pl++)
          {
            for (i=0;i<k-1;i++)
              Tres[i]*=coval;
            Tres+= degree;
          } /* endfor */
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
          /* code for plus_a_a */
          Tres = *T[res];
          Targ1 = *T[arg1];
          Targ2 = *T[arg2];
          T0[res]=T0[arg1]+T0[arg2];
          for (pl=0;pl<p;pl++)
          {
            for (i=0;i<k-1;i++)
              Tres[i]=Targ1[i]+Targ2[i];
            Tres+= degree;
            Targ1+= degree;
            Targ2+= degree;
          } /* endfor */
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
          /* code for min_a_a */
          Tres = *T[res];
          Targ1 = *T[arg1];
          Targ2 = *T[arg2];
          T0[res]=T0[arg1]-T0[arg2];
          for (pl=0;pl<p;pl++)
          {
            for (i=0;i<k-1;i++)
              Tres[i]=Targ1[i]-Targ2[i];
            Tres+= degree;
            Targ1+= degree;
            Targ2+= degree;
          } /* endfor */
        } /* endfor */
        break;
      case dot_av_av:
        loc1_v   = get_locint_f();
        loc2_v   = get_locint_f();
        size     = get_locint_f();
        res = get_locint_f();
        Tres = *T[res];
        T0[res]=0;
        for (pl=0;pl<p;pl++)
        {
          for (i=0;i<k-1;i++)
            Tres[i]=0.0;
          Tres+=degree;
        } /* endfor */
        for (l=0;l<size;l++)
        {
          arg2 = loc2_v + l;
          arg1 = loc1_v + l;
          /* code for mult_a_a  */
          Targ1 = *T[arg1];
          Targ2 = *T[arg2];
          Tres = *T[res];
          T0[res]+=T0[arg1]*T0[arg2];
          for (pl=0;pl<p;pl++)
          {
            for (i=0;i<k-1;i++)
            {
	      x=Tres[i];
              x+=T0[arg1]*Targ2[i];
	      for (j=0;j<i;j++)
	        x+=Targ1[j]*Targ2[i-j-1];
              x+=Targ1[i]*T0[arg2];
	      Tres[i]=x;
            } /* endfor */
            Tres+= degree;
            Targ1+= degree;
            Targ2+= degree;
          } /* endfor */
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
          /* code for mult_d_a */
          Tres = *T[res];
          Targ = *T[arg];
          T0[res]=T0[arg]*coval;
          for (pl=0;pl<p;pl++)
          {
            for (i=0;i<k-1;i++)
              Tres[i]=Targ[i]*coval;
            Tres+= degree;
            Targ+= degree;
          } /* endfor */
        } /* endfor */
        break;
      case div_av_a:
        loc1_v   = get_locint_f();
        arg2     = get_locint_f();
        size     = get_locint_f();
        result_v = get_locint_f();
        for (l=0;l<size;l++)
        {
          arg1 = loc1_v   + l;  /* Location of the right-hand-side vector[l] */
          res = result_v + l;  /* Location of the result */
          /* code for div_a_a */
          Tres = *T[res];
          Targ1 = *T[arg1];
          Targ2 = *T[arg2];
          divs   = 1/T0[arg2];
          T0[res]=T0[arg1]*divs;
          for (pl=0;pl<p;pl++)
          {
            for (i=0;i<k-1;i++)
            {
	      x=Targ1[i]*divs;
              x+=T0[res]*(-Targ2[i]*divs);
	      z[i]=-Targ2[i]*divs;
	      for (j=0;j<i;j++)
	        x+=Tres[j]*z[i-j-1];
	      Tres[i]=x;
            } /* endfor */
            Tres+= degree;
            Targ1+= degree;
            Targ2+= degree;
          } /* endfor */
        } /* endfor */
        break;
      case mult_av_a:
        loc1_v   = get_locint_f();
        arg2     = get_locint_f();
        size     = get_locint_f();
        result_v = get_locint_f();
        for (l=0;l<size;l++)
         {
          arg1 = loc1_v   + l; /* Location of the right-hand-side vector[l] */
          res = result_v + l; /* Location of the result */
          /* code for mult_a_a */
          Tres = *T[res];
          Targ1 = *T[arg1];
          Targ2 = *T[arg2];
          T0[res]=T0[arg1]*T0[arg2];
          for (pl=0;pl<p;pl++)
          {
            for (i=0;i<k-1;i++)
            {
	      x=T0[arg1]*Targ2[i];
	      for (j=0;j<i;j++)
	        x+=Targ1[j]*Targ2[i-j-1];
              x+=Targ1[i]*T0[arg2];
	      Tres[i]=x;
            } /* endfor */
            Tres+= degree;
            Targ1+= degree;
            Targ2+= degree;
          } /* endfor */
        } /* endfor */
        break;
      case mult_a_av:
        loc1_v   = get_locint_f();
        arg1     = get_locint_f();
        size     = get_locint_f();
        result_v = get_locint_f();
        for (l=0;l<size;l++)
        {
          arg2= loc1_v+l;      /* Location of right hand side vectore[l]  */
          res = result_v + l; /* Location of the result */
          /* code for mult_a_a */
          Tres = *T[res];
          Targ1 = *T[arg1];
          Targ2 = *T[arg2];
          T0[res]=T0[arg1]*T0[arg2];
          for (pl=0;pl<p;pl++)
          {
            for (i=0;i<k-1;i++)
            {
              x=T0[arg1]*Targ2[i];
              for (j=0;j<i;j++)
                x+=Targ1[j]*Targ2[i-j-1];
              x+=Targ1[i]*T0[arg2];
              Tres[i]=x;
            } /* endfor */
            Tres+= degree;
            Targ1+= degree;
            Targ2+= degree;
          } /* endfor */
        } /* endfor */
        break;

#ifdef conditional
      case cond_assign:
        arg = get_locint_f();
        loc1 = get_locint_f();
        loc2 = get_locint_f();
        res = get_locint_f();
        Tres = *T[res];
        Targ = *T[arg];
        Tloc1 = *T[loc1];
        Tloc2 = *T[loc2];
        if (T0[arg]>0)
        {
          T0[res]=T0[loc1];
        } /* endif */
        else
        {
          T0[res]=T0[loc2];
        } /* endif */
        for( l = 0; l < p; l++)
        {
          if (T0[arg]>0)
          {
            for (i=0;i<k-1;i++)
              Tres[i]=Tloc1[i];
          } /* endif */
          else
          {
            for (i=0;i<k-1;i++)
              Tres[i]=Tloc2[i];
          } /* endelse */
          Tres += degree;
          Tloc1+= degree;
          Tloc2+= degree;
        } /* endfor */
        break;
      case cond_assign_s:
        arg = get_locint_f();
        loc1 = get_locint_f();
        res = get_locint_f();
        Tres = *T[res];
        Targ = *T[arg];
        Tloc1 = *T[loc1];
        if (T0[arg]>0)
	{
          T0[res]=T0[loc1];
        } /* endif */
        for(l = 0; l < p; l++)
        {
          if (T0[arg]>0)
          {
            for (i=0;i<k-1;i++)
              Tres[i]=Tloc1[i];
          } /* endif */
          Tres += degree;
          Tloc1+= degree;
        } /* endfor */
        break;
      case subscript:
        loc1_v=get_locint_f();  /* Base */
        loc1=get_locint_f();/* pointer to the variable containing the offset */
        res=get_locint_f();
        Tres=*T[res];
        arg=loc1_v+(int)(T0[loc1]);
        Targ = *T[arg];
        T0[res]=T0[arg];
        for( l = 0; l < p; l++)
        { 
          for (i=0;i<k-1;i++)
            Tres[i]=Targ[i];
          Tres += degree;
          Targ += degree; 
        } /* endfor */
        break;
      case subscript_l:
        loc1_v=get_locint_f();  /* Base */
	loc1=get_locint_f();
        arg=loc1_v+(int)(T0[loc1]);
        res=get_locint_f();
        Tres = *T[res];
        Targ = *T[arg];
        T0[arg]=T0[res];
        for( l = 0; l < p; l++)
        {
          for (i=0;i<k-1;i++)
            Targ[i]=Tres[i];
         Tres += degree;
         Targ += degree; 
        } /* endfor */
        break;
      case subscript_ld:
        loc1_v=get_locint_f();  /* Base */
        loc1=get_locint_f();  /* pointer to the variable containing the offset */
        res = loc1_v+(int)(T0[loc1]);
        Tres = *T[res];
        T0[res]=get_val_f();
        for( l = 0; l < p; l++)
        {
          for (i=0;i<k-1;i++)
            Tres[i]=0;
          Tres += degree;
        } /* endfor */
        break;
      case m_subscript:
        loc1_v=get_locint_f();  /* Base */
        loc1=get_locint_f();/* pointer to the variable containing the offset */
        size=get_locint_f();
        result=get_locint_f();
        for (l=0;l<size;l++)
        {
          res=result+l;
          Tres=*T[res];
          arg=loc1_v+(int)(T0[loc1])*size+l;
          Targ = *T[arg];
          T0[res]=T0[arg];
          for(pl=0;pl<p;pl++)
          { 
            for (i=0;i<k-1;i++)
              Tres[i]=Targ[i];
            Tres += degree;
            Targ += degree; 
          } /* endfor */
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
          Tres = *T[res];
          Targ = *T[arg+l];
          T0[res]=T0[arg+l];
          for(pl=0;pl<p;pl++)
          {
            for (i=0;i<k-1;i++)
              Tres[i]=Targ[i];
           Tres += degree;
           Targ += degree; 
          } /* endfor */
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
          Tres = *T[res];
          T0[res]=d[l];
          for(pl=0;pl<p;pl++)
          {
            for (i=0;i<k-1;i++)
              Tres[i]=0;
            Tres += degree;
          } /* endfor */
        } /* endfor */
        break;
#endif

      case death_not:
        loc1=get_locint_f();
        loc1=get_locint_f();
  	break;
      default:
	/* Die here, we screwed up */
	printf("ADOL-C error: Fatal error in hov_forward for op %d\n",operation);
	break;
	
      } /* endswitch */
      operation=get_op_f();
    }  /* endwhile */
  end_sweep();
} /* end hov_forward */



#ifdef __cplusplus
}
#endif

