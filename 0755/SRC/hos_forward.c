/*
   ----------------------------------------------------------------
   File hos_forward.c of ADOL-C version 1.6 as of January 1,   1995
   ----------------------------------------------------------------
   Contains the routine hos_forward (higher-order-scalar forward 
   mode).

*/

#ifdef __cplusplus
extern "C" {
#endif

#include "dvlparms.h" /* Developers Parameters */

/* Necessary Includes */

#include "usrparms.h"
#include "oplate.h"
#include "taputil1.h"
#include "taputil2.h"
#include "taputil3.h"
#include "tayutil.h"


/* Local Static Variables */

static double** T;

static short tag;

static int for_location_cnt;
static int dep_cnt;
static int ind_cnt;
static int degree;
static int deaths;

/****************************************************************************/
/* Higher Order Scalar version of the forward mode.                         */
/****************************************************************************/
void hos_forward(short tnum,         /* tape id */
		 int depcheck,       /* consistency chk on # of dependents */
		 int indcheck,       /* consistency chk on # of independents */
		 int gdegree,        /* highest derivative degree */
		 int keep,           /* flag for reverse sweep */
		 double **argument,  /* independant variable values */
		 double **taylors)   /* matrix of coifficient vectors */
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

  int k, i, j ,l, loop,rloc;
  double r0;
  locint arg,res;
  int indexi =0; 
  int indexd =0;
  double  *Targ,*Tres,*Tqo,*Trloc2,*Targ1,*Targ2; 
  double coval,divs;
  int buffer; 
  double x,y;
  static int fax, kax;
  static double *Tdum, *z; 
  int taylbuf;
  int even;
  double* Td;

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
  deaths = tape_stats[3];
  buffer =  tape_stats[4];

  set_buf_size(buffer);


  
  if ((depcheck != dep_cnt)||(indcheck != ind_cnt))
    {
      printf("\n ADOL-C error: forward sweep on tape %d  aborted!\n",tag);
      printf("Number of dependent and/or independent variables passed to forward is\ninconsistant with number recorded on tape %d \n",tag);
      exit (-1);
    }
  
   if (k compsize kax || for_location_cnt compsize fax)
    {
      free((char*) Tdum);
      free((char **) T);
      free((char *) z); 
      Tdum = (double*) malloc(for_location_cnt*k*sizeof(double));
      z = (double *)malloc(sizeof(double)*k);
      Td=Tdum;
      T = (double **)malloc(sizeof(double*)*for_location_cnt);
      for (i=0;i<for_location_cnt;i++)
	{
	  T[i]=Td;
	  Td += k;
	}
      kax= k;
      fax = for_location_cnt;
    }
  taylbuf = keep*sizeof(revreal)*deaths;
  if (keep)
    {
      taylbuf = 4+taylbuf/(1+taylbuf/bufsize);
      taylor_begin(taylbuf,T,keep-1);
    }
  
  /* Initialize the Forward Sweep */
  init_for_sweep(tag);

  operation=get_op_f();
  while (operation !=end_of_tape)
    { /* Switch Statement to perform operations */
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
	
	/* Scalar Operations */

      case int_adb_a:	
        Targ = T[get_locint_f()];
        Tres = T[get_locint_f()];
	for (i=0;i<k;i++)
	  Tres[i]=Targ[i];
	break;
      case assign_a:	
        Targ = T[get_locint_f()];
        res = get_locint_f();
        Tres = T[res];
	if (keep)
	  write_taylor(res,keep);
	for (i=0;i<k;i++)
	  Tres[i]=Targ[i];
	break;
      case int_adb_d:
	Tres = T[get_locint_f()];
	Tres[0]=get_val_f();
	for (i=1;i<k;i++)
	  Tres[i]=0;
	break;
      case assign_d:
        res = get_locint_f();
	Tres = T[res];
	if (keep)
	  write_taylor(res,keep);
	Tres[0]=get_val_f();
	for (i=1;i<k;i++)
	  Tres[i]=0;
	break;
      case assign_ind:
        res = get_locint_f();
	Tres = T[res];
	if (keep)
	  write_taylor(res,keep);
	for (i=0;i<k;i++){
	  Tres[i]=argument[indexi][i];}
	++indexi;
	break;
      case assign_dep:
	Tres = T[get_locint_f()];
	if (taylors != 0 )  
	  for (loop=0;loop<k;loop++)
	    taylors[indexd][loop] = Tres[loop];
	indexd++;
	break;
      case eq_plus_d:
        res=get_locint_f();
	Tres = T[res];
	if (keep)
	  write_taylor(res,keep);
	Tres[0]+= get_val_f();
	break;
      case eq_plus_a: 
	Targ = T[get_locint_f()];
        res=get_locint_f();
	Tres = T[res];
	if (keep)
	  write_taylor(res,keep);
	for (i=0;i<k;i++)
	  Tres[i]+=Targ[i];
	break;
      case eq_min_d:
	res=get_locint_f();
	Tres = T[res];
	if (keep)
	  write_taylor(res,keep);
	Tres[0]-=get_val_f();
	break;
      case eq_min_a:
	Targ = T[get_locint_f()];
        res=get_locint_f();
	Tres = T[res];
	if (keep)
	  write_taylor(res,keep); 
	for (i=0;i<k;i++)
	  Tres[i] -= Targ[i];
	break;
      case eq_mult_d: 
        res=get_locint_f();
        coval = get_val_f();
	Tres = T[res];
	if (keep)
	  write_taylor(res,keep); 
	for (i=0;i<k;i++)
	  Tres[i]*=coval;
	break;
      case eq_mult_a: /* Convolution */
	Targ = T[get_locint_f()];
        res=get_locint_f();
	Tres = T[res];
	if (keep)
	  write_taylor(res,keep); 
	for (i=k-1;i>=0;i--)
	  {
	    x=0;
	    for (j=0;j<=i;j++)
	      x+=Tres[j]*Targ[i-j];
	    Tres[i]=x;
	  }
	break;
      case plus_a_a:
	Targ1 = T[get_locint_f()];
	Targ2 = T[get_locint_f()];
	Tres = T[get_locint_f()];
	for (i=0;i<k;i++)
	  Tres[i]=Targ1[i]+Targ2[i];
	break;
      case plus_d_a:
	Targ = T[get_locint_f()];
	Tres = T[get_locint_f()];
	Tres[0]=Targ[0]+ get_val_f();
	for (i=1;i<k;i++)
	  Tres[i]=Targ[i];
	break;
      case min_a_a:
        Targ1 = T[get_locint_f()];
        Targ2 = T[get_locint_f()];
        Tres = T[get_locint_f()];
	for (i=0;i<k;i++)
	  Tres[i]=Targ1[i]-Targ2[i];
	break;
      case min_d_a:
        Targ = T[get_locint_f()];
        Tres = T[get_locint_f()];
	Tres[0]=get_val_f() -Targ[0];
	for (i=1;i<k;i++)
	  Tres[i] = -Targ[i];
	break;
      case mult_a_a:
        Targ1 = T[get_locint_f()];
        Targ2 = T[get_locint_f()];
        Tres = T[get_locint_f()];
	for (i=0;i<k;i++)
	  {
	    x=0;
	    for (j=0;j<=i;j++)
	      x+=Targ1[j]*Targ2[i-j];
	    Tres[i]=x;
	  }
	break;
      case mult_d_a:
        Targ = T[get_locint_f()];
        Tres = T[get_locint_f()];
	coval =get_val_f();
	for (i=0;i<k;i++)
	  Tres[i]=Targ[i]*coval;
	break;
      case div_a_a:
        Targ1 = T[get_locint_f()];
        Targ2 = T[get_locint_f()];
        Tres = T[get_locint_f()];
	divs   = 1/Targ2[0];
	for (i=0;i<k;i++)
	  {
	    x=Targ1[i]*divs;
	    z[i]=-Targ2[i]*divs;
	    for (j=0;j<i;j++)
	      x+=Tres[j]*z[i-j];
	    Tres[i]=x;
	  }
	break;
      case div_d_a:
        Targ = T[get_locint_f()];
        Tres = T[get_locint_f()];
        coval =get_val_f();
	divs = 1/Targ[0];
	Tres[0]=coval*divs;
	for (i=1;i<k;i++)
	  {
	    x=0;
	    z[i]=-Targ[i]*divs;
	    for (j=0;j<i;j++)
	      x+=Tres[j]*z[i-j];
	    Tres[i] = x;
	  }
	break;
      case sqrt_op:
        Targ = T[get_locint_f()];
	Tres = T[get_locint_f()];
	Tres[0]=sqrt(Targ[0]);
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
          r0 = 0.5/Tres[0];
        } /* end else */
        even =0;
        for (i=1;i<k;i++) {
          x=0;
          for (j=1;2*j<i;j++)
            x+=Tres[j]*Tres[i-j];
          x *= 2;
          if(even) x += Tres[i/2]*Tres[i/2];
          even = !even;
          Tres[i]=r0*(Targ[i]-x);
        } /* end for */
	break;
      case abs_val:
        Targ = T[get_locint_f()];
        Tres = T[get_locint_f()];
        x=0.0;
        for (i=0;i<k;i++)
        {
            if((x==0.0) && (Targ[i]!=0.0))
            {
              if (Targ[i]<0.0)
                x = -1.0;
              else
                x=1.0;
            }
            Tres[i] = x*Targ[i];
        }
        break;

      case exp_op:
	Targ = T[get_locint_f()];
	Tres = T[get_locint_f()];
	Tres[0]=exp(Targ[0]);
	for (i=1;i<k;i++)
	  { double x=0;

	    z[i] = i*Targ[i];
	    for (j=0;j<i;j++)
	      x+= Tres[j]*z[i-j];
	    Tres[i]=x/i;
	  }
	break;
      case sin_op:
        coval=get_val_f();
        Tqo = T[get_locint_f()];
        Targ = T[get_locint_f()];
        Tres = T[get_locint_f()];
	Targ[0]=cos(Tqo[0]);
	Tres[0] = sin(Tqo[0]);
	for (i=1;i<k;i++)
	  {
	    x=0;y=0;
	    z[i] = i*Tqo[i];
	    for (j=0;j<i;j++)
	      {
		x+= Targ[j]*z[i-j]; 
		y+= Tres[j]*z[i-j];
	      }
	    Targ[i] = -y/i;
	    Tres[i]=x/i;
	  }
	break;
      case cos_op: 
        coval=get_val_f();
        Tqo = T[get_locint_f()];
        Targ = T[get_locint_f()];
        Tres = T[get_locint_f()];
	Targ[0]=sin(Tqo[0]);
	Tres[0]=cos(Tqo[0]);
	for (i=1;i<k;i++)
	  {
	    x=0;y=0;
	    z[i] = i*Tqo[i];
	    for (j=0;j<i;j++)
	      {
		x+= Targ[j]*z[i-j]; 
		y+= Tres[j]*z[i-j];
	      }
	    Targ[i] = y/i;
	    Tres[i] = -x/i;
	  }
	break;
      case asin_op:
	coval=get_val_f();
        Tqo = T[get_locint_f()];
        Targ = T[get_locint_f()];
        Tres = T[get_locint_f()];
        Tres[0]=asin(Tqo[0]);
	for (i=1;i<k;i++)
	  {
	    x=0;
	    z[i] = i*Tqo[i];
	    for (j=0;j<i;j++)
	      x+= Targ[j]*z[i-j]; 
	    Tres[i]=x/i;
	  }
	break;
      case acos_op:
        coval=get_val_f();
        Tqo = T[get_locint_f()];
        Targ = T[get_locint_f()];
        Tres = T[get_locint_f()];
        Tres[0]=acos(Tqo[0]);
	for (i=1;i<k;i++)
	  {
	    x=0;
	    z[i] = i*Tqo[i];
	    for (j=0;j<i;j++)
	      x+= Targ[j]*z[i-j];
	    Tres[i]=x/i;
	  }
	break;
      case atan_op:
        coval=get_val_f();
        Tqo = T[get_locint_f()];
        Targ = T[get_locint_f()];
        Tres = T[get_locint_f()];
	Tres[0]=atan(Tqo[0]);
	for (i=1;i<k;i++)
	  {
	    x=0;
	    z[i] = i*Tqo[i];
	    for (j=0;j<i;j++)
	      x+= Targ[j]*z[i-j]; 
	    Tres[i]=x/i;
	  }
	break;
      case gen_quad:
        Tqo = T[get_locint_f()];
        Targ = T[get_locint_f()];
        Tres = T[get_locint_f()];
	Tres[0]=get_val_f();
	for (i=1;i<k;i++)
	  {
	    x=0;
	    z[i] = i*Tqo[i];
	    for (j=0;j<i;j++)
	      x+= Targ[j]*z[i-j];
	    Tres[i]=x/i;
	  }
	break;
      case log_op:
	Targ = T[get_locint_f()];
	Tres = T[get_locint_f()];
	divs = 1.0/Targ[0];
	Tres[0]=log(Targ[0]);
	for (i=1;i<k;i++)
	  {
	    x=0;
	    for (j=1;j<i;j++)
	      x+= Targ[i-j]*z[j]; 
	    Tres[i]=(Targ[i]-x/i)*divs;
	    z[i] = i*Tres[i];
	  }
	break;
      case pow_op:
        coval = get_val_f();
	Targ = T[get_locint_f()];
	Tres = T[get_locint_f()];
	Tres[0] = pow(Targ[0],coval);
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
          r0 = 1.0/Targ[0];
        } /* end else */
	for (i=1;i<k;i++)
	  {
	    x = 0 ;
	    y = coval*i;
	    for (j=0;j<i;j++)
	      {
		x += Tres[j]*Targ[i-j]*y;
		y -= coval+1;
	      }
	    Tres[i] = x*r0/i;
	  }
	break;
	
      case death_not:
        loc1 = get_locint_f();
        loc2 = get_locint_f();
        if (keep) 
        {
          do write_taylor(loc2,keep);
          while(loc1 < loc2-- ); 
        }
	break;
	
      /* Vector Operations */
      
      case int_av_av:
        loc1_v = get_locint_f();
        size = get_locint_f();
        result_v= get_locint_f();
	for (l=0;l<size;l++)
	  {
	    Tres=T[result_v + l];   /* Location of left-hand-side  */
	    Targ=T[loc1_v + l];     /* Location of right-hand-side */
	    /* code for int_adb_a */
	    for (i=0;i<k;i++)
	      Tres[i]=Targ[i];
	  }
	break;
      case assign_dv:
        size = get_locint_f();
        loc1_v = get_locint_f();
        d = get_val_v_f(size);
	for (l=0;l<size;l++)
	  {
	    res = loc1_v + l;           /* Location of left-hand-side */
	    Tres = T[res];
	    if (keep)
	      write_taylor(res,keep);
	    Tres[0]= d[l];
	    for (i=1;i<k;i++)
	      Tres[i]=0;
	  }
	break;
      case assign_av:
        loc2_v = get_locint_f();
        size = get_locint_f();
        loc1_v = get_locint_f();
	for (l=0;l<size;l++)
	  {
	    res = loc1_v + l;             /* Location of left-hand-side  */
	    Targ= T[loc2_v + l];  
            Tres=T[res];
	    /* code for assign_a */
	    if (keep)
	      write_taylor(res,keep);
	    for (i=0;i<k;i++)
	      Tres[i]=Targ[i];
	  }
	break;
      case assign_indvec:
        size = get_locint_f();
        loc1_v = get_locint_f();
	for (l=0;l<size;l++)
	  {
	    res = loc1_v + l;             /* Location of the left-hand-side */
	    /* code for assign_ind */
	    Tres = T[res];
	    if (keep)
	      write_taylor(res,keep);
	    for (i=0;i<k;i++){
	      Tres[i]=argument[indexi][i];}
	    ++indexi;
	  }
	break;
      case assign_depvec:
        size = get_locint_f();
        loc1_v = get_locint_f();
	for (l=0;l<size;l++)
	  {
	    Tres=T[loc1_v + l];             /* Location of the left-hand-side */
	    /* code for assign_dep */
	    if (taylors != 0 )  
	      for (loop=0;loop<k;loop++)
		taylors[indexd][loop] = Tres[loop];
	    indexd++;
	  }
	break;
      case eq_min_av:
        loc1_v = get_locint_f();
        size = get_locint_f();
        result_v = get_locint_f();
	for (l=0;l<size;l++)
	  {
	    res = result_v + l;            /* Location of left-hand-side  */
	    /* code for eq_min_a */ 
	    Tres = T[res];
	    Targ = T[loc1_v + l];
	    if (keep)
	      write_taylor(res,keep); 
	    for (i=0;i<k;i++)
	      Tres[i] -= Targ[i];
	  }
	break;
      case eq_plus_av:
        loc1_v = get_locint_f();
        size = get_locint_f();
        result_v = get_locint_f();
	for (l=0;l<size;l++)
	  {
	    res = result_v + l;            /* Location of left-hand-side  */
	    /* code for eq_plus_a */
	    Tres = T[res];
	    Targ = T[loc1_v   + l];
	    if (keep)
	      write_taylor(res,keep);
	    for (i=0;i<k;i++)
	      Tres[i]+=Targ[i];
	  }
	break;
      case eq_mult_av_a:
        Targ = T[get_locint_f()];
        size = get_locint_f();
        result_v = get_locint_f();
	for (l=0;l<size;l++)
	  {
	    res = result_v + l;      /* Location of the left-hand-side  */
	    /* code for eq_mult_a*/
	    Tres = T[res];
	    if (keep)
	      write_taylor(res,keep); 
	    for (i=k-1;i>=0;i--)
	      {
		x=0;
		for (j=0;j<=i;j++)
		  x+=Tres[j]*Targ[i-j];
		Tres[i]=x;
	      }
	  }
	break;
      case eq_mult_av_d:
        size = get_locint_f();
        result_v = get_locint_f();
        coval = get_val_f();
	for (l=0;l<size;l++)
	  {
	    res = result_v + l;      /* Location of the left-hand-side  */
	    /* stored_val = fixed;         Location on the right-hand-side */
	    Tres = T[res];
	    if (keep)
	      write_taylor(res,keep); 
	    for (i=0;i<k;i++)
	      Tres[i]*=coval;
	  }
	break;
      case plus_av_av:
        loc1_v   = get_locint_f();
        loc2_v   = get_locint_f();
        size     = get_locint_f();
        result_v = get_locint_f();
	for (l=0;l<size;l++)
	  {
	    Targ2= T[loc2_v    + l];       /* Location of var 2  */
	    Targ1 = T[loc1_v    + l];       /* Location of var 1  */
	    Tres = T[result_v  + l];       /* Location of result */
	    for (i=0;i<k;i++)
	      Tres[i]=Targ1[i]+Targ2[i];
	  }
	break;
      case sub_av_av:
        loc1_v   = get_locint_f();
        loc2_v   = get_locint_f();
        size     = get_locint_f();
        result_v = get_locint_f();
	for (l=0;l<size;l++)
	  {
            Targ2= T[loc2_v    + l];       /* Location of var 2  */
            Targ1 = T[loc1_v    + l];       /* Location of var 1  */
            Tres = T[result_v  + l];       /* Location of result */
	    for (i=0;i<k;i++)
	      Tres[i]=Targ1[i]-Targ2[i];
	  }
	break;
      case dot_av_av:
        loc1_v   = get_locint_f();
        loc2_v   = get_locint_f();
        size     = get_locint_f();
        Tres = T[get_locint_f()];
	for (i=0;i<k;i++)
	  Tres[i]=0.0;
	for (l=0;l<size;l++)
	  {
	    Targ1 = T[loc2_v + l];
	    Targ2 = T[loc1_v + l];
	    /* code for mult_a_a  */
	    for (i=0;i<k;i++)
	      {
		x=Tres[i];
		for (j=0;j<=i;j++)
		  x+=Targ1[j]*Targ2[i-j];
		Tres[i]=x;
	      }
	  }
	break;
      case mult_d_av:
        loc1_v = get_locint_f();
        size = get_locint_f();
        result_v = get_locint_f();
        coval =get_val_f();
	for (l=0;l<size;l++)
	  {
	    Targ=T[loc1_v   + l];   /* Location on the right-hand-side */
	    Tres =T[result_v + l];   /* location of the result */
	    for (i=0;i<k;i++)
	      Tres[i]=Targ[i]*coval;
	  }
	break;
      case div_av_a:
        loc1_v   = get_locint_f();
        Targ2 =T[ get_locint_f()];
        size     = get_locint_f();
        result_v = get_locint_f();
	for (l=0;l<size;l++)
	  {
	    Targ1= T[loc1_v   + l];  /* Location of right-hand-side vector[l] */
	    Tres = T[result_v + l];  /* Location of the result */
	    divs   = 1/Targ2[0];
	    for (i=0;i<k;i++)
	      {
		x=Targ1[i]*divs;
		z[i]=-Targ2[i]*divs;
		for (j=0;j<i;j++)
		  x+=Tres[j]*z[i-j];
		Tres[i]=x;
	      }
	  }
	break;
      case mult_av_a:
        loc1_v   = get_locint_f();
        Targ2=T[ get_locint_f()];
        size     = get_locint_f();
        result_v = get_locint_f();
	for (l=0;l<size;l++)
	  {
	    Targ1=T[loc1_v   + l]; /* Location of right-hand-side vector[l] */
	    Tres=T[ result_v + l]; /* Location of the result */
	    for (i=0;i<k;i++)
	      {
		x=0;
		for (j=0;j<=i;j++)
		  x+=Targ1[j]*Targ2[i-j];
		Tres[i]=x;
	      }
	  }
	break;
      case mult_a_av:
        loc1_v   = get_locint_f();
        Targ1 = T[get_locint_f()];
        size     = get_locint_f();
        result_v = get_locint_f();
	for (l=0;l<size;l++)
	  {
	    Targ2=T[loc1_v+l];   /* Location of right hand side vector[l]  */
	    Tres=T[result_v + l]; /* Location of the result */
	    for (i=0;i<k;i++)
	      {
		x=0;
		for (j=0;j<=i;j++)
		  x+=Targ1[j]*Targ2[i-j];
		Tres[i]=x;
	      }
	  }
	break;

#ifdef conditional
      case cond_assign:
        Targ=T[get_locint_f()];
        Targ1=T[get_locint_f()];
        Targ2=T[get_locint_f()];
        res=get_locint_f();
        Tres=T[res];
        if (keep)
          write_taylor(res,keep);
	if (*Targ>0)
	  {
	    for (i=0;i<k;i++)
	      Tres[i]=Targ1[i];
	  }
	else
	  {
	    for (i=0;i<k;i++)
	      Tres[i]=Targ2[i];
	  }
	break;
      case cond_assign_s:
        Targ =T[get_locint_f()];
        Targ1 =T[get_locint_f()];
        res = get_locint_f();
        Tres=T[res];
        if (keep)
          write_taylor(res,keep);
	if (*Targ>0)
	{
	  for (i=0;i<k;i++)
	    Tres[i]=Targ1[i];
	}
	break;
      case subscript:	
        loc1_v=get_locint_f();  /* Base */
        loc1=get_locint_f();  /* Offset */
        result=get_locint_f();
	rloc = result;
	Trloc2 = T[loc1_v+(int)T[loc1][0]];
	if (keep)
	  write_taylor(rloc,keep);
	for (i=0;i<k;i++)
	  T[rloc][i]=Trloc2[i];
	break;
      case subscript_l:	
        loc1_v=get_locint_f();  /* Base */
        loc1=get_locint_f();  /* Offset */
        result=get_locint_f();
	rloc = loc1_v+(int)T[loc1][0];
	Trloc2 = T[result];
	if (keep)
	  write_taylor(rloc,keep);
	for (i=0;i<k;i++)
	  T[rloc][i]=Trloc2[i];
	break;
      case subscript_ld:	
        stored_val = get_val_f();
        loc1_v=get_locint_f();  /* Base */
        loc1=get_locint_f();  /* Offset */
	res = loc1_v+(int)T[loc1][0];
	Tres = T[res];
	coval = stored_val;
	Tres[0]=coval;
	for (i=1;i<k;i++)
	  Tres[i]=0;
	break;
      case m_subscript:	
        loc1_v=get_locint_f();  /* Base */
        loc1=get_locint_f();  /* Offset */
        size=get_locint_f();
        result=get_locint_f();
        for (l=0;l<size;l++)
        {
          arg=loc1_v+(int)T[loc1][0]*size+l;
          res=result+l;
	  Targ = T[arg];
          Tres = T[res];
	  for (i=0;i<k;i++)
	    Tres[i]=Targ[i];
        }
	break;
      case m_subscript_l:	
        loc1_v=get_locint_f();  /* Base LHS */
        loc1=get_locint_f();  /* Offset LHS */
        size=get_locint_f();
        arg=get_locint_f(); /* RHS */
        for (l=0;l<size;l++)
        {
          res = loc1_v+(int)T[loc1][0]*size+l;
	  Targ = T[arg+l];
          Tres = T[res];
	  if (keep)
	    write_taylor(res,keep);
	  for (i=0;i<k;i++)
	    Tres[i]=Targ[i];
        }
	break;
      case m_subscript_ld:	
        loc1_v=get_locint_f();  /* Base */
        loc1=get_locint_f();  /* Offset */
        loc2=get_locint_f(); /* offset in the vector itself */
        size=get_locint_f();
        d = get_val_v_f(size);
        for (l=0;l<size;l++)
        {
	  res = loc1_v+(int)T[loc1][0]*size+l+loc2;
	  Tres = T[res];
	  Tres[0]=d[l]; 
	  for (i=1;i<k;i++)
	    Tres[i]=0;
        }
	break;
#endif	

	
      default:
	/* Die here, we screwed up */
	printf("ADOL-C error: Fatal error in hos_forward with operation %d\n",
	       operation);
	exit(-1);
	break;
	
      } /* end of switch */

      /* Read the next operation */

      operation=get_op_f();
    } /* end of while */

  if (keep) taylor_close(taylbuf,depcheck,indcheck);

  end_sweep();

} /* end of hos_forward */



#ifdef __cplusplus
}
#endif

