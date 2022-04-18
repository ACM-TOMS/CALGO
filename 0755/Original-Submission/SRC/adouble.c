/*
   --------------------------------------------------------------
   File adouble.c of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
   Adouble.c contains that definitions of procedures used to defined various
   badouble, adub, and adouble operations. These operations actually 
   have two purposes.
   The first purpose is to actual compute the function, just as the same
   code written for double precision (single precision - complex - interval) 
   arithmetic would.  The second purpose is to write a transcript of the
   computation for the reverse pass of automatic differentiation.
*/

/* Local Includes */

#include "adouble.h"
#include "oplate.h"

/* Include these routines that are written in straight C. */

extern "C" {
#include "taputil1.h"
}

/* Include Files */

#include <string.h>
#include <stdio.h>

void condassign(double &res, const double &cond, const double &arg1, const double &arg2)
{
 res=cond ? arg1 : arg2;
}

void condassign(double &res, const double &cond, const double &arg1)
{
  res=cond ? arg1 : res;
}

/* Global vars */

double* store;
int trace_flag =0;

static locint maxloc = sizeof(locint) ==2 ? 65535 : 2147483647;
static locint current_top = 0; // = largest live location + 1
static locint location_cnt = 0 ;    // = maximal # of lives so far
static locint maxtop = 0; // = current size of store
static locint maxtop2;
static locint dealloc = 0; // = # of locations to be freed
static locint deminloc = 0;  // = lowest loc to be freed 

/*------------------------------------------------------*/
/* The first several routines are for memory management */
/*------------------------------------------------------*/

/*
  Return the next free location in "adouble" memory
*/
locint next_loc()
{ 
/* First deallocate dead adoubles if they form a contiguous tail: */
#ifdef overwrite
  if (dealloc && dealloc+deminloc == current_top)    
    {          
     if(trace_flag) write_death(deminloc, current_top - 1);
     current_top = deminloc ;
     dealloc =0; 
     deminloc = maxloc;
     }
#endif
  if ( current_top == location_cnt) ++location_cnt ;
  if ( location_cnt > maxtop ) 
     {
      maxtop2 = ++maxtop*2 > maxloc ? maxloc : 2*maxtop;
      if (maxtop2 == maxloc )
	 {
         printf("ADOL-C fatal error  !! \n ");
         printf("maximal number of live active variables exceeded\n\n");
	 printf("Possible remedies :\n\n ");
       	 printf(" 1. Use more automatic local variables and \n");
	 printf("    allocate/deallocate adoubles on free store\n");
	 printf("     in a strictly last in first out fashion\n\n");
	 printf(" 2. Extend the range by redefining the type of \n");
         printf("    locint from unsigned short or int to int or long\n");
	 exit(-3);
	 }
      else
	 {
	 maxtop = maxtop2;
	 if(maxtop == 2)
	  {
	   store = (double *)malloc(maxtop*sizeof(double));
	   deminloc = maxloc;
	   }
         else
           {
	 store = (double *)realloc((char *)store,maxtop*sizeof(double));
	   }
	 if( store == 0) 
	   {
           printf("ADOL-C fatal error !! \n");
  	   printf("Failure to reallocate storage for adouble values\n");
	 printf("Possible remedies :\n\n ");
       	 printf(" 1. Use more automatic local variables and \n");
	 printf("    allocate/deallocate adoubles on free store\n");
	 printf("    in a strictly last in first out fashion\n");
	 printf(" 2. Extend the range by redefining the type of\n");
         printf("    locint to unsigned short or int to int or long\n");
	 printf(" 3. Enlarge your system stacksize limit\n");
	   exit(-3);
	   }
	}
      }
  return current_top++ ;
}


locint next_loc(int size)
{ 
/* First deallocate dead adoubles if they form a contiguous tail: */

#ifdef overwrite
  if (dealloc && dealloc+deminloc == current_top)    
    {          
     if(trace_flag) write_death(deminloc, current_top - 1);
     current_top = deminloc ;
     dealloc =0; 
     deminloc = maxloc;
     }
#endif
  if ( (current_top+size) >= location_cnt) location_cnt=current_top+size+1;
  while ( location_cnt > maxtop ) 
     {
      maxtop2 = ++maxtop*2 > maxloc ? maxloc : 2*maxtop;
      if (maxtop2 == maxloc )
	 {
         printf("ADOL-C fatal error !! \n ");
         printf("maximal number of live active variables exceeded\n\n");
	 printf("Possible remedies :\n\n ");
       	 printf(" 1. Use more automatic local variables and \n");
	 printf("    allocate/deallocate adoubles on free store\n");
	 printf("     in a strictly last in first out fashion\n\n");
	 printf(" 2. Extend the range by redefining the type of \n");
         printf("    locint from unsigned short or int to int or long\n");
	 exit(-3);
	 }
      else
	 {
	 maxtop = maxtop2;
	 if(maxtop == 2)
	  {
	   store = (double *)malloc(maxtop*sizeof(double));
	   deminloc = maxloc;
	   }
         else
           {
	     /* Allocate the storage */
	     double *temp;
	     temp = (double *)malloc(maxtop*sizeof(double));
	      
	     /* Copy over storage */
	     for (int i=0; i<current_top; i++)
	       temp[i]=store[i];

	      free((char*) store);
    
	      store = temp;
	   }
	 if( store == 0) 
	   {
           printf("ADOL-C fatal error !! \n");
  	   printf("Failure to reallocate storage for adouble values\n");
	 printf("Possible remedies :\n\n ");
       	 printf(" 1. Use more automatic local variables and \n");
	 printf("    allocate/deallocate adoubles on free store\n");
	 printf("    in a strictly last in first out fashion\n");
	 printf(" 2. Extend the range by redefining the type of\n");
         printf("    locint to unsigned short or int to int or long\n");
	 printf(" 3. Enlarge your system stacksize limit\n");
	   exit(-3);
	   }
	}
      }

#ifdef DEBUG
  printf ("ADOL-C debug: Top is: %d\n ",current_top+size);
#endif
  locint return_val=current_top;
  current_top+=size;
  return return_val;
}


/*
	Free a location in "adouble" memory.  
*/
inline void free_loc(locint old_loc)
{
  ++dealloc;
  if (old_loc < deminloc) deminloc = old_loc ;
}

void free_loc(int old_loc, int size)
{
  dealloc+=size;
  if (old_loc < deminloc) deminloc = old_loc ;
}
  

void take_stock() 
{ 
  for(int i =0; i< current_top; i++)
       { 
	double storei =store[i]; // Avoid I/O of NaN's !
        if (storei == storei) write_int_assign_d(i,storei);
	}
  trace_flag = 1;
}

locint keep_stock()
{
  if (current_top > 0) write_death(0,current_top - 1);
  trace_flag = 0;
  return location_cnt;
}


/*----------------------------------------------------------------*/
/* The remaining routines define the badouble,adub,and adouble    */
/* routines.                                                      */
/*----------------------------------------------------------------*/

/* Basic constructors */
/*
adub::adub(double y)
{
  location = next_loc();
  store[location] = y;
  if (trace_flag) write_int_assign_d(location,y);
}
*/
adouble::adouble()
{
  location = next_loc();
  
}

adouble::adouble(double y)
{
  location = next_loc();
  store[location] = y;
  if (trace_flag) write_int_assign_d(location,y);
}

adouble::adouble(const adouble& a)
{
  location = next_loc();
  store[location]=store[a.location];
  if (trace_flag) write_int_assign_a(location,a.location);
}


adouble::adouble(const adub& a)
{
  location = next_loc();
  store[location]=store[a.loc()];
  if (trace_flag) write_int_assign_a(location,a.loc());
}

/* Destructors */

#ifdef overwrite
adouble::~adouble()
{
   ++dealloc; 
   if (location < deminloc) deminloc = location;
}

adub::~adub()
{
  ++dealloc; 
  if (location < deminloc) deminloc = location;
}

#ifdef conditional
asub::~asub()
{
  ++dealloc;
  if (location < deminloc)
    deminloc = location;
}
#endif

#endif

/*
  Member function returns the location of this adouble.
*/
locint badouble::loc() const 
{
  return location;
}

/*
  double returns the true floating point value of an adouble variable
*/
double value(const badouble& x) 
{
  return store[x.location];
}

/*
  Define what it means to assign an adouble variable a constant value.
*/
badouble& badouble::operator = (double y) 
{  
  if (trace_flag) write_assign_d(location,y);
  store[location]=y;
  return *this;
}   

/*
  Define what it means to assign an adouble variable an independent value.
*/

badouble& badouble::operator <<= (double y) 
{  
  if (trace_flag) write_assign_ind(location); 
  store[location]=y;
  return *this;
}   

/*
  Define what it means to assign a float variable a dependent adouble value.
*/
badouble& badouble::operator >>= (double& y) 
{  
  if (trace_flag) write_assign_dep(location);
  y = double (store[location]);
  return *this;
}   

/*
  Define what it means to assign an Badouble variable an Badouble value.
  Optionally trace this action.
*/
badouble& badouble::operator = (const badouble& x) 
{  
  if (trace_flag) write_assign_a(location,x.location);
  store[location]=store[x.location];
  return *this;
}   
badouble& badouble::operator = (const adub& a)
{
  if (trace_flag) write_assign_a(location,a.location);
  store[location]=store[a.location] ;
  return *this;
}

/*
  Define define what it means to output an adouble value. 
  No tracing of this action
*/

ostream& operator << (ostream& out,const badouble& y)
{
  return out << store[y.location] << "(a)" ;
}

/*
  Define define what it means to input adouble value.
*/

istream& operator >> (istream& in,const badouble& y)
{
  double temp;
  in >> temp;
  store[y.location]=temp;
  if (trace_flag) write_assign_d(y.location,temp);
  return in;
}
  

adub adouble::operator++(int)  /* postfix increment */
{
  locint locat = next_loc();
  store[locat]=store[location];
  if (trace_flag) write_assign_a(locat,location);
  if (trace_flag) write_d_same_arg(eq_plus_d,location,1.0);
  store[location]++;
  return locat ;
}

adub adouble::operator--(int) /* postfix decrement */
{
  locint locat = next_loc();
  store[locat]=store[location];
  if (trace_flag) write_assign_a(locat,location);
  if (trace_flag) write_d_same_arg(eq_min_d,location,1.0);
  store[location]--;
  return locat ;
}

badouble& adouble::operator++() /* prefix increment */
{
  if (trace_flag) write_d_same_arg(eq_plus_d,location,1.0);
  store[location]++;
  return *this;
}

badouble& adouble::operator--() /* prefix decrement */
{
  if (trace_flag) write_d_same_arg(eq_min_d,location,1.0);
  store[location]--;
  return *this;
}


/*
  Adding a floating point to an adouble. Optionally trace this action.
*/
badouble& badouble::operator += (double y) 
{  
  if (trace_flag) write_d_same_arg(eq_plus_d,location,y);
  store[location]+=y;
  return *this; 
} 


/*
  Subtracting a floating point from an adouble. Optionally trace this
  activity.
*/
badouble& badouble::operator -= (double y) 
{  
  if (trace_flag) write_d_same_arg(eq_min_d,location,y);
  store[location]-=y;
  return *this;
}

/*
  Add an adouble to another adouble. Optionally trace this action.
*/
badouble& badouble::operator += (const badouble& y) 
{  
  if (trace_flag) write_a_same_arg(eq_plus_a,location,y.location);
  store[location]+=store[y.location];
  return *this;
}

/*
  Subtract an adouble from another adouble. Optionally trace this execution.
*/
badouble& badouble::operator -= (const badouble& y) 
{  
  if (trace_flag) write_a_same_arg(eq_min_a,location,y.location);
  store[location]-=store[y.location];
  return *this;
}

/* 
  Multiply an adouble by a float. Optionally trace this execution.
*/
badouble& badouble::operator *= (double y) 
{  
  if (trace_flag) write_d_same_arg(eq_mult_d,location,y);
  store[location]*=y;
  return *this;
}

/*
  Multiply one badouble by another. Optional trace.
*/
badouble& badouble::operator *= (const badouble& y) 
{ 
  if (trace_flag) write_a_same_arg(eq_mult_a,location,y.location);
  store[location]*=store[y.location];
  return *this;
}

badouble& badouble::operator /= (double y) 
{  
  *this = *this/y;
  return *this;
}

badouble& badouble::operator /= (const badouble& y) 
{ 
  *this = *this*(1.0/y);
  return *this;
}

/*   The Not Equal Operator (!=)      */

int operator != (const badouble& u,const badouble& v)
{
  return store[u.location] != store[v.location];
}

int operator != (double u,const badouble& v)
{
  return u != store[v.location];
}

int operator != (const badouble& v,double u)
{
  return store[v.location] != u;
}

/*   The Equal Operator (==)      */

int operator == (const badouble& u,const badouble& v)
{
  return store[u.location] == store[v.location];
}

int operator == (double u,const badouble& v)
{
  return u == store[v.location];
}

int operator == (const badouble& v,double u)
{
  return store[v.location] == u;
}

/*   The Less than or Equal Operator (<=)      */

int operator <= (const badouble& u,const badouble& v)
{
  return store[u.location] <= store[v.location];
}

int operator <= (double u,const badouble& v)
{
  return u <= store[v.location];
}

int operator <= (const badouble& v,double u)
{
  return store[v.location] <= u;
}

/*   The Greater than or Equal Operator (>=)      */

int operator >= (const badouble& u,const badouble& v)
{
  return store[u.location] >= store[v.location];
}

int operator >= (double u,const badouble& v)
{
  return u >= store[v.location];
}

int operator >= (const badouble& v,double u)
{
  return store[v.location] >= u;
}

/*   The Greater than Operator (>)      */

int operator > (const badouble& u,const badouble& v)
{
  return store[u.location] > store[v.location];
}

int operator > (double u,const badouble& v)
{
  return u > store[v.location];
}

int operator > (const badouble& v,double u)
{
  return store[v.location] > u;
}

/*   The Less than Operator (<)      */

int operator < (const badouble& u,const badouble& v)
{
  return store[u.location] < store[v.location];
}

int operator < (double u,const badouble& v)
{
  return u < store[v.location];
}

int operator < (const badouble& v,double u)
{
  return store[v.location] < u;
}

/* 
  Adding two badoubles.  NOTE: calculates address of temporary, and returns
  an adub.
*/

adub operator + (const badouble& x, const badouble& y) 
{ 
  locint locat = next_loc();
  store[locat]= store[x.location]+store[y.location];
  if (trace_flag) write_two_a_rec(plus_a_a,locat,x.location,y.location);
  return locat;
}

/*
  Adding a badouble and a double. Optional trace. Temporary assignment.
*/

adub operator + (double x, const badouble& y) 
{ 
  locint locat = next_loc();
  store[locat]= x+store[y.location];
  if (trace_flag) write_args_d_a(plus_d_a,locat,x,y.location);
  return locat;
}

adub operator + (const badouble& y, double x) 
{ 
  locint locat = next_loc();
  store[locat]= x+store[y.location];
  if (trace_flag) write_args_d_a(plus_d_a,locat,x,y.location);
  return locat;
}

/*
  Subtraction of two badoubles. Optional Trace. Temporary used.
*/

adub operator - (const badouble& x, const badouble& y)
{  
  locint locat = next_loc();
  store[locat]=store[x.location]-store[y.location];
  if (trace_flag) write_two_a_rec(min_a_a,locat,x.location,y.location);
  return locat;
}

/*
  Subtract a badouble from a double. Optional trace. Temporary used. 
*/

adub operator - (double x, const badouble& y)
{
  locint locat = next_loc();
  store[locat]=x-store[y.location];
  if (trace_flag) write_args_d_a(min_d_a,locat,x,y.location);
  return locat; 
}


/*
  Multiply two badouble numbers.  Optional trace. Use temporary.
*/

adub operator * (const badouble& x, const badouble& y)
{ 
  locint locat = next_loc();
  store[locat]=store[x.location]*store[y.location];
  if (trace_flag) write_two_a_rec(mult_a_a,locat,x.location,y.location);
  return locat;
}

/* 
  Multiply a badouble by a double.  Optional Trace. 
*/ 

adub operator * (double x, const badouble& y)
{ 
  locint locat = next_loc();
  store[locat]=x*store[y.location];
  if (trace_flag) write_args_d_a(mult_d_a,locat,x,y.location);
  return locat;
}

/*
   Divide a badouble by an badouble.
*/

adub operator / (const badouble& x, const badouble& y)
{  
  locint locat = next_loc();
  store[locat] = store[x.location]/store[y.location];
  if (trace_flag) write_two_a_rec(div_a_a,locat,x.location,y.location);
  return locat;
}

/*
  Division double - badouble.
*/

adub operator / (double x, const badouble& y)
{ 
  locint locat = next_loc();
  store[locat]= x/store[y.location];
  if (trace_flag) write_args_d_a(div_d_a,locat,x,y.location);
  return locat;
}

/* 
  Compute exponential of badouble.  
*/

adub exp (const badouble& x)
{  
  locint locat = next_loc();
  store[locat]=exp(store[x.location]);
  if (trace_flag) write_single_op(exp_op,locat,x.location);
  return locat; 
}

/*
  Compute logarithm of badouble. Optional Trace. Use temporary. 
*/

adub log (const badouble& x)
{ 
  locint locat = next_loc();
  store[locat]=log(store[x.location]);
  if (trace_flag) write_single_op(log_op,locat,x.location);
  return locat; 
}

/*
  Compute sqrt of adouble. Optional Trace. Use temporary. 
*/

adub sqrt (const badouble& x)
{ 
  locint locat = next_loc();
  store[locat]=sqrt(store[x.location]);
  if (trace_flag) write_single_op(sqrt_op,locat,x.location);
  return locat; 
}

/*
  Compute sin of badouble. Optional trace.
  Note:Sin and Cos are always evaluated together
*/
adub sin (const badouble& x) 
{ 
  locint locat = next_loc(); 
  store[locat]=sin(store[x.location]);
  adouble y;
  store[y.location]=cos(store[x.location]);
  if (trace_flag) write_quad(sin_op,locat,x.location,y.location);
  return locat;
}

/*
  Compute cos of badouble. Optional trace.
*/

adub cos (const badouble& x)
{ 
  locint locat = next_loc();
  store[locat]=cos(store[x.location]);
  adouble y;
  store[y.location]=sin(store[x.location]);
  if (trace_flag) write_quad(cos_op, locat,x.location,y.location);
  return locat; 
}

adub tan (const badouble& x) 
{ 
  return sin(x)/cos(x);
}

/*
  Asin value. -- really a quadrature. Use temporary. Optional trace.
*/

adub asin (const badouble& x)
{
  locint locat = next_loc();
  adouble y = 1.0/sqrt(1.0-x*x);
  store[locat]=asin(store[x.location]);
  if (trace_flag) write_quad(asin_op,locat,x.location,y.location);
  return locat; 
}

/*
  Acos value. -- really a quadrature. Use temporary. Optional trace.
*/

adub acos (const badouble& x)
{
  locint locat = next_loc();
  adouble y = -1.0/sqrt(1.0-x*x);
  store[locat]=acos(store[x.location]);
  if (trace_flag) write_quad(acos_op,locat,x.location,y.location);
  return locat; 
}

/*
  Atan value. -- really a quadrature. Use temporary. Optional trace.
*/

adub atan (const badouble& x)
{
  locint locat = next_loc();
  adouble y= 1.0/(1.0+x*x);
  store[locat]=atan(store[x.location]);
  if (trace_flag) write_quad(atan_op,locat,x.location,y.location);
  return locat; 
}

adub atan2(const badouble& y,const badouble& x)
{
  const double pihalf = asin(1.0);
  double vx = value(x);
  double vy = value(y);
  double sy = (vy > 0) ? 1.0 : -1.0 ;
  if(fabs(vx) > fabs(vy))
    {
      if(vx>0)
	return atan(y/x);
      else
	return atan(y/x)+sy*2*pihalf;
    } 
  else
    {
      if(vy !=0) 
	return sy*pihalf-atan(x/y);
      else 
      {  /* nodifferentiable case */
        locint locat=next_loc();
        store[locat]=0.0;
        if(trace_flag) write_int_assign_d(locat,0.0);
	return locat;  
      }
    }
}

/*
  power value. -- adouble/double .
*/

adub pow (const badouble& x,double y)
{
  locint locat = next_loc();
  store[locat] = pow(store[x.location],y);
  if (trace_flag) write_args_d_a(pow_op,locat,y,x.location);
  return locat;
}

/*
  power value. --- adouble/adouble .
*/

adub pow (const badouble& x,const badouble& y)  
{
  double vx = store[x.location];
  if(vx>0)
    return exp(y*log(x));
  else
    {
      double vy = store[y.location];
      if(vx < 0 || vy >= 0)
	{
	  cout << "ADOL-C message: exponent of negative basis deactivated \n";
	  return pow(x,value(y));
	}
      else 
	cout << "ADOL-C message: negative exponent and zero basis deactivated \n";
      locint locat=next_loc();
      store[locat]=pow(vx,vy);
      if(trace_flag) write_int_assign_d(locat,store[locat]);
      return locat;
    }   
}

/*
  log base 10 of an adouble value.
*/

adub log10 (const badouble& x) 
{
  return log(x)/log(10);
}

/*
  Hyperbolic Sine of an adouble value.
*/

adub sinh (const badouble& x) 
{
  adouble temp= exp(x);
  return  0.5*(temp-1/temp);
}

/*
  Hyperbolic Cosine of an adouble value.
*/

adub cosh (const badouble& x) 
{
  adouble temp= exp(x);
  return 0.5*(temp+1/temp);
}

/*
  Hyperbolic Tangent of an adouble value.
*/

adub tanh (const badouble& x) 
{
  adouble temp= exp(2*x);
  return (temp-1)/(temp+1);
}

/*
  Ceiling Function (Note: This function is nondifferentiable)
*/

adub ceil (const badouble& x) 
{
  locint locat=next_loc();
  store[locat]=ceil(store[x.location]);
  if(trace_flag) write_int_assign_d(locat,store[locat]);
  return locat;
}

/*
  Floor Function (Note: This function is nondifferentiable)
*/

adub floor (const badouble& x) 
{
  locint locat=next_loc();
  store[locat]=floor(store[x.location]);
  if(trace_flag) write_int_assign_d(locat,store[locat]);
  return locat;
}

#ifdef Inverse_hyperbolics

/*
  Asinh value. -- really a quadrature. Use temporary. Optional trace.
*/

adub asinh (const badouble& x)
{
  locint locat = next_loc();
  adouble y= 1.0/sqrt(1.0+x*x);
  store[locat]=asinh(store[x.location]);
  if (trace_flag) write_quad(gen_quad,locat,x.location,y.location);
  return locat; 
} 

/*
  Acosh value. -- really a quadrature. Use temporary. Optional trace.
*/

adub acosh (const badouble& x)
{
  locint locat = next_loc();
  adouble y= 1.0/sqrt(1.0-x*x);
  store[locat]=acosh(store[x.location]);
  if (trace_flag) write_quad(gen_quad,locat,x.location,y.location);
  return locat; 
}

/*
  Atanh value. -- really a quadrature. Use temporary. Optional trace.
*/

adub atanh (const badouble& x)
{
  locint locat = next_loc();
  adouble y= 1.0/(1.0-x*x);
  store[locat]=atanh(store[x.location]);
  if (trace_flag) write_quad(gen_quad,locat,x.location,y.location);
  return locat; 
}
#endif

/*
  Fabs Function (Note: This function is also nondifferentiable at x=0)
*/
adub fabs (const badouble& x)
{
  locint locat = next_loc();
  store[locat] = fabs(store[x.location]);
  if (trace_flag) write_single_op(abs_val,locat,x.location);
  return locat;
}

/*
  max and min functions 
*/

adub max (const badouble& x, const badouble& y)
{
  return (0.5*(x+y+fabs(x-y)));
}

adub min (const badouble& x, const badouble& y)
{
  return (0.5*(x+y-fabs(x-y)));
}


/*
  Ldexp Function. 
*/

adub ldexp (const badouble& x,int exp) 
{
  return x*ldexp(1.0,exp);
}


/*  The error function erf, enable if your math.h contains erf(double) */

adub erf(const badouble& x) 
{
  locint locat = next_loc();
  adouble q= exp(-x*x);
  store[locat]=erf(store[x.location]);
  if (trace_flag) write_quad(gen_quad,locat,x.location,q.location);
  return locat;
} 

/* Macro for user defined quadratures, example myquad is below */

#define extend_quad(func,integrand)\
adouble func (const badouble& arg)\
{  adouble temp; \
    adouble val; \
    integrand; \
    store[temp.location]=func(store[arg.location]); \
    if (trace_flag) \
	write_quad(gen_quad,temp.location,arg.location,val.location);\
    return temp; }

double myquad(double& x)
{
  double res;
  res = log(x);
  return res;
}

/* This defines the natural logarithm as a quadrature */

extend_quad(myquad,val = 1/arg)


/* ADDITIONAL ASSIGNMENTS */

/*
  Define what it means to assign an adouble variable a float value.
*/

adouble& adouble::operator = (double y) 
{  
   if (trace_flag) write_assign_d(location,y);
   store[location]=y;
   return *this;
}  

adouble& adouble::operator = (const badouble& x) 
{  
   if (trace_flag) write_assign_a(location,x.loc());
   store[location]=store[x.loc()];
   return *this;
}   

adouble& adouble::operator = (const adouble& x) 
{  
   if (trace_flag) write_assign_a(location,x.location);
   store[location]=store[x.location];
   return *this;
}   

badouble& badouble::operator = (const adouble& x) 
{  
   if (trace_flag) write_assign_a(location,x.location);
   store[location]=store[x.location];
   return *this;
}   

adouble& adouble::operator = (const adub& a)
{
  if (trace_flag) write_assign_a(location,a.loc());
  store[location]=store[a.loc()] ;
  return *this;
}


#ifdef conditional
void condassign(adouble &result, const adouble &arg1, const adouble &r1, const adouble &r2) 
{
  if (trace_flag)
    write_condassign(result.location,arg1.location,r1.location,
		     r2.location);

  if (store[arg1.location]>0)
    store[result.location]=store[r1.location];
  else
    store[result.location]=store[r2.location];
}

void condassign(adouble &result, const adouble &arg1, const adouble &r1) 
{
  if (trace_flag)		
    write_condassign2(result.location,arg1.location,r1.location);
  
  if (store[arg1.location]>0)
    store[result.location]=store[r1.location];
}

asub::asub(locint start, locint index)
{
#ifdef DEBUG
  printf("ADOL-C debug: Constructing an asub with 2 arguments\n");
#endif
  
  base=start;
  offset=index;

  location=next_loc();

  store[location]=store[base+(int)store[offset]];
  
  if (trace_flag)
    write_associating_value(subscript,location,base,offset);
}

along& along::operator = (int y) 
{  
   if (trace_flag) write_assign_d(location,y);
   store[location]=y;
   return *this;
}   

asub& asub::operator <<= (double y) 
{  
  if (trace_flag)
    {
      write_assign_ind(location);
      write_associating_value(subscript_l,location,base,offset);
    }

  store[base+(int)store[offset]]=y;

  return *this;
}   

asub& asub::operator = (const adub& a)
{
  if (trace_flag)
      write_associating_value(subscript_l,a.loc(),base,offset);
  store[base+(int)store[offset]]=store[a.loc()] ;
  return *this;
}

asub& asub::operator = (double x)
{
  if (trace_flag)
    write_associating_value_ld(subscript_ld,x,base,offset);
    store[base+(int)store[offset]]=x;
  return *this;
}

asub& asub::operator = (const badouble& x) 
{  
  if (trace_flag)
    write_associating_value(subscript_l,x.loc(),base,offset);
  store[base+(int)store[offset]]=store[x.loc()];
  return *this;
}   

along& along::operator = (const badouble& x) 
{  
  if (trace_flag)
    write_assign_a(location,x.loc());
  store[location]=store[x.loc()];
  return *this;
}

along& along::operator = (const along& x) 
{  
  if (trace_flag)
    write_assign_a(location,x.location);
  store[location]=store[x.location];
  return *this;
}

along& along::operator = (const adub& a)
{
  if (trace_flag)
    write_assign_a(location,a.loc());
  store[location]=store[a.loc()] ;
  return *this;
}

along::along(int y)
{
  store[location] = y;
  if (trace_flag) write_int_assign_d(location,y);
}

along::along(const along& a)
{
  store[location]=store[a.location];
  if (trace_flag) write_int_assign_a(location,a.location);
}

along::along(const adub& a)
{
  store[location]=store[a.loc()];
  if (trace_flag) write_int_assign_a(location,a.loc());
}
#endif



/* end of adouble.c */


