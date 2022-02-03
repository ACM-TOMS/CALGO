/*
  -----------------------------------------------------------------------
  file avector.c of ADOL-C version 1.6 as of January 1,   1995            
  -----------------------------------------------------------------------
  Avector.c contains the necessary routines for vector operations       
  that are defined in avector.h.  Note: avector.h is already included   
  in adouble.h, and hence does not need to be included here.            
*/

/* Local Includes */

#include "adouble.h"
#include "oplate.h"

/* The "write" routines for the vector operations are written in straight C */
extern "C" {
#include "taputil1.h"
}

/* Extra Include Files */

#include <string.h>
#include <stdio.h>
#include <stdarg.h>

/* Global vars and routines from adouble.c */

extern double* store;
extern int trace_flag;
extern locint next_loc(int size);
extern locint next_loc();
extern locint free_loc(int,int);

/* ----- Start of vector operations ----- */

/* ----- ACTIVE VECTOR SECTION ----- */
  
adoublev::adoublev(int n)
{
#ifdef DEBUG
  printf("ADOL-C debug:Declaring active vector\n");
#endif
  
  size=n;
  start_loc=next_loc(size); 
}



adoublev::adoublev(const adoublev &op1)
{
#ifdef DEBUG
  printf("ADOL-C debug:Declaring active vector and initializing from adoublev\n");
#endif

  size=op1.size;
  start_loc=next_loc(size);

  for (int i=0; i<size; i++)
    store[start_loc+i]=store[op1.start_loc+i];
  
  if (trace_flag) write_intvec_assign_av(size,start_loc,op1.start_loc);
    
}

adoublev::adoublev(const adubv& a)
{
  start_loc = next_loc(a.sz());
  for (int i=0; i<a.sz(); i++)
    store[start_loc+i]=store[a.loc()+i];
  if (trace_flag) write_intvec_assign_av(a.sz(),start_loc,a.loc());
}


adoublev& adoublev::operator = (double* y) 
{
  if (trace_flag) write_assign_vec_dv(size,start_loc,y);

  for (int i=0; i<size; i++)
    store[start_loc+i]=y[i];
  
  return *this;
}


adoublev& adoublev::operator = (double y)
{
#ifdef DEBUG
  printf("ADOL-C debug:In adoublev=double\n");
#endif


  for (int i=0; i<size; i++)
      store[start_loc+i]=y;

  if (trace_flag)
    write_assign_vec_dv(size,start_loc,store+start_loc);
  return *this;
}



adoublev& adoublev::operator = (const badoublev& x) 
{  
  if (trace_flag) write_assign_av(size,start_loc,x.loc());
  for (int i=0; i<size; i++)
    store[start_loc+i]=store[x.loc()+i];
  return *this;
}

adoublev& adoublev::operator = (const adoublev& x) 
{  
  if (trace_flag) write_assign_av(size,start_loc,x.start_loc);
  for (int i=0; i<size; i++)
    store[start_loc+i]=store[x.start_loc+i];
  return *this;
}

adoublev& adoublev::operator = (const adubv& a)
{
  if (trace_flag) write_assign_av(size,start_loc,a.loc());
  for (int i=0; i<size; i++)
    store[start_loc+i]=store[a.loc()+i];
  return *this;
}

#ifdef overwrite
 adoublev::~adoublev()
{
#ifdef DEBUG
  printf("ADOL-C debug:Destructing active vector\n");
#endif
 
  free_loc(start_loc,size);
}

adubv::~adubv()
{
  free_loc(start_loc,size);
}

adoublem::~adoublem()
{
#ifdef DEBUG
  printf("ADOL-C debug:Destructing active matrix\n");
#endif
  delete[] index;
}
#endif

ostream& operator << (ostream& out,const badoublev &op1)
{
  out << "(";

  for (int i=0; i<op1.size-1; i++)
    out << store[op1.start_loc+i] << ", ";

  out << store[(op1.start_loc+op1.size)-1] << ")(a)";

  return out;
}  


badouble badoublev::operator[](int i) const  
{
  /* Used so can access the vector like an array with the [] */
  /* Check if out of range */
  if (i<0 || i>=size)
    {
      printf ("ADOL-C error: adoublev index out of range.\n");
      exit(-3);
    }

  return start_loc+i;
}

badoublev& badoublev::operator=(const badoublev &op1)
{
#ifdef DEBUG
  printf("ADOL-C debug:In badoublev = badoublev\n");
#endif

  if (trace_flag)
    write_assign_av(size,start_loc,op1.start_loc);
  
  for (int i=0; i<op1.size; i++)
    store[start_loc+i]=store[op1.start_loc+i];
  
  return *this;
}

badoublev& badoublev::operator= (const adubv &op1)
{
  if (trace_flag)
    write_assign_av(size,start_loc,op1.start_loc);
  
  for (int i=0; i<op1.size; i++)
    store[start_loc+i]=store[op1.start_loc+i];

  return *this;
}

badoublev& badoublev::operator = (const adoublev& x) 
{  
  if (trace_flag) write_assign_av(size,start_loc,x.start_loc);

  for (int i=0; i<size; i++)
    store[start_loc+i]=store[x.start_loc+i];
  
  return *this;
}   



/*
Define what it means to assign an adouble vector an independent float vector.
*/
adoublev& adoublev::operator <<= (double* y) 
{
#ifdef DEBUG
  printf("ADOL-C debug:IND EQ double*\n");
#endif

  if (trace_flag) write_assign_indvec(size,start_loc,y);
  for (int i=0; i<size; i++)
    store[(start_loc)+i]=y[i];
  return *this;
}   

/*
  Define what it means to assign a float vector a dependent adouble vector.
*/
adoublev& adoublev::operator >>= (double* y) 
{  
#ifdef DEBUG
  printf("ADOL-C debug:DEP EQ double* operator\n");
#endif

  if (trace_flag) write_assign_depvec(size,start_loc);
  for (int i=0; i<size; i++)
    y[i] = double (store[(start_loc)+i]);
  return *this;
}   

/* New operators here */




badoublev& badoublev::operator -= (const badoublev& y) 
{ 
  if (trace_flag)
    write_av_same_arg(eq_min_av,size,start_loc,y.start_loc);
  
  for (int i=0; i<size; i++)
    store[start_loc+i]-=store[y.start_loc+i];
  
  return *this;
}



badoublev& badoublev::operator += (const badoublev& y) 
{ 
  if (trace_flag)
    write_av_same_arg(eq_plus_av,size,start_loc,y.start_loc);
  
  for (int i=0; i<size; i++)
    store[start_loc+i]+=store[y.start_loc+i];
  
  return *this;
}

badoublev& badoublev::operator *= (double y) 
{
  if (trace_flag)
    write_samearg_av_d(eq_mult_av_d,size,start_loc,y);

  for (int i=0; i<size; i++)
    store[start_loc+i]*=y;

  return *this;
}

badoublev& badoublev::operator *= (const badouble& y) 
{
  int loc1,i;
  loc1 = y.loc();

  if (trace_flag)
    write_av_same_arg(eq_mult_av_a,size,start_loc,loc1);
  for (i=0; i<size; i++)
    store[start_loc+i]*=store[loc1];
  return *this;
}

badoublev& badoublev::operator/=(double y) 
{  
  *this = *this/y;
  return *this;
}

badoublev& badoublev::operator/=(const badouble& y) 
{ 
 *this = *this*(1.0/y);
  return *this;
}

adubv operator+(const badoublev &op1,const badoublev &op2)
{
  int i;
  locint size;
  locint start_loc;
  if (op1.size!=op2.size)
    {
      printf("ADOL-C error: Can not add vectors as not same size\n");
      exit(-3);
    }
  else
    {
      size=op1.size;
      start_loc=next_loc(size);
      for (i=0; i<size; i++)
	store[start_loc+i]=store[op1.start_loc+i]+store[op2.start_loc+i];
      
      if (trace_flag) write_two_av_rec(plus_av_av,size,start_loc,
				       op1.start_loc,op2.start_loc);
      
    }
  return adubv(start_loc,size);
}



adubv operator*(const badoublev &op1, double n)
{
  int i;
  locint size=op1.size;
  locint start_loc=next_loc(size);

  for (i=0; i<size; i++)
    store[start_loc+i]=store[op1.start_loc+i]*n;
  
  if (trace_flag) write_args_d_av(mult_d_av,size,start_loc,n,op1.start_loc);

  return adubv(start_loc,size);
}

adubv operator*(double n, const badoublev &op1)
{
  int i;
  locint size=op1.size;
  locint start_loc=next_loc(size);

  for (i=0; i<size; i++)
    store[start_loc+i]=store[op1.start_loc+i]*n;
  
  if (trace_flag)
    write_args_d_av(mult_d_av,size,start_loc,n,op1.start_loc);
    
  return adubv(start_loc,size);
}

adub operator*(const badoublev &op1,const badoublev &op2)
{
  int i;
  double x=0;
  locint locat=next_loc();
  
  if (op1.size!=op2.size)
    {
      printf("ADOL-C error: Can not take dot product, vectors are not same size\n");
      exit(-3);
    }
  else
    {
      for (i=0; i<op1.size; i++)
	x+=store[op1.start_loc+i]*store[op2.start_loc+i];
      
      store[locat]=x;
      
      if (trace_flag)
	write_two_av_rec(dot_av_av,op1.size,locat,
			 op1.start_loc,op2.start_loc);
      
    }
  return locat;
}

adubv operator / (const badoublev &x, const badouble &y)
{
  int i;
  int loc1 = y.loc();
  int size=x.size;
  locint start_loc=next_loc(size);

  for (i=0; i<size; i++)
    store[start_loc+i]=store[x.start_loc+i]*(1.0/store[loc1]);

  if (trace_flag)
    write_av_a_rec(div_av_a,size,start_loc,x.start_loc,loc1);
  
  return adubv(start_loc,size);
}

adubv operator-(const badoublev &op1,const badoublev &op2)
{
  int i;
  locint size;
  locint start_loc;
  if (op1.size!=op2.size)
    {
      printf("ADOL-C error: Can not add vectors as not same size\n");
      exit(-3);
    }
  else
    {
      size=op1.size;
      start_loc=next_loc(size);
      for (i=0; i<size; i++)
	store[start_loc+i]=store[op1.start_loc+i]-store[op2.start_loc+i];
      
      if (trace_flag) write_two_av_rec(sub_av_av,size,start_loc,
				       op1.start_loc,op2.start_loc);
      
    }
  return adubv(start_loc,size);
}




  
adubv operator*(const badoublev &op1, const badouble &n)
{
  int loc1 = n.loc();
  int size=op1.size;
  locint start_loc=next_loc(size);
  int i;

  for (i=0; i<size; i++)
    store[start_loc+i]=store[op1.start_loc+i]*store[loc1];

  if (trace_flag) write_av_a_rec(mult_av_a,size,start_loc,op1.start_loc,loc1);
    
  return adubv(start_loc,size);
}

adubv operator*(const badouble &n, const badoublev &op1)
{
  int loc1 = n.loc();
  int size=op1.size;
  locint start_loc=next_loc(size);
  int i;

  for (i=0; i<size; i++)
    store[start_loc+i]=store[loc1]*store[op1.start_loc+i];

  if (trace_flag)
    write_av_a_rec(mult_a_av,size,start_loc,op1.start_loc,loc1);

  return adubv(start_loc,size);
}

/* Active matrix section */
adoublem::adoublem(int row, int col)
{
  m=row;
  n=col;
  index = new adoublev[m];
  for (int i=0; i<m; i++)
  {
    index[i].size=n;
    index[i].start_loc=next_loc(n);
  }
}

adoublev& adoublem::operator[](int i) 
{
  if (i<0 || i>=m)
    {
      printf ("ADOL-C error: adoublem index out of range.\n");
      exit(-3);
    }
  return (index[i]);
}



#ifdef conditional
asub badoublev::operator[](const along &i)  const
{
#ifdef DEBUG
  printf("ADOL-C debug:In along overloaded []\n");
#endif
  
  /* Used so can access the vector like an array with the [] */
  /* Check if out of range */
  if ((i<0) || (i>=size))
    printf ("ADOL-C warning:: adoublev index out of range.\n");

  return asub(start_loc,i.loc());
}



/* asubv operations <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

#ifdef overwrite
asubv::~asubv()
{
#ifdef DEBUG
  printf("ADOL-C debug:Destructing active subscript vector\n");
#endif

  free_loc(start_loc,size);

};
#endif


asubv::asubv(adoublev* start, locint index)
{
#ifdef DEBUG
  printf("ADOL-C debug: Constructing an asubv with 3 arguments\n");
#endif
  begin=(start[0]).loc(); /* start of matrix */
  base=(start[(int)store[index]]).loc(); /* start of the i-th row */
  offset=index;
  size=(start[(int)store[index]]).sz(); /* size of the row-vector */
  start_loc=next_loc(size);
  for(int i=0;i<size;i++)
    store[start_loc+i]=store[base+i];
  if (trace_flag)
    write_associating_vector(m_subscript,start_loc,begin,offset,size);
}


asubv adoublem::operator[](const along &i)
{
#ifdef DEBUG
  printf("ADOL-C debug: In along overloaded []\n");
#endif
  /* Used so can access the vector like an array with the [] */
  /* Check if out of range */
  if (i<0 || i>=n)
    printf ("ADOL-C warning:: adoublem index out of range.\n");
  return asubv(index,i.loc());
}


asubv& asubv::operator = (const adubv& a)
{
  if (trace_flag)
    write_associating_vector(m_subscript_l,a.loc(),begin,offset,size);
  for(int i=0;i<size;i++)
    store[base+i]=store[a.loc()+i] ;
  return *this;
}

asubv& asubv::operator = (const badoublev& x)
{ 
  if (trace_flag)
    write_associating_vector(m_subscript_l,x.loc(),begin,offset,size);
  for(int i=0;i<size;i++)
    store[base+i]=store[x.loc()+i];
  return *this;
}  

asubv& asubv::operator <<= (double* y)
{
  if (trace_flag)
    {
      write_assign_indvec(size,start_loc,y);
      write_associating_vector(m_subscript_l,start_loc,begin,offset,size);
    }
  for(int i=0;i<size;i++)
    store[base+i]=y[i];
  return *this;
}

asubv& asubv::operator = (double* x)
{
  if (trace_flag)
    write_associating_vector_ld(x,begin,offset,size);
  for(int i=0;i<size;i++)
    store[base+i]=x[i];
  return *this;
}



#endif




adubv operator+ (const badoublev& x)
{
  return x*(1.0);
}
  
adubv operator- (const badoublev& x)
{
  return x*(-1.0);
}
