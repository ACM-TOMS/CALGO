/*
  ---------------------------------------------------------------
  file adouble.h of ADOL-C version 1.6 as of January 1,   1995 

  Included in:
              adouble.c
              avector.c
	      all ADOL-C applications programs.
  
  -----------------------------------------------------------------

  adouble.h contains the basis for the class of adouble
  included here are all the possible functions defined on
  the adouble class.  Notice that, as opposed to ealier versions,
  both the class adub and the class adouble are derived from a base
  class (badouble).  See below for further explanation.
*/

/*      D I S C L A I M E R 

The ADOL-C Software is provided under the following disclaimer:

NO WARRANTY.  The software was created in the course of a research
endeavor. It is not a commercial package.  The present version is
still in development, and is distributed "AS IS, WITH ALL DEFECTS."
By using the software, each user agrees to assume all responsibility
for any and all such use.  The authors and Argonne National Laboratory
are not aware that the software or the use thereof infringe any
proprietary right belonging to a third party.  However, NO WARRANTY,
CONDITION, OR REPRESENTATION OF ANY KIND, EXPRESS OR IMPLIED, is made
about the software, including without limitation any warranty of title,
noninfringement, merchantability, or fitness for a particular purpose,
by the authors or their affiliated institutions.

NO CONSEQUENTIAL DAMAGES.  Independent of the foregoing disclaimer
of warranties, each person that uses the software thereby agrees, that
NEITHER ARGONNE NATIONAL LABORATORY NOR THE AUTHORS OR THEIR AFFILIATED
INSTITUTIONS SHALL BE LIABLE FOR ANY INCIDENTAL OR CONSEQUENTIAL DAMAGES
IN CONNECTION WITH THE USE OF THE SOFTWARE, INCLUDING WITHOUT LIMITATION
LOST PROFITS OR INJURY TO BUSINESS, WHETHER OR NOT ARGONNE NATIONAL
LABORATORY, AND THE AUTHORS AND THEIR AFFILIATED INSTITUTIONS KNOW OR
HAVE REASON TO KNOW OF THE POSSIBILITY OF SUCH DAMAGES.

INDEMNITY.  Each person that uses the software thereby agrees, to
indemnify and defend Argonne National Laboratory and the authors and
their affiliated institutions, or any of them, against any loss, expense,
claim, damage, or liability of any kind arising from or connected with
their respective uses of the software, and to hold them or any of them
harmless from any of the same, WHETHER OR NOT ARISING IN WHOLE OR IN PART
FROM THE NEGLIGENCE OR GROSS NEGLIGENCE OF ARGONNE NATIONAL LABORATORY OR
ANY OF THE AUTHORS OR THEIR AFFILIATED INSTITUTIONS.

SUPPORT. Each person that uses this software understands that the software
is not supported by the authors or by their affiliated institutions.
*/


#include "usrparms.h"
#include "dvlparms.h" /* Developers Parameters */

#include <stream.h>

/*
  Notice that the purpose of the class adub is merely to avoid the generation
  and recording of an extra return adouble for each elementary 
  operation and function call. The same result can be achieved much
  more elegantly with GNUs named return variables, which would also 
  achieve the desired last in first out pattern for adouble construction 
  and destruction.
*/

class adouble;
class adub;
class badouble;
#ifdef conditional
class badoublev;
class adoublev;
class adubv;
/* class doublev; */


void condassign(double &res, const double &cond, const double &arg1, const double &arg2);
void condassign(double &res, const double &cond, const double &arg1);

inline double max(const double &x, const double &y){
  return (0.5*(x+y+fabs(x-y))); }
inline double min(const double &x, const double &y){
  return (0.5*(x+y-fabs(x-y))); }



#endif

/* 
   The class badouble contains the basic definitions for 
   the arithmetic operations, comparisons, etc. 
   This is a basic class from which the adub and adouble are 
   derived.  Notice that the constructors/destructors for 
   the class badouble are of the trivial variety.  This is the
   main difference among badoubles, adubs, and adoubles.
*/

class badouble{
  friend class badoublev;
 protected:
  locint location;
  badouble(){};
  badouble(const badouble& a){location = a.location;};
  badouble(locint lo){location = lo;};

 public:

#ifdef conditional

  friend void condassign(adouble &result, const adouble &arg1,
                         const adouble &r1, const adouble &r2);
  friend void condassign(adouble &result, const adouble &arg1,
                         const adouble &r1);
#endif

  locint loc() const;
  friend double value(const badouble&);
  badouble& operator >>= (double&);
  badouble& operator <<= (double);
  badouble& operator = (double);
  badouble& operator = (const badouble&);
  badouble& operator = (const adub&);
  badouble& operator = (const adouble&);
  badouble& operator += (double);
  badouble& operator += (const badouble&);
  badouble& operator -= (double y);
  badouble& operator -= (const badouble&);
  badouble& operator *= (double);
  badouble& operator *= (const badouble&);
  badouble& operator /= (double);
  badouble& operator /= (const badouble&);

  friend ostream& operator << (ostream&, const badouble&);
  friend istream& operator >> (istream&, const badouble&);

  friend int operator != (const badouble&,const badouble&);
  friend int operator != (double,const badouble&);
  friend int operator != (const badouble&,double);
  friend int operator == (const badouble&,const badouble&);
  friend int operator == (double,const badouble&);
  friend int operator == (const badouble&,double);
  friend int operator >= (const badouble&,const badouble&);
  friend int operator >= (double,const badouble&);
  friend int operator >= (const badouble&,double);
  friend int operator <= (const badouble&,const badouble&);
  friend int operator <= (double,const badouble&);
  friend int operator <= (const badouble&,double);
  friend int operator > (const badouble&,const badouble&);
  friend int operator > (double,const badouble&);
  friend int operator > (const badouble&,double);
  friend int operator < (const badouble&,const badouble&);
  friend int operator < (double,const badouble&);
  friend int operator < (const badouble&,double);

  /* End of Comparision Operators */

  inline friend adub operator + (const badouble& x); //{return x + 0.0 ;} ; 
  friend adub operator + (const badouble&,const badouble&); 
  friend adub operator + (double, const badouble&); 
  friend adub operator + (const badouble&, double); 
  inline friend adub operator - (const badouble& x ,double y); //{return (-y)+x;}; 
  friend adub operator - (const badouble&,const badouble&); 
  friend adub operator - (double, const badouble&); 
  inline friend adub operator - (const badouble& x); //{return 0.0 - x;}; 
  friend adub operator * (const badouble&,const badouble&); 
  friend adub operator * (double, const badouble& ); 
  inline friend adub operator * (const badouble& x, double y); //{return y*x;}; 
  inline friend adub operator / (const badouble& x, double y); // {return (1.0/y)*x;}; 
  friend adub operator / (const badouble&,const badouble&); 
  friend adub operator / (double,const badouble&); 
  friend adub exp (const badouble&); 
  friend adub log (const badouble&); 
  friend adub sqrt (const badouble&);
  friend adub sin (const badouble&); 
  friend adub cos (const badouble&);
  friend adub tan (const badouble&);
  friend adub asin (const badouble&);
  friend adub acos (const badouble&);
  friend adub atan (const badouble&); 
  friend adub atan2 (const badouble&,const badouble&); 
  friend adub pow (const badouble&,double);
  friend adub pow (const badouble&,const badouble&);
  friend adub log10 (const badouble&);

  /* User defined version of logarithm to test extend_quad macro */

  friend adouble myquad(const badouble&);

  /* Additional ANSI C standard Math functions Added by DWJ on 8/6/90 */

  friend adub sinh (const badouble&);
  friend adub cosh (const badouble&);
  friend adub tanh (const badouble&);
  friend adub ceil (const badouble&);
  friend adub floor (const badouble&);
  friend adub asinh (const badouble&);
  friend adub acosh (const badouble&);
  friend adub atanh (const badouble&);
  friend adub fabs (const badouble&);
  friend adub max (const badouble&, const badouble&);
  friend adub min (const badouble&, const badouble&);
  friend adub ldexp (const badouble&,int);
  friend adub frexp (const badouble&,int*);
  friend adub erf (const badouble&);
  /* End of ANSI C Additions */
};

/* The derived classes */
/* 
   The class Adub
   ---- Basically used as a temporary result.  The address for an
        adub is usually generated within an operation.  That address
        is "freed" when the adub goes out of scope (at destruction time).
   ---- operates just like a badouble, but it has a destructor defined for it.
*/

class adub:public badouble{
  friend class adouble;
 protected:
  adub(locint lo):badouble(lo){};
  adub():badouble(0){
      cout << "ADOL-C error: illegal default construction of adub variable\n" ;
      exit(-2);
         };
  adub(double):badouble(0){
      cout << "ADOL-C error: illegal  construction of adub variable from double\n" ;
      exit(-2);
         };
 public:
  friend double value(const badouble&);
  friend ostream& operator << (ostream&, const badouble&);
  friend istream& operator >> (istream&, const badouble&);

  friend int operator != (const badouble&,const badouble&);
  friend int operator != (double,const badouble&);
  friend int operator != (const badouble&,double);
  friend int operator == (const badouble&,const badouble&);
  friend int operator == (double,const badouble&);
  friend int operator == (const badouble&,double);
  friend int operator >= (const badouble&,const badouble&);
  friend int operator >= (double,const badouble&);
  friend int operator >= (const badouble&,double);
  friend int operator <= (const badouble&,const badouble&);
  friend int operator <= (double,const badouble&);
  friend int operator <= (const badouble&,double);
  friend int operator > (const badouble&,const badouble&);
  friend int operator > (double,const badouble&);
  friend int operator > (const badouble&,double);
  friend int operator < (const badouble&,const badouble&);
  friend int operator < (double,const badouble&);
  friend int operator < (const badouble&,double);

  /* End of Comparision Operators */

  friend adub operator + (const badouble& x); //{return x + 0.0 ;} ; 
  friend adub operator + (const badouble&,const badouble&); 
  friend adub operator + (double, const badouble&); 
  friend adub operator + (const badouble&, double); 
  friend adub operator - (const badouble& x ,double y); //{return (-y)+x;}; 
  friend adub operator - (const badouble&,const badouble&); 
  friend adub operator - (double, const badouble&); 
  friend adub operator - (const badouble& x); //{return 0.0 - x;}; 
  friend adub operator * (const badouble&,const badouble&); 
  friend adub operator * (double, const badouble& ); 
  friend adub operator * (const badouble& x, double y); //{return y*x;}; 
  friend adub operator / (const badouble& x, double y); // {return (1.0/y)*x;}; 
  friend adub operator / (const badouble&,const badouble&); 
  friend adub operator / (double,const badouble&); 
  friend adub exp (const badouble&); 
  friend adub log (const badouble&); 
  friend adub sqrt (const badouble&);
  friend adub sin (const badouble&); 
  friend adub cos (const badouble&);
  friend adub tan (const badouble&);
  friend adub asin (const badouble&);
  friend adub acos (const badouble&);
  friend adub atan (const badouble&); 
  friend adub atan2 (const badouble&,const badouble&); 
  friend adub pow (const badouble&,double);
  friend adub pow (const badouble&,const badouble&);
  friend adub log10 (const badouble&);

  /* User defined version of logarithm to test extend_quad macro */

  friend adouble myquad(const badouble&);

  /* Additional ANSI C standard Math functions Added by DWJ on 8/6/90 */

  friend adub sinh (const badouble&);
  friend adub cosh (const badouble&);
  friend adub tanh (const badouble&);
  friend adub ceil (const badouble&);
  friend adub floor (const badouble&);
  friend adub asinh (const badouble&);
  friend adub acosh (const badouble&);
  friend adub atanh (const badouble&);
  friend adub fabs (const badouble&);
  friend adub max (const badouble&, const badouble&);
  friend adub min (const badouble&, const badouble&);
  friend adub ldexp (const badouble&,int);
  friend adub frexp (const badouble&,int*);
  friend adub erf (const badouble&);

/*  removed 1/95 
  friend adub operator*(const badoublev &op1, const doublev &op2);
  friend adub operator*(const doublev &op1, const badoublev &op2);
*/
  friend adub operator*(const badoublev &op1, double* op2);
  friend adub operator*(double* op1, const badoublev &op2);
  friend adub operator*(const badoublev &op1, const badoublev &op2);
  /* excluded because g++ warnings
  friend adub operator*(const badoublev &op1, const doublev &op2);
  friend adub operator*(const doublev &op1, const badoublev &op2);
  */

#ifdef overwrite
  ~adub();
#endif
  };
/*
  The class adouble.
  ---Derived from badouble.  Contains the standard constructors/destructors.
  ---At construction, it is given a new address, and at destruction, that
     address is freed.
*/
class adouble:public badouble{
 public:
  adouble(const adub& a);
  adouble(const adouble&);
  adouble();
  adouble(double);
  adub operator++(int);
  adub operator--(int);
  badouble&  operator++();
  badouble&  operator--();

#ifdef overwrite
  ~adouble();
#endif
  adouble& operator = (double);
  adouble& operator = (const badouble&);
  adouble& operator = (const adouble&);
  adouble& operator = (const adub&);
};

inline adub operator + (const badouble& x){return x + 0.0 ;} 
inline adub operator - (const badouble& x ,double y){return (-y)+x;}
inline adub operator - (const badouble& x){return 0.0 - x;}
inline adub operator * (const badouble& x, double y){return y*x;}
inline adub operator / (const badouble& x, double y){return (1.0/y)*x;}


/* Supporting routines. */

void trace_on(short,int&); 
void trace_on(short);

void trace_off(int);
void trace_off();

/* Include the definition of the vector classes */

#ifdef conditional
class asub:public badouble
{
  locint base,offset;
 public:
  asub(locint start,locint index);
#ifdef overwrite
  ~asub();
#endif
  asub& operator <<= (double);
  asub& operator = (double x);
  asub& operator = (const adub&);
  asub& operator = (const badouble&);
};

class along:public adouble
{
 public:
  along(const adub& a);
  along(const along&);
  along(){};
  along(int);
  along& operator = (int);
  along& operator = (const badouble&);
  along& operator = (const along&);
  along& operator = (const adub&);
};
#endif

#include "avector.h"






