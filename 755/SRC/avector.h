/* 
  --------------------------------------------------------------------
  file avector.h of ADOL-C version 1.6 as of January 1,   1995 
  Included in
              adouble.h, and hence ---->
	                                 adouble.c
					 avector.c
					 All ADOL-C applications.
  --------------------------------------------------------------------

  Avector.h defines classes of vectors and matrices.

  badoublev  --> class of basic active vectors. 
  adubv      --> class of temporary active vectors.
                 (derived from badoublev.  Contains copy constructors, 
                  destructors.)
  adoublev   --> class of active vectors. (derived from badoublev,
                 contains the standard constructors and destructors.
		 
*/

/* ----- Forward declaration ----- */
class badoublev;
class adoublev;  
class adubv;
/* removed 1/95 
class doublev;
*/
class err_retu;
class asubv;


class err_retu
{
char* message;
public:
err_retu(char* x){printf("%s \n",x);};
};


/* ----- DECLARATION OF VECTOR CLASSES ----- */     
/* Passive vectors and matrices */


/* REMOVED 1/95  */


/* ----- End of passive section ----- */

class badoublev
{
 protected:
  locint start_loc;  /* Starting location of vector in store */
  int size;          /* Size of the vector */
  badoublev(){};
  badoublev(int lo, int sz){start_loc = lo; size=sz;};
  badoublev(const badoublev& a){start_loc = a.start_loc; size=a.size;};
  
 public:

  /* Access functions */
  int sz() const {return size;}  /* Get the size of the vector */
  int loc() const {return start_loc;}  /* Get the size of the vector */

#ifdef conditional
  asub operator[](const along&) const;
#endif


/* excluded before 1/95 
  badoublev& operator >>= (doublev& );
  badoublev& operator <<= (doublev& );
  badoublev& operator >>= (double* );
  badoublev& operator <<= (double* );
*/

  badouble operator[](int) const;  /* Can access component like an array */

  badoublev& operator+=(const badoublev&);
  badoublev& operator-=(const badoublev&);
  badoublev& operator*=(double);
  badoublev& operator/=(double);
/* removed 1/95
  badoublev& operator-=(const doublev&);
  badoublev& operator+=(const doublev&);
*/
  badoublev& operator-=(double*);
  badoublev& operator+=(double*);
  badoublev& operator*=(const badouble& );
  badoublev& operator/=(const badouble& );
  friend adubv operator/(const badoublev &op1, const badouble &n);
  inline friend adubv operator/(const badoublev &op1, double n);
/*  removed 1/95
  badoublev& operator= (const doublev&);
*/
  badoublev& operator= (const badoublev&);
  badoublev& operator= (const adubv &y);
  badoublev& operator= (const adoublev &y);

  friend ostream& operator << (ostream&, const badoublev&);

  friend adubv operator+ (const badoublev &x);
  friend adubv operator- (const badoublev &x);
  
  /* overload operations */
  friend adubv operator+(const badoublev &op1,const badoublev &op2);
  friend adubv operator-(const badoublev &op1,const badoublev &op2);
  friend adubv operator*(const badoublev &op1, double n);
  friend adubv operator*(double n, const badoublev &op1);
  friend adub operator*(const badoublev &op1, const badoublev &op2);

  /* overloaded for interaction of constant and active vectors */
/* removed 1/95
  friend adubv operator+(const badoublev &op1, const doublev &op2);
  friend adubv operator+(const doublev &op1, const badoublev &op2);
*/
  friend adubv operator+(const badoublev &op1, double* op2);
  friend adubv operator+(double* op2, const badoublev &op1);
/* removed 1/95
  friend adubv operator-(const badoublev &op1, const doublev &op2);
  friend adubv operator-(const doublev &op1, const badoublev &op2);
*/
  friend adubv operator-(const badoublev &op1, double* op2);
  friend adubv operator-(double* op1, const badoublev &op2);
/* removed 1/95
  friend adub operator*(const badoublev &op1, const doublev &op2);
  friend adub operator*(const doublev &op1, const badoublev &op2);
*/
  friend adub operator*(const badoublev &op1, double* op2);
  friend adub operator*(double* op1, const badoublev &op2);
  
  /* overloaded for interaction of active scalars and active vectors */
/* removed 1/95
  friend adubv operator/(const doublev &op1, const badouble &n);
*/
  friend adubv operator*(const badoublev &op1, const badouble &n);
  friend adubv operator*(const badouble &n, const badoublev &op1);
  /* excluded operations */
  err_retu operator>>=(double) {return("ADOL-C error: illegal argument combination for operator >>=\n"); };
  err_retu operator<<=(double) {return("ADOL-C error: illegal argument combination for operator <<=\n"); };
  err_retu operator+= (double) {return("ADOL-C error: illegal argument combination for operator +=\n"); };
  err_retu operator-= (double) {return("ADOL-C error: illegal argument combination for operator -=\n"); };
  inline friend err_retu operator+(const badoublev,double) {return("ADOL-C error: illegal argument combination for operator +\n"); };
  inline friend err_retu operator-(const badoublev,double) {return("ADOL-C error: illegal argument combination for operator -\n"); };
  inline friend err_retu operator+(double,const badoublev) {return("ADOL-C error: illegal argument combination for operator +\n"); };
  inline friend err_retu operator-(double,const badoublev) {return("ADOL-C error: illegal argument combination for operator -\n"); };
};

class adubv:public badoublev{
  adubv(int lo,int sz){start_loc=lo;size=sz;};
/* removed 1/95
  adubv(doublev&);
*/
  adubv():badoublev(0,0){
      cout << "ADOL-C error: illegal default construction of adub variable\n" ;
      exit(-2);
         };

 public:    
/* removed 1/95
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
/* removed 1/95
  friend adubv operator+(const badoublev &op1, const doublev &op2);
  friend adubv operator+(const doublev &op1, const badoublev &op2);
  friend adubv operator-(const badoublev &op1, const doublev &op2);
  friend adubv operator-(const doublev &op1, const badoublev &op2);
  friend adubv operator/(const doublev &op1, const badouble &n);
  friend adubv operator*(const doublev &op1, const badouble &n);
  friend adubv operator*(const badouble &n, const doublev &op1);
*/
  friend adubv operator/(const badoublev &op1, const badouble &n);
  inline friend adubv operator/(const badoublev &op1, double n);
  friend adubv operator+ (const badoublev &x);
  friend adubv operator- (const badoublev &x);
  friend adubv operator+(const badoublev &op1,const badoublev &op2);
  friend adubv operator-(const badoublev &op1,const badoublev &op2);
  friend adubv operator*(const badoublev &op1, double n);
  friend adubv operator*(double n, const badoublev &op1);
  /* excluded because g++ warnings
  friend adubv operator+(const badoublev &op1, const doublev &op2);
  friend adubv operator+(const doublev &op1, const badoublev &op2);
  */
  friend adubv operator+(const badoublev &op1, double* op2);
  friend adubv operator+(double* op2, const badoublev &op1);
  /* excluded because g++ warnings
  friend adubv operator-(const badoublev &op1, const doublev &op2);
  friend adubv operator-(const doublev &op1, const badoublev &op2);
  */
  friend adubv operator-(const badoublev &op1, double* op2);
  friend adubv operator-(double* op1, const badoublev &op2);
  /* excluded because g++ warnings
  friend adubv operator/(const doublev &op1, const badouble &n);
  */
  friend adubv operator*(const badoublev &op1, const badouble &n);
  friend adubv operator*(const badouble &n, const badoublev &op1);
#ifdef overwrite
  ~adubv();
#endif
};

class adoublev:public badoublev
{
  friend class adoublem;
  adoublev(){};
 public:
  adoublev(const adubv& a);
  adoublev(const adoublev&);
  adoublev(int sz);
//  adoublev(int n, double *values);
/* removed 1/95
  adoublev(doublev&);
*/
#ifdef overwrite
  ~adoublev();
#endif
/* removed 1/95
  adoublev& operator= (const doublev &y);
*/
  adoublev& operator= (const badoublev&);
  adoublev& operator= (const adoublev&);
  adoublev& operator= (const adubv&);
  adoublev& operator= (double y);
  adoublev& operator= (double* y);
/* removed 1/95
  adoublev& operator >>= (doublev& );
  adoublev& operator <<= (doublev& );
*/
  adoublev& operator >>= (double* );
  adoublev& operator <<= (double* );
};


inline adubv operator / (const badoublev& x, double y){return (1.0/y)*x;}


/* Active matrix class */
class adoublem
{
  int n, m;          /* Size of the matrix */
  adoublev *index;     /* So each row is an adoublev */
 public:
  adoublem(int n, int m);
  ~adoublem();
  adoublev& operator[](int i);  /* Can access component like an array */
#ifdef conditional
   asubv operator[](const along&);
#endif
};

#ifdef conditional
class asubv:public badoublev
{
  locint base,offset,begin;
  public:
  asubv(adoublev* start, locint index);
#ifdef overwrite
    ~asubv();
#endif
/* removed 1/95
  asubv& operator <<= (doublev&);
  asubv& operator = (doublev);
*/
  asubv& operator <<= (double*);
  asubv& operator = (double*);
  asubv& operator = (const adubv&);
  asubv& operator = (const badoublev&);
 };
#endif


