//****************************************************************************
//  r2d2lri.cpp
//  Implementation file for the DoubleIntegral class.
//
//  Code written by Ian Robinson and Michael Hill.
//  Last updated: 23 April 2001
//
//  Refer to the header file r2d2lri.h and to the files r2d2lri_sample.cpp and
//  README.txt for comments on, and examples of, the use of this class.
//
//****************************************************************************

#include "r2d2lri.h"

// CONSTRUCTORS

//****************************************************************************
DoubleIntegral::DoubleIntegral(double req_rel_acc, double req_abs_acc,
                               unsigned max_evals)

// Default values for all three of the parameters are provided in the prototype
// for this constructor.
{
	rel_tol = req_rel_acc;
   abs_tol = req_abs_acc;
   max_points = min(max_evals, DEFAULT_MAX_POINTS);
   xdir = NOT_CHECKED;
   ydir = NOT_CHECKED;
   max_iters = DEFAULT_C_SIZE - 1;
   c = new double [DEFAULT_C_SIZE];
   c[0] = 0.0;
   e = new double [DEFAULT_C_SIZE-7];
}

//****************************************************************************
DoubleIntegral::DoubleIntegral(double A, double B, BOUNDARY_FUNCTION G,
                               BOUNDARY_FUNCTION H, INTEGRAND_FUNCTION F,
                               double req_rel_acc, double req_abs_acc,
                               unsigned max_evals)

// The prototype for this constructor provides default values for the last
// three parameters, req_rel_acc, req_abs_acc and max_evals.
{
   a = A;  b = B;
   g = G;  h = H;
   f = F;
   rel_tol = req_rel_acc;
   abs_tol = req_abs_acc;
   max_points = min(max_evals,DEFAULT_MAX_POINTS);
   set_region_type();
   set_initial_values();
   max_iters = DEFAULT_C_SIZE - 1;
   c = new double [DEFAULT_C_SIZE];
   c[0] = 0.0;
   e = new double [DEFAULT_C_SIZE-7];
}

// COPY CONSTRUCTOR

//****************************************************************************
DoubleIntegral::DoubleIntegral(const DoubleIntegral& I)

// Copies the contents of all private data members of I into the corresponding
// private data members of the new object.
{
   unsigned i;

   a = I.a;
   b = I.b;
   g = I.g;
   h = I.h;
   f = I.f;
   rel_tol = I.rel_tol;
   abs_tol = I.abs_tol;
   max_points = I.max_points;
   xdir = I.xdir;
   ydir = I.ydir;
   x_transform = I.x_transform;
   y_transform = I.y_transform;
   cubature = I.cubature;
   err_est = I.err_est;
   fun_vals = I.fun_vals;
   ifail = I.ifail;
   max_iters = I.max_iters;
   c = new double [max_iters + 1];
   for (i = 0; i < max_iters; i++)
     c[i] = I.c[i];
   e = new double [max_iters - 6];
   for (i = 0; i < max_iters - 7; i++)
     e[i] = I.e[i];
}

// DESTRUCTOR

//****************************************************************************
DoubleIntegral::~DoubleIntegral()
{
   delete [] c;
   delete [] e;
}

// ASSIGNMENT OPERATOR

//****************************************************************************
void DoubleIntegral::operator =(const DoubleIntegral &I)
{
   unsigned i;
   double *temp;

   a = I.a;
   b = I.b;
   g = I.g;
   h = I.h;
   f = I.f;
   rel_tol = I.rel_tol;
   abs_tol = I.abs_tol;
   max_points = I.max_points;
   xdir = I.xdir;
   ydir = I.ydir;
   x_transform = I.x_transform;
   y_transform = I.y_transform;
   cubature = I.cubature;
   err_est = I.err_est;
   fun_vals = I.fun_vals;
   ifail = I.ifail;
   if (max_iters != I.max_iters)
   {
     temp = new double [I.max_iters + 1];
     delete [] c;
     c = temp;
     temp = new double [I.max_iters - 6];
     delete [] e;
     e = temp;
     max_iters = I.max_iters;
   }
   for (i = 0; i < max_iters; i++)
     c[i] = I.c[i];
   for (i = 0; i < max_iters - 7; i++)
     e[i] = I.e[i];
}

// MODIFICATION MEMBER FUNCTIONS

// The calls to the procedure set_initial_values() in these modification
// functions may appear to be redundant since this initializing procedure is
// called again when the new integral thus created is evaluated.  However,
// the values of the cubature, the error estimate and the number of function
// values have been re-initialized as soon as any aspect of the integral is
// changed so that any possible misinterpretation or misuse of the residual
// values of these variables is avoided.

//****************************************************************************
void DoubleIntegral::set_outer_interval(double A, double B)
{
   a = A;
   b = B;
   xdir = set_interval_type(infinite(a), infinite(b));
   x_transform = set_transform_type(xdir);
   set_initial_values();
}

//****************************************************************************
void DoubleIntegral::set_inner_interval(BOUNDARY_FUNCTION G,
                                        BOUNDARY_FUNCTION H)
{
   g = G;
   h = H;
   ydir = set_interval_type(infinite((*g)(a)), infinite((*h)(a)));
   y_transform = set_transform_type(ydir);
   set_initial_values();
}

//****************************************************************************
void DoubleIntegral::set_integrand(INTEGRAND_FUNCTION F)
{
   f = F;
   set_initial_values();
}

//****************************************************************************
void DoubleIntegral::set_rel_tol(double req_rel_acc)
{
   rel_tol = req_rel_acc;
   set_initial_values();
}

//****************************************************************************
void DoubleIntegral::set_abs_tol(double req_abs_acc)
{
   abs_tol = req_abs_acc;
   set_initial_values();
}

//****************************************************************************
void DoubleIntegral::set_max_evals(unsigned max_evals)
{
   max_points = max_evals;
   set_initial_values();
}

//****************************************************************************
void DoubleIntegral::set_new_integral(double A, double B, BOUNDARY_FUNCTION G,
                                      BOUNDARY_FUNCTION H, INTEGRAND_FUNCTION F)

// Note: rel_tol, abs_tol and max_points remain unchanged.
{
   a = A;  b = B;
   g = G;  h = H;
   f = F;
   set_region_type();
   set_initial_values();
}

// ACCESS FUNCTIONS

//****************************************************************************
void DoubleIntegral::get_outer_interval(double& A, double& B) const
{
   A = a;
   B = b;
}

//****************************************************************************
void DoubleIntegral::get_inner_interval(BOUNDARY_FUNCTION& G,
                                        BOUNDARY_FUNCTION& H) const
{
   G = g;
   H = h;
}

//****************************************************************************
void DoubleIntegral::get_integrand(INTEGRAND_FUNCTION& F) const
{
   F = f;
}

//****************************************************************************
double DoubleIntegral::get_rel_tol() const
{
   return rel_tol;
}

//****************************************************************************
double DoubleIntegral::get_abs_tol() const
{
   return abs_tol;
}

//****************************************************************************
unsigned DoubleIntegral::get_max_evals() const
{
   return max_points;
}

//****************************************************************************
double DoubleIntegral::value() const
{
   return cubature;
}

//****************************************************************************
double DoubleIntegral::rel_err_est() const
{
   return err_est;
}

//****************************************************************************
double DoubleIntegral::abs_err_est() const
{
   if (err_est == -1.0)
     return -1.0;
   return cubature == 0.0 ? err_est : err_est*fabs(cubature);
}

//****************************************************************************
int DoubleIntegral::evals() const
{
   return fun_vals;
}

//****************************************************************************
int DoubleIntegral::error_flag() const
{
   return ifail;
}

// INTEGRAL EVALUATION FUNCTIONS

//****************************************************************************
double DoubleIntegral::evaluate()
{
   integrate();
   return cubature;
}

//****************************************************************************
double DoubleIntegral::evaluate(double req_rel_acc)
{
   rel_tol = req_rel_acc;
   integrate();
   return cubature;
}

//****************************************************************************
double DoubleIntegral::evaluate(double& rel_err_est, int& flag, int& evals)
{
   integrate();
   rel_err_est = err_est;
   flag = ifail;
   evals = fun_vals;
   return cubature;
}

//****************************************************************************
double DoubleIntegral::evaluate(double A, double B, BOUNDARY_FUNCTION G,
                                BOUNDARY_FUNCTION H, INTEGRAND_FUNCTION F,
                                double& rel_err_est, int& flag,
                                double req_rel_acc, double req_abs_acc,
                                unsigned max_evals)

// The prototype for this function provides default values for the last three
// parameters, req_rel_acc, req_abs_acc and max_evals.
{
   a = A;  b = B;
   g = G;  h = H;  f = F;
   rel_tol = req_rel_acc;
   abs_tol = req_abs_acc;
   max_points = max_evals;
   set_region_type();
   integrate();
   flag = ifail;
   rel_err_est = err_est;
   return cubature;
}

// PRIVATE MEMBER FUNCTIONS

//****************************************************************************
void DoubleIntegral::set_region_type()

// Characterizes the integration region as the product of two one-dimensional
// intervals, each of which may be finite, semi-infinite or infinite, and sets
// the appropriate transformation flags.
{
   BOOL a_infinite = infinite(a);
   BOOL b_infinite = infinite(b);
   double r;                                 // A "random" value
   
   if (a_infinite)
      if (b_infinite)
         r = 2.149783;                       // Must be (-INFINITY, INFINITY)
      else
         r = b + 6.451372;                   // Must be [b, INFINITY)
   else
      if (b_infinite)
         r = a - 6.451372;                   // Must be (-INFINITY, a]
      else
         r = (a + b)/2.130684;               // Must be [a, b]
      
   xdir = set_interval_type(a_infinite, b_infinite);
   ydir = set_interval_type(infinite((*g)(r)), infinite((*h)(r)));
   x_transform = set_transform_type(xdir);
   y_transform = set_transform_type(ydir);
}

//****************************************************************************
INTERVAL DoubleIntegral::set_interval_type(BOOL lower_inf, BOOL upper_inf) const

// Returns a description of a one-dimensional integration interval:
// FINITE, UPPER_INFINITE ([a,infinity)), LOWER_INFINITE ((-infinity,b])
// or INFINITE((-infinity,infinity)).
{
   if (lower_inf)
      if (upper_inf)
         return INFINITE;
      else
         return LOWER_INFINITE;
   else
      if (upper_inf)
         return UPPER_INFINITE;
      else
         return FINITE;
}

//****************************************************************************
TRANSFORM DoubleIntegral::set_transform_type(INTERVAL direction) const

// Returns a value corresponding to the type of transformation to be applied
// before attempting to evaluate the integral.
{
   if (direction == NOT_CHECKED || direction == FINITE)
     return NONE;
   else
     return RATIONAL;
}

//****************************************************************************
void DoubleIntegral::set_initial_values()

// Initializes cubature, err_est, fun_vals and ifail.
{
   cubature = 0.0;
   err_est = -1.0;
   fun_vals = 0;
   ifail = -1;
}

//****************************************************************************
void DoubleIntegral::integrate()

// Calls the core integration routine to evaluate cubature, err_est, ifail and
// fun_vals.

// If the integration region is the plane, a radial transformation is first
// applied (subject to there being a likelihood of success with this trans-
// formation).  If a satisfactory result is not achieved, then the integral
// is recomputed using a rational transformation in each direction.  If a
// satisfactory result is still not achieved, then a logarithm transformation
// is applied in each direction.  The best of the results obtained using one,
// two or three of these transformations is accepted.

// If the integration region is semi-infinite, a rational transformation is
// first applied in the infinite direction(s).  If a satisfactory result is not
// achieved, then the integral is recomputed using a logarithm transformation in
// the infinite direction(s), with the better of the two results being accepted.

{
   int    fun_vals1 = 0;
   double cubature1, err_est1 = HUGE_VAL;
   BOOL   repeat;

   // First, check that the integration region has been characterized
   if (xdir == NOT_CHECKED || ydir == NOT_CHECKED)
   {
      cerr << "Integration region not fully specified.  Integration aborted."
           << endl;
      exit(1);
   }

   radial_transform = FALSE;

   // If the integration region is the plane, ...
   if (xdir == INFINITE && ydir == INFINITE)
   {
     // ... compute the optimal radius value for the radial transform
     r = optimal_radius(*f);

     if (r <= 0.0 || r > 6.3)
     {
       // If the radial transform is unlikely to produce a good
       // result, then try the rational transform instead
       fun_vals1 = 36;
     }
     else
     {
       // Try the radial transform
       radial_transform = TRUE;
       twoPIr2 = TWO_PI*r*r;
       core();
       fun_vals = 2*fun_vals + 36;
       if (ifail != 0)
       {
         //If it doesn't work, save values and try the rational transform
         fun_vals1 = fun_vals;
         cubature1 = cubature;
         err_est1 = err_est;
         radial_transform = FALSE;
       }
     }

     // If the radial transform was skipped or didn't work, ...
     if (!radial_transform)
     {
       // ... apply the rational transform
       core();
       fun_vals += fun_vals1;
       if (ifail != 0)
       {
         // If it doesn't work, save values and try the logarithm transform
         if (err_est < err_est1)
         {
           cubature1 = cubature;
           err_est1 = err_est;
         }
         fun_vals1 = fun_vals;
         x_transform = LOGARITHM;
         y_transform = LOGARITHM;

         // Logarithm transform for integration over the plane
         core();
         fun_vals += fun_vals1;
         if (ifail != 0)
         {
           // If it doesn't work, return the best result
           if (err_est > err_est1)
           {
             cubature = cubature1;
             err_est = err_est1;
           }
         }

         // Reset transform variables
         x_transform = RATIONAL;
         y_transform = RATIONAL;
       }
     }
   }
   else
   {
     // Straightforward evaluation for finite or semi-infinite region
     core();

     // If not successful ...
     if (ifail != 0)
     {
       // ... check if region is semi-infinite, ...
       repeat = FALSE;
       if (x_transform == RATIONAL)
       {
         x_transform = LOGARITHM;
         repeat = TRUE;
       }
       if (y_transform == RATIONAL)
       {
         y_transform = LOGARITHM;
         repeat = TRUE;
       }

       // ... and, if so, save values and try the logarithm transform
       if (repeat)
       {
         cubature1 = cubature;
         err_est1 = err_est;
         fun_vals1 = fun_vals;
         core();
         fun_vals += fun_vals1;

         // If it doesn't work, return the best result
         if (ifail != 0)
         {
           if (err_est > err_est1)
           {
             cubature = cubature1;
             err_est = err_est1;
           }
         }

         // Reset transform variables
         if (x_transform == LOGARITHM)
           x_transform = RATIONAL;
         if (y_transform == LOGARITHM)
           y_transform = RATIONAL;
       }
     }
   }
}

//****************************************************************************
double DoubleIntegral::optimal_radius(INTEGRAND_FUNCTION f)

// Returns an appropriate value for the radius of the circle into which to map
// the plane when using the radial transformation.  If no appropriate value
// can be determined, the value -1 is returned.
//
// The methodology is as implemented in the code associated with Cubpack++ 
// (R. Cools, D. Laurie and L. Pluym, "Algorithm 764: Cubpack++: A C++ Package
// for Automatic Two-Dimensional Cubature", ACM Trans on Math Soft, Vol 23,
// No. 1, March, 1997, pp. 1-15).
{

   int    i,k,              // Loop and array indices
          index_start = 0;  // Loop index

   double Ia,               // Used to hold successive annulus cubatures
          Iabs[7],    	    // For the absolute values of the annulus cubatures
          Iabs_sum = 0.0,   // Sum of the absolute cubatures for the annuli
          Iabs_sum2,        // Iabs_sum/2
          Iabs_temp = 0.0,  // Partial sums of Iabs[k]
          r;					 // Used in computation of the final optimal value

   // Calculating Iabs[k]

   for (k = 0; k < 7; k++)
   {
     Ia = 0.0;
     for (i = index_start; i < index_end[k]; i++)
     {
       Ia += (*f)(x_coord[i],y_coord[i]);
     };
     Ia *= mod_CH_weight[k];
     Iabs[k] = fabs(Ia);
     Iabs_sum += Iabs[k];
     index_start = index_end[k];
   };

   Iabs_sum2 = Iabs_sum/2.0;

   // Finding r: the first step

   k = -1;
   while (Iabs_temp <= Iabs_sum2)
   {
      k++;
      r = radius[k];
   	Iabs_temp += Iabs[k];
   };

   // The second step: finding the bit "left over"

   r = p[k] + (Iabs_temp - Iabs_sum2)/Iabs[k]*q[k];

   if (r > 0.0 && r < PI)
   {
     return sqrt(log(PI/r));
   }
   else
   {
     return -1.0;
   };
}

//****************************************************************************
void DoubleIntegral::core()

// Evaluates cubature, err_est, ifail and fun_vals.
{
   double l2, l4, s1, s2;

   BOOL   convergence = FALSE,
        limit_reached = FALSE,
             roundoff = FALSE;

   int      j = 1,                 // Index for the cubature array c
            n = 2,                 // Level of the displacement lattice
            m = SEED,              // Number of points in the current cubature
        two_m = 2*m,
       four_m = 4*m;

   // Compute initial cubature using the seed lattice

   set_initial_values();
   l4 = displaced_fibonacci(2,0,0);
   c[1] = l4/m;

   // Generate successive lattice rule approximations until convergence
   // is achieved or an abnormal condition is detected

   LATTICE lattice = L3;
   while (!convergence && !limit_reached && !roundoff)
   {
      switch (lattice)
      {
      case L3: // Generate L^3(2m) displaced lattice approximation
               s2 = generate_set(3,n);
               l2 = l4 + s2;
               c[++j] = l2/two_m;
               if (j > 2)
                  set_termination_flags(j,limit_reached,convergence,roundoff);
               lattice = L1;
               break;

      case L1: // Generate L^1(2m) displaced lattice approximation
               s1 = generate_set(1,n);
               l2 = l4 + s1;
               c[++j] = l2/two_m;
               set_termination_flags(j,limit_reached,convergence,roundoff);
               lattice = L2;
               break;

      case L2: // Generate L^2(2m) and L(4m) approximations
               l2 = l4 + generate_set(2,n);
               c[++j] = l2/two_m;
               l4 = l2 + s1 + s2;
               c[++j] = l4/four_m;
               set_termination_flags(j,limit_reached,convergence,roundoff);
               lattice = L3;

               // Update m and n
               m = four_m;
               two_m = 2*m;
               four_m = 4*m;
               n = 2*n;
       }
   }

   // Set cubature, err_est, ifail and fun_vals

   if ((j == 3) && approx_178)
      cubature = c[2];
   else
      cubature = c[j];
   err_est = max(SAFETY,err_est);
   if (convergence) ifail = 0;
   else if (roundoff) ifail = 1;
   else if (limit_reached) ifail = 2;
   switch (lattice)
   {
   case L1: fun_vals += two_m;
            break;
   case L2: fun_vals += 3*m;
            break;
   case L3: fun_vals += m;
   }
}

//***************************************************************************
double DoubleIntegral::generate_set(int set, int n)

// Allows generation of the x and y coordinates and corresponding weights for
// the displacement lattices.  n refers to the level of the Fibonacci lattice
// displacement process; set specifies which of the three displacement lattices
// L^1(2n), L^2(2n) or L^3(2n) is to be generated.
{
   int  p1, q1;
   double sum = 0.0;
   double s, t, u = 0.0;    // Used in the compensated Kahan summation 

   switch (set)
   {
      case 1: p1 = 0;       // L^1(2n): Even, Odd
              q1 = 1;
              break;

      case 2: p1 = 1;       // L^2(2n): Odd, Even
              q1 = 0;
              break;

      case 3: p1 = 1;       // L^3(2n): Odd, Odd
              q1 = 1;
   }

   for (int p = p1; p < n; p = p + 2)
      for (int q = q1; q < n; q = q + 2)
      {
         // Compensated Kahan summation
         s = displaced_fibonacci(n,p,q) - u;
         t = sum + s;
         u = (t - sum) - s;
         sum = t;
      }
      
   return sum;
}

//****************************************************************************
double DoubleIntegral::displaced_fibonacci(int n, int p, int q)

// Returns the weighted sum of function values to be used in the displaced
// Fibonacci rule. Only simple integer arithmetic is necessary to generate
// the references in the static arrays coordinate[] and weight[]. (Subsequent
// division by the number of points in the lattice yields the cubature
// approximation.)
{
   int denom = n*SEED,
       nprev = n*PREV_FIB_NUMBER,
       xnum = p*SEED - n,
       ynum = q*SEED - nprev;

   int    x_index, y_index;     // Indices of coordinate[] and weight[]
   double x, y, w, w1, w2;      // Values of coordinates and weights
   double sum = 0.0;            // Weighted sum of function values returned
   int    MULT = POINTS/denom;  // Used in the calculation of x_index & y_index 
   double rcx, rsx;             // Used when for the radial transformation
   double s, t, u = 0.0;        // Used in the compensated Kahan summations

   // Generate the displaced Fibonacci sum

   for (int counter = 0; counter < SEED; counter++)
   {
      // x and y are generated as for a rank 1 lattice rule

      xnum += n;
      x_index = (xnum % denom)*MULT;

      ynum += nprev;
      y_index = (ynum % denom)*MULT;

      w = weight[x_index] * weight[y_index];

      // Apply transformations

      if (w != 0)
      {
        y = coordinate[y_index];
        if (radial_transform)
        {
          if (y != 0.0)
          {
            rcx = r*cos2pix[x_index];
            rsx = r*sin2pix[x_index];

            // Compensated Kahan summation

            s = twoPIr2*w*(y*(*f)(y*rcx,y*rsx) + (*f)(rcx/y,rsx/y)/(y*y*y)) - u;
            t = sum + s;
            u = (t - sum) - s;
            sum = t;
          }
          else
            fun_vals--;
        }
        else
        {
          x = coordinate[x_index];
          map_to_ab(x,w1);
          if (w1 != 0)
          {
            map_to_gh(x,y,w2);
            if (w2 != 0)
            {
              // Compensated Kahan summation

              s = w*w1*w2*(*f)(x,y) - u;
              t = sum + s;
              u = (t - sum) - s;
              sum = t;
            }
            else
              fun_vals--;
          }
          else
            fun_vals--;
        }
      }
      else
        fun_vals--;
   }
   return sum;
}

//****************************************************************************
void DoubleIntegral::map_to_ab(double& x, double& w) const

// Transformation to [a,b] in the x direction.
// Precondition: x in [0,1].   Post-condition: x in [a,b].
{
   double z, z1, z2;     // Temporary variables

   switch (xdir)
   {
      case FINITE:          w = b - a;
                            x = a + w*x;
                            break;

      case UPPER_INFINITE:  if (x <= 0.0)
                               w = 0.0;
                            else
                            {
                               if (x_transform == RATIONAL)
                               {
                                  z = 1.0/x;
                                  w = z/x;
                                  x = a + z - 1.0;
                               }
                               else // x_transform == LOGARITHM
                               {
                                  if (x >= 1.0)
                                     w = 0.0;
                                  else
                                  {
                                     z = 1.0 - x;
                                     w = 1.0/z;
                                     x = a - log(z);
                                  }
                               }
                            }
                            break;

      case LOWER_INFINITE:  if (x <= 0.0)
                               w = 0.0;
                            else
                            {
                               if (x_transform == RATIONAL)
                               {
                                  z = 1.0/x;
                                  w = z/x;
                                  x = b - z + 1.0;
                               }
                               else // x_transform == LOGARITHM
                               {
                                 w = 1.0/x;
                                 x = b + log(x);
                               }
                            }
                            break;

      case INFINITE:        if (x <= 0.0 || x >= 1.0)
                               w = 0.0;
                            else
                            {
                               if (x_transform == RATIONAL)
                               {
                                  z1 = 1.0/(1.0 - x);
                                  z2 = 1.0/x;
                                  w = z1*z1 + z2*z2;
                                  x = z1 - z2;
                               }
                               else // x_transform == LOGARITHM
                               {
                                  z = 1.0/(1.0 - x);
                                  w = z/x;
                                  x = log(x*z);
                               }
                            }
   }
}

//****************************************************************************
void DoubleIntegral::map_to_gh(double x, double& y, double& w) const

// Transformation to [g(x),h(x)] in the y direction.
// Precondition: x in [a,b], y in [0,1].
// Post-condition: x in [a,b], y in [g(x),h(x)].
{
   double z, z1, z2;     // Temporary variables

   switch (ydir)
   {
      case FINITE:          z = (*g)(x);
                            w = (*h)(x) - z;
                            y = z + w*y;
                            break;

      case UPPER_INFINITE:  if (y <= 0.0)
                               w = 0.0;
                            else
                            {
                               if (y_transform == RATIONAL)
                               {
                                  z = 1.0/y;
                                  w = z/y;
                                  y = (*g)(x) + z - 1.0;
                               }
                               else // y_transform == LOGARITHM
                               {
                                  if (y >= 1.0)
                                     w = 0.0;
                                  else
                                  {
                                     z = 1.0 - y;
                                     w = 1.0/z;
                                     y = (*g)(x) - log(z);
                                  }
                               }
                            }
                            break;

      case LOWER_INFINITE:  if (y <= 0.0)
                               w = 0.0;
                            else
                            {
                               if (y_transform == RATIONAL)
                               {
                                  z = 1.0/y;
                                  w = z/y;
                                  y = (*h)(x) - z + 1.0;
                               }
                               else // y_transform == LOGARITHM
                               {
                                  w = 1.0/y;
                                  y = (*h)(x) + log(y);
                               }
                            }
                            break;

      case INFINITE:        if (y <= 0.0 || y >= 1.0)
                               w = 0.0;
                            else
                            {
                               if (y_transform == RATIONAL)
                               {                            
                                  z1 = 1.0/(1.0 - y);
                                  z2 = 1.0/y;
                                  w = z1*z1 + z2*z2;
                                  y = z1 - z2;
                               }
                               else // y_transform == LOGARITHM
                               {
                                  z = 1.0/(1.0 - y);
                                  w = z/y;
                                  y = log(y*z);
                               }
                            }
   }
}

//****************************************************************************
void DoubleIntegral::set_termination_flags(int j, BOOL& limit,
                                           BOOL& cnvrgnce, BOOL& round)

// Sets err_est to the error estimate for the current approximation and
// returns values for limit, cnvrgnce and round according to, respectively,
// whether the maximum number of points has been reached, the current
// approximation satisfies the accuracy requirement or rounding error has
// been detected in the current approximation.
{
   err_est = error_estimate(j);
   limit = j == max_iters;
   cnvrgnce = err_est <= max(rel_tol, abs_tol*fabs(c[j]));
   if (j > 12)
      round = diff(j-4,j-8) < 20.0*diff(j,j-4);
}

//****************************************************************************
double DoubleIntegral::error_estimate(int j)

// Returns an estimate of the error in c[j] (the current approximation).
{
   double d1,d2;

   switch(j)
   {
      case 3:  // 267 points
               d1 = fabs(c[3] - c[1]);
               d2 = fabs(c[2] - c[1]);
               if (d1 > d2)
               {
                  approx_178 = FALSE;
                  return c[3]==0.0 ? 4.0*d1 : 4.0*d1/fabs(c[3]);
               }
               else
               {
                  approx_178 = TRUE;
                  return c[2]==0.0 ? 4.0*d2 : 4.0*d2/fabs(c[2]);
               }

      case 5:  // 356 points
               e[1] = diff(5,1);
               conservative = (min_diff(5)!=0.0) && (e[1]/min_diff(5) < 60.0);
               if (conservative)
                  e[4] = max(e[1],max_diff(5));
               else
                  e[4] = go_fer(e[1]);
               return e[4];

      case 6:  // 712 points
               e[2] = diff(6,2);
               if (conservative || (e[4] > e[2]/10.0))
                  return max(e[2],diff(6,5));
               else
                  return go_fer(e[2]);

      case 7:  // 1068 points
               e[3] = diff(7,3);
               if (conservative || (e[4] > e[3]/10.0))
                  return max(e[3],diff(7,5),diff(7,6));
               else
                  return go_fer(e[3]);

      case 9:  // 1424 points
               e[4] = diff(9,5);
               if (((e[4] <= e[1]*e[1]) && (e[4] <= 0.5E-3)) || (e[4] <= 1E-9))
                  e[7] = pow(e[4],1.4);
               else
                  e[7] = max(e[4], max_diff(9));
               return e[7];

     case 10:  // 2848 points
               e[2] = diff(6,2);
               e[5] = diff(10,6);
               if ((((e[5] <= e[2]*e[2]) && (e[5] <= 0.5E-3))
                   || (e[5] <= 1E-9)) && (e[7] <= e[5]/10.0))
                  return max(pow(e[5],1.4), diff(10,9));
               else
                  return max(e[5], diff(10,9));

     case 11:  // 4272 points
               e[3] = diff(7,3);
               e[6] = diff(11,7);
               if ((((e[6] <= e[3]*e[3]) && (e[6] <= 0.5E-3))
                   || (e[6] <= 1E-9)) && (e[7] <= e[6]/10.0))
                  return max(pow(e[6],1.4), diff(11,9));
               else
                  return max(e[6], diff(11,9), diff(11,10));

     case 13:  // 5696 points
               e[7] = diff(13,9);
               if ((e[7] <= pow(e[4],1.5)) && (e[7] < 0.5E-6))
                  return pow(e[7],1.4);
               else
                  return max(e[7], max_diff(13));

     case 14:  // 11392 points
               e[8] = diff(14,10);
               if ((e[8] <= pow(e[5],1.5)) && (e[8] < 0.5E-6))
                  return max(pow(e[8],1.4), diff(14,13));
               else
                  return 2.0*max(e[8], diff(14,13));

     case 15:  // 17088 points
               e[9] = diff(15,11);
               if ((e[9] <= pow(e[6],1.5)) && (e[9] < 0.5E-6))
                  return max(pow(e[9],1.4), diff(15,13));
               else
                  return max(e[9], diff(15,14), diff(15,13));

     case 17:  // 22784 points
               e[10] = diff(17,13);
               return max(e[10], max_diff(17));

     default:  return diff(j,j-4);
   }
}

//****************************************************************************
double DoubleIntegral::diff(int j1, int j2) const

// Returns the relative difference between c[j1] and c[j2].
// Precondition: j1 > j2 > 0.
{
   return c[j1]==0.0 ? fabs(c[j1] - c[j2]) : fabs((c[j1] - c[j2])/c[j1]);
}

//****************************************************************************
double DoubleIntegral::max_diff(int j) const

// Returns the maximum relative difference between c[j] and c[j-k], k = 1,2,3.
// Precondition: j > 3.
{
   double curr_max = fabs(c[j] - c[j-1]);
   for (int k = j-2; k > j-4; k--)
   {
      double d = fabs(c[j] - c[k]);
      if (d > curr_max)
         curr_max = d;
   }
   return c[j]==0.0 ? curr_max : curr_max/fabs(c[j]);
}

//****************************************************************************
double DoubleIntegral::min_diff(int j) const

// Returns the minimum relative difference between c[j] and c[j-k], k = 1,2,3.
// Precondition: j > 3.
{
   double curr_min = fabs(c[j] - c[j-1]);
   for (int k = j-2; k > j-4; k--)
   {
      double d = fabs(c[j] - c[k]);
      if (d < curr_min)
         curr_min = d;
   }
   return c[j]==0.0 ? curr_min : curr_min/fabs(c[j]);
}

//****************************************************************************
double DoubleIntegral::go_fer(double x) const

// Let s be the maximum integer s.t. x <= 0.5*10^(-s).  If x is interpreted as
// the relative error in a number, then that number is correct to at least s
// significant figures.  go_fer returns x^1.4 if s > 3 and x^1.8, otherwise.
// Precondition: 0 < x < 1.
{
   if (x == 0.0)
      return 0.0;
   else
   {
      int s = (int)(NEG_LOG2 - log10(x));
      if (s > 3)
         return pow(x,1.4);
      else
         return pow(x,1.8);
   }
}

// FRIEND FUNCTIONS

//****************************************************************************
BOOL infinite(double x)

// Returns TRUE if x is +infinity or -infinity, and FALSE otherwise.
{
   return (fabs(x) >= INFINITY);
}

//****************************************************************************
double max(double x, double y)

// Returns maximum of x and y.
{
   return x > y ? x : y;
}

//****************************************************************************
double max(double x, double y, double z)

// Returns maximum of x, y and z.
{
   if (x > y)
      return x > z ? x : z;
   else
      return y > z ? y : z;
}

//****************************************************************************
double max(int x, int y)

// Returns maximum of x and y.
{
   return x > y ? x : y;
}

//****************************************************************************
unsigned min(unsigned x, unsigned y)

// Returns minimum of x and y.
{
   return x < y ? x : y;
}

