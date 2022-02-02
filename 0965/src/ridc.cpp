#include "ridc.h"

#include <stdio.h>

using namespace std;

///////////////////////////////////////////////////////
// Function Name: ridc_fe
// Usage: main explicit ridc loop 
//
///////////////////////////////////////////////////////
// Assumptions: none
//
///////////////////////////////////////////////////////
// Inputs:
//	   int order, order of integrator
//	   ode, ode definition + problem parameters
//         sol, initial condition
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//	    sol, solution at finest level, final time
//
///////////////////////////////////////////////////////
// Functions Called: corr_fe, integration matrices,
//    
//
///////////////////////////////////////////////////////
// Called By: none
//
///////////////////////////////////////////////////////
void ridc_fe(ODE * ode, int order, double *sol) {

  // unpack param
  double dt = ode->dt;
  double ti = ode->ti;
  int neq = ode->neq;
  int Nt = ode->nt;

  // setup memory stencil for RIDC
  int * stencilSize = new int[order];
  for (int p=0; p<(order-1); p++) {
    stencilSize[p] = p+2;
  }
  stencilSize[order-1] = 1;

  // memory footprint for quadrature approximations, for all levels
  double *** u = new double **[order];
  double *** f = new double **[order];
  for (int p=0; p<order; p++) {
    u[p] = new double*[stencilSize[p]];
    f[p] = new double*[stencilSize[p]];
    for (int i=0; i<stencilSize[p]; i++) {
      u[p][i] = new double[neq];
      f[p][i] = new double[neq];
    }
  }

  // this is a temporary variable for storing the solution at the new
  // time step (for each level)
  double **unew = new double*[order];
  double **fnew = new double*[order];
  for (int j=0;j<order;j++) {
    unew[j] = new double[neq];
    fnew[j] = new double[neq];
  }

  // setting variables for integration matrices
  double ***S = new double**[order];
  for (int i=0;i<order;i++) {
    S[i] = new double*[i+1];
    for (int j=0; j<i+1;j++) {
      S[i][j] = new double[i+1];
    }
  }

  // initialize integration matrices
  for (int j=1; j<order; j++) {
    integration_matrices(j+1,S[j]);
  }

  // get initial condition
  for (int i = 0; i<neq;i++) {
    u[0][0][i] = sol[i];
  }

  // compute f[0][0]
  ode->rhs(ti,u[0][0],f[0][0]);
  
  // copy initial condition for each RIDC level
  for (int j=1;j<order;j++){
    for (int i=0; i<neq; i++) {
      u[j][0][i]=u[0][0][i];
      f[j][0][i]=f[0][0][i];
    }
  }

  // variable containing flag for whether to shift memory in stencil
  bool *doMemShift = new bool[order];

  // flag for whether the node is able to compute the next step for
  // each correction level
  bool *doComputeStep = new bool[order];

  // counter for where the node is for each correction level
  int *ii = new int[order];

  for (int j=0; j<order; j++) {
    ii[j] = 0;
  }


  ///////////////////////////////////////////////////////////////////
  // START-UP ROUTINES
  ///////////////////////////////////////////////////////////////////

  // after startnum, all the nodes should be computing
  int startnum = (order*(order+1)/2-1);

  // special case
  if (order==1) {
    startnum = 1;
  }

  int **filter = new int*[order];
  for (int p=0; p<order; p++) {
    filter[p] = new int[startnum];
  }
  // set all of them to zero
  for (int p = 0; p<order; p++) {
    for (int q = 0; q<startnum; q++) {
      filter[p][q] = 0;
    }
  }

  if (order==1) {
    filter[0][0] = 1;
  } else {

    int q = 0;

    for (int p = 1; p<order; p++) {

      // all except level p step
      for (int pp = 0; pp<p; pp++) {
	filter[pp][q] = 1;
      }

      q++;

      // one step
      for (int count=0; count<p-1; count++) {
	filter[p][q] = 1;
	q++;
      }

      // all step, including level p
      for (int pp = 0; pp<=p; pp++) {
	filter[pp][q] = 1;
      }
      q++;
    } // loop over order
  } // create filter


  ///////////////////////////////////////////////////////////////////
  // MAIN TIME LOOP
  ///////////////////////////////////////////////////////////////////

  int i = 0;	// counter
  while(1) {
    {
      ///////////////////////////////////////////////////////////////////
      // PIPELINE UPDATE
      ///////////////////////////////////////////////////////////////////

#pragma omp parallel for
      for (int p=0; p<order; p++) {

	doComputeStep[p] = false;
	doMemShift[p] = false;

	if (  (ii[p] < Nt)  &&						\
	      ( (i >= (startnum-1)) || (filter[p][i] == 1) )		\
	      )
	{
	  doComputeStep[p] = true;
	}

	if (doComputeStep[p]) {
	  // Now do a step of the pth Corrector/Predictor
	  int index_into_array = min(ii[p], stencilSize[p]-1);
	  // compute current time
	  double t = ti + (ii[p])*dt;
	  if (p==0) {

	    // prediction step
	    ode->step(t,u[0][index_into_array],unew[0]);

	    // populate memory stencil for f (needed for order > 1)
	    ode->rhs(t+dt,unew[0],fnew[0]);

	  } else {
	    int index_sth = min(ii[p], p-1);

	    //correction step!
	    corr_fe(ode,u[p][index_into_array], f[p-1], S[p],
                   index_sth, p+1,t, unew[p]);

	    // populate memory stencil if not the last correction loop
	    if (p < (order-1) ) {
	      ode->rhs(t+dt,unew[p],fnew[p]);
	    }
	  }
	  ii[p]++;

	  if (ii[p] >= stencilSize[p]) {
	    doMemShift[p] = true;
	  }
	}  // end if doComputeStep[p]
      } // end p
#pragma omp barrier // ** BARRIER 1 **
    } // end pragma for loop
    
    
    // memshifts
    for (int p=0; p<order; p++) {
      if (doComputeStep[p]) {
	if (!(doMemShift[p])) {
	  int index_into_array = ii[p];
	  double *temp;

	  temp = u[p][index_into_array];
	  u[p][index_into_array] = unew[p];
	  unew[p] = temp;

	  temp = f[p][index_into_array];
	  f[p][index_into_array] = fnew[p];
	  fnew[p] = temp;
	} else {
	  double * temp;
	  temp = u[p][0];
	  for (int k=0; k<(stencilSize[p]-1); k++) {
	    u[p][k] = u[p][k+1];
	  }
	  u[p][stencilSize[p]-1] = unew[p];
	  unew[p] = temp;

	  temp = f[p][0];
	  for (int k=0; k<(stencilSize[p]-1); k++) {
	    f[p][k] = f[p][k+1];
	  }
	  f[p][stencilSize[p]-1] = fnew[p];
	  fnew[p] = temp;
	} // end doMemShift
      } // end doComputeStep
    } // end p

    if (ii[order-1] == Nt) {
      break;
    }
    i++;

  } // end while loop

  ///////////////////////////////////////////////////////////////////
  // END: MAIN TIME LOOP
  ///////////////////////////////////////////////////////////////////


  int p = order-1;
  for (int i =0; i<neq; i++) {
    sol[i] = u[p][stencilSize[p]-1][i];
  }

  // clear memory
  delete [] doMemShift;
  delete [] doComputeStep;
  delete [] ii;

  for (int p=0; p<order; p++) {
    delete [] unew[p];
    delete [] fnew[p];
  }
  delete [] unew;
  delete [] fnew;

  for (int p=0; p<order; p++)
    delete [] filter[p];
  delete [] filter;

  for (int p=0; p<order; p++) {
    for (int j =0; j<stencilSize[p]; j++) {
      delete [] u[p][j];
      delete [] f[p][j];
    }
    delete [] u[p];
    delete [] f[p];
  }

  delete [] u;
  delete [] f;

  for (int p=0; p<order; p++) {
    for (int j =0; j<p+1; j++) {
      delete [] S[p][j];
    }
    delete [] S[p];
  }
  delete [] S;

  delete [] stencilSize;
} // end ridc_fe


///////////////////////////////////////////////////////
// Function Name: ridc_be
// Usage: main ridc loop
//
///////////////////////////////////////////////////////
// Assumptions: none
//
///////////////////////////////////////////////////////
// Inputs:
//	   int order, order of integrator
//	   ode, ode definition + problem parameters
//         sol, initial condition
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//	    sol, solution at finest level, final time
//
///////////////////////////////////////////////////////
// Functions Called: corr_be, integration matrices,
//           initial condition
//
///////////////////////////////////////////////////////
// Called By: none
//
///////////////////////////////////////////////////////
void ridc_be(ODE * ode, int order, double *sol) {

  // unpack param
  double dt = ode->dt;
  double ti = ode->ti;
  int neq = ode->neq;
  int Nt = ode->nt;

  // setup memory stencil for RIDC
  int * stencilSize = new int[order];
  for (int p=0; p<(order-1); p++) {
    stencilSize[p] = p+2;
  }
  stencilSize[order-1] = 1;

  // memory footprint for quadrature approximations, for all levels
  double *** u = new double **[order];
  double *** f = new double **[order];
  for (int p=0; p<order; p++) {
    u[p] = new double*[stencilSize[p]];
    f[p] = new double*[stencilSize[p]];
    for (int i=0; i<stencilSize[p]; i++) {
      u[p][i] = new double[neq];
      f[p][i] = new double[neq];
    }
  }

  // this is a temporary variable for storing the solution at the new
  // time step (for each level)
  double **unew = new double*[order];
  double **fnew = new double*[order];
  for (int j=0;j<order;j++) {
    unew[j] = new double[neq];
    fnew[j] = new double[neq];
  }

  double ***S = new double**[order];
  for (int i=0;i<order;i++) {
    S[i] = new double*[i+1];
    for (int j=0; j<i+1;j++) {
      S[i][j] = new double[i+1];
    }
  }

  // initialize integration matrices
  for (int j=1; j<order; j++) {
    integration_matrices(j+1,S[j]);
  }

  // get initial condition
  for (int i = 0; i<neq;i++) {
    u[0][0][i] = sol[i];
  }
  
  // compute f[0][0]
  ode->rhs(ti,u[0][0],f[0][0]);
  
  // copy initial condition for each RIDC level
  for (int j=1;j<order;j++){
    // copy initial condition for each RIDC level
    for (int j=1;j<order;j++){
      for (int i = 0; i<neq;i++) {
	u[j][0][i]=u[0][0][i];
	f[j][0][i]=f[0][0][i];
      }
    }
  }

  
  // flag for whether to shift memory in stencil
  bool *doMemShift = new bool[order];

  // flag for whether the node is able to compute the next step for
  // each correction level
  bool *doComputeStep = new bool[order];

  // counter for where the node is for each correction level
  int *ii = new int[order];

  for (int j=0; j<order; j++) {
    ii[j] = 0;
  }

  ///////////////////////////////////////////////////////////////////
  // START-UP ROUTINES
  ///////////////////////////////////////////////////////////////////

  // after startnum, all the nodes should be computing
  int startnum = (order*(order+1)/2-1);

  // special case
  if (order==1) {
    startnum = 1;
  }

  int **filter = new int*[order];
  for (int p=0; p<order; p++) {
    filter[p] = new int[startnum];
  }
  // set all of them to zero
  for (int p = 0; p<order; p++) {
    for (int q = 0; q<startnum; q++) {
      filter[p][q] = 0;
    }
  }

  if (order==1) {
    filter[0][1] = 1;
  } else {

    int q = 0;

    for (int p = 1; p<order; p++) {

      // all except level p step
      for (int pp = 0; pp<p; pp++) {
	filter[pp][q] = 1;
      }

      q++;

      // one step
      for (int count=0; count<p-1; count++) {
	filter[p][q] = 1;
	q++;
      }

      // all step, including level p
      for (int pp = 0; pp<=p; pp++) {
	filter[pp][q] = 1;
      }
      q++;
    } // loop over order
  } // create filter


  ///////////////////////////////////////////////////////////////////
  // MAIN TIME LOOP
  ///////////////////////////////////////////////////////////////////

  int i = 0;	// counter
  while(1) {
    {
      ///////////////////////////////////////////////////////////////////
      // PIPELINE UPDATE
      ///////////////////////////////////////////////////////////////////

#pragma omp parallel for
      for (int p=0; p<order; p++) {

	doComputeStep[p] = false;
	doMemShift[p] = false;

	if (  (ii[p] < Nt)  &&						\
	      ( (i >= (startnum-1)) || (filter[p][i] == 1) )		\
	      )
	{
	  doComputeStep[p] = true;
	}

	if (doComputeStep[p]) {
	  // Now do a step of the pth Corrector/Predictor
	  int index_into_array = min(ii[p], stencilSize[p]-1);
	  // compute current time
	  double t = ti + (ii[p])*dt;
	  if (p==0) {

  
	    // prediction step
	    ode->step(t,u[0][index_into_array],unew[0]);
	    
	    // populate memory stencil for f
	    ode->rhs(t+dt,unew[0],fnew[0]);
	  } else {
	    int index_sth = min(ii[p], p-1);

	    //correction step!
	    corr_be(ode,u[p][index_into_array], f[p-1], S[p],
                    index_sth, p+1,t,unew[p]);

	    // populate memory stencil if not the last correction loop
	    if (p < (order-1) ) {
	      ode->rhs(t+dt,unew[p],fnew[p]);
	    }
	  }
	  ii[p]++;

	  if (ii[p] >= stencilSize[p]) {
	    doMemShift[p] = true;
	  }
	}  // end if doComputeStep[p]
      } // end p
#pragma omp barrier // ** BARRIER 1 **
    } // end pragma for loop

    
    // memshifts
    for (int p=0; p<order; p++) {
      if (doComputeStep[p]) {
	if (!(doMemShift[p])) {
	  int index_into_array = ii[p];
	  double *temp;

	  temp = u[p][index_into_array];
	  u[p][index_into_array] = unew[p];
	  unew[p] = temp;

	  temp = f[p][index_into_array];
	  f[p][index_into_array] = fnew[p];
	  fnew[p] = temp;
	} else {
	  double * temp;
	  temp = u[p][0];
	  for (int k=0; k<(stencilSize[p]-1); k++) {
	    u[p][k] = u[p][k+1];
	  }
	  u[p][stencilSize[p]-1] = unew[p];
	  unew[p] = temp;

	  temp = f[p][0];
	  for (int k=0; k<(stencilSize[p]-1); k++) {
	    f[p][k] = f[p][k+1];
	  }
	  f[p][stencilSize[p]-1] = fnew[p];
	  fnew[p] = temp;
	} // end doMemShift
      } // end doComputeStep
    } // end p

    if (ii[order-1] == Nt) {
      break;
    }
    i++;
  } // end while loop

  ///////////////////////////////////////////////////////////////////
  // END: MAIN TIME LOOP
  ///////////////////////////////////////////////////////////////////

  int p = order-1;
  for (int i =0; i<neq; i++) {
    sol[i] = u[p][stencilSize[p]-1][i];
  }

  // clear memory
  delete [] doMemShift;
  delete [] doComputeStep;
  delete [] ii;

  for (int p=0; p<order; p++) {
    delete [] unew[p];
    delete [] fnew[p];
  }
  delete [] unew;
  delete [] fnew;

  for (int p=0; p<order; p++)
    delete [] filter[p];
  delete [] filter;

  for (int p=0; p<order; p++) {
    for (int j =0; j<stencilSize[p]; j++) {
      delete [] u[p][j];
      delete [] f[p][j];
    }
    delete [] u[p];
    delete [] f[p];
  }

  delete [] u;
  delete [] f;

  for (int p=0; p<order; p++) {
    for (int j =0; j<p+1; j++) {
      delete [] S[p][j];
    }
    delete [] S[p];
  }
  delete [] S;

  delete [] stencilSize;
} // end ridc_be


///////////////////////////////////////////////////////
// Function lagrange_coeff
// Usage: generates coefficients for Lagrange interpolatory
//	  polynomials, L_{n,i}(x).  Used to generate weight w_i
//
///////////////////////////////////////////////////////
// Assumptions: nodes x are contained within [0,1], distinct
//
///////////////////////////////////////////////////////
// Inputs:
//	   x, (quadrature nodes)
//	   i, the i^{th} Lagrange polynomial
//	   Nx, number of quadrature nodes
//
///////////////////////////////////////////////////////
// Outputs: (By reference)
//	    L, coefficients for the Lagrange intepolatory
//	    polynomial.	 L is a vector of elements such that
//	    p(x) = L(0) + L(1)*x + L(2)*x^2 + ...
//
///////////////////////////////////////////////////////
// Functions Called: none
//
///////////////////////////////////////////////////////
// Called By: integration_matrix
//
///////////////////////////////////////////////////////
void lagrange_coeff(double *x, int Nx, int i, double *L) {

  int k=0; // counter
  int j; // counter
  int l; // counter
  double denom=1; // denominator for the coefficients

  // allocate memory, temporary storage variables
  double * c = new double[Nx];

  // set p(x) = 1
  c[0] = 1;
  for (j=1; j<Nx ;j++) {
    c[j] = 0.0;
  }

  for (j=0; j<Nx; j++) {
    if (j != i) {
      // set p(x) = (x-x[j])*p(x); store in L
      L[0] = -x[j]*c[0];
      for (l=1; l<=k ; l++) {
	L[l] = c[l-1] - x[j]*c[l];
      }
      k++;
      L[k] = 1;

      // update denominator
      denom = denom*(x[i]-x[j]);

      // copy from L back to c
      for (l=0; l<=k ; l++) {
	c[l] = L[l];
      }
    }
  }

  // scale coefficients correctly with denominator
  for (l=0; l < Nx; l++) {
    L[l] = L[l]/denom;
  }

  delete [] c;
} // end lagrange_coeff


///////////////////////////////////////////////////////
// Function Name: get_quad_weight
// Usage: generates quadrature weight, \int(L_{n,i}(x),x=a..b)
//
///////////////////////////////////////////////////////
// Assumptions: none
//
///////////////////////////////////////////////////////
// Inputs:
//	   a,b: range of integration (0->y)
//	   Nx:	number of quadrature nodes
//	   L[i]: coefficients for Lagrange poly,
//	     L[0] + L[1]x + L[2]x^2 + ...
//
///////////////////////////////////////////////////////
// Outputs:
//	    quadrature weight,
//	     w[i] = int(L_{n,i}(x),t=a..b)
//
///////////////////////////////////////////////////////
// Functions Called: none
//
///////////////////////////////////////////////////////
// Called By: integration_matrix
//
///////////////////////////////////////////////////////
double get_quad_weight(double *L, int Nx, double a, double b) {

  int j; // counter

  double answer=0;
  double an=1;
  double bn=1;

  for (j=0; j<Nx ; j++ ) {
    bn *= b;
    an *= a;
    answer += L[j]*(bn-an)/(j+1);
  }

  return answer;
} // end get_quad_weight


///////////////////////////////////////////////////////
// Function Name: integration_matrix
// Usage: generates integration matrices
//
///////////////////////////////////////////////////////
// Assumptions: none
//
///////////////////////////////////////////////////////
// Inputs:
//	   j, number of quadrature nodes
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//	     S = int(L_{n,i}(x),t=a..y[j])
//	     stored as S[j,i]
//
///////////////////////////////////////////////////////
// Functions Called: lagrange_coeff, get_quad_weight
//		     init_unif_nodes
//
///////////////////////////////////////////////////////
// Called By: ridc_fe, ridc_be
//
///////////////////////////////////////////////////////
void integration_matrices(int N, double **S) {

  int i; // counter
  int j; // counter

  double * x = new double [N];

  // allocate memory, temporary storage variables
  double * L = new double [N];

  init_unif_nodes(x,N,0,1);

  for (i=0; i<N; i++) {
    lagrange_coeff(x,N,i,L);
    for (j=0; j<N; j++) {
      S[j][i] = get_quad_weight(L,N,0,x[j]);
    }
  }

  delete [] x;
  delete [] L;
} // end integration_matrices


///////////////////////////////////////////////////////
// Function Name: init_unif_nodes
// Usage: initialize uniformly spaced quadrature nodes
//
///////////////////////////////////////////////////////
// Assumptions: creates uniformly spaced nodes
//
///////////////////////////////////////////////////////
// Inputs: Nx,	number of quadrature nodes
//	   a,b, domain where nodes are located
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//	    location of nodes, x
//
///////////////////////////////////////////////////////
// Functions Called: none
//
///////////////////////////////////////////////////////
// Called By: integration_matrices
//
///////////////////////////////////////////////////////
void init_unif_nodes(double *x, int Nx, double a, double b) {
  double h; // size of interval
  int j; //counter
  h = (b-a)/(Nx-1.0);
  for (j=0; j<Nx; j++){
    x[j] = a + j*h;
  }
} // end init_unif_nodes


///////////////////////////////////////////////////////
// Function Name: corr_fe
// Usage: computes one correction step
//
///////////////////////////////////////////////////////
// Assumptions: none
//
///////////////////////////////////////////////////////
// Inputs:
//	   uold, sol at prev step, current level
//	   fprev, sol' from previous level
//	   S, integration matrix
//	   index,
//	   level, which correction equation are we solving?
//	   t, time
//	   ode, ode definition + problem parameters
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//	    unew, new solution vector
//
///////////////////////////////////////////////////////
// Called By: ridc_fe
//
///////////////////////////////////////////////////////
void corr_fe(ODE * ode,
	     double * uold,
	     double ** fprev,
	     double ** S,
	     int index, int level,
	     double t,
	     double * unew) {

  ode->step(t,uold,unew);

  int neq = ode->neq;
  double dt = ode->dt;

  // forward euler update
  for (int j=0; j<neq; j++) {
    unew[j] -= dt*(fprev[index][j]);
    for (int i=0;i<level;i++) {
      unew[j] = unew[j] + dt*(S[index+1][i]-S[index][i])*fprev[i][j]*(level-1);
    }
  }
} // end corr_fe


///////////////////////////////////////////////////////
// Function Name: corr_be
// Usage: computes one correction step
//
///////////////////////////////////////////////////////
// Assumptions: none
//
///////////////////////////////////////////////////////
// Inputs:
//	   uold, sol at prev step, current level
//	   fprev, sol' from previous level
//	   S, integration matrix
//	   index,
//	   level, which correction equation are we solving?
//	   t, time
//	   ode, ode definition + problem parameters
//
///////////////////////////////////////////////////////
// Outputs: (by reference)
//	    unew, new solution vector
//
///////////////////////////////////////////////////////
// Called By: ridc_be
//
///////////////////////////////////////////////////////
void corr_be(ODE * ode,
	     double * uold,
	     double ** fprev,
	     double ** S,
	     int index, int level,
	     double t,
	     double * unew) {

  int neq = ode->neq;
  double dt = ode->dt;
  
  // backwards euler update
  for (int j=0; j<neq; j++) {
    uold[j] -= dt*(fprev[index+1][j]);
    for (int i=0;i<level;i++) {
      uold[j] += dt*(S[index+1][i]-S[index][i])*fprev[i][j]*(level-1);
    }
  }

  ode->step(t,uold,unew);
  

} // end corr_be
