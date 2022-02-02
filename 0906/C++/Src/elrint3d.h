#ifndef _ELRINT3D_H_
#define	_ELRINT3D_H_

/***************************************************************************************

   'elrint3d' is an automatic cubature routine designed to compute a numerical
   approximation to an integral in one of the following two forms:
 
		   b   d   f
   (1)    I   I   I   f(x,y,z) dz dy dx,
		   a   c   e

   where a, b, c, d, e and f are all constants (possibly INFINITY) or
 
		   b   h(x)   q(x,y)
   (2)    I   I      I       f(x,y,z) dz dy dx,
		   a   g(x)   p(x,y)
 
   where a and b are constants (again, possibly INFINITY) and g(x), h(x), p(x,y) 
   and q(x,y) are functions (whose values may be constant, including INFINITY). 
 
   The user must provide the integrand (either as a class or a function pointer), 
   the integration domain (as constants or function pointers) and the absolute
   and/or relative accuracy desired in the approximation.  The default maximum
   number of function evaluations permitted is 421,376, but a smaller maximum
   may be specified by the user.
 
   The routine will terminate when it believes it has achieved either the requested 
   absolute error or requested relative error (whichever is the easier) or when it 
   determines that neither accuracy can be obtained (due to accumulation of rounding
   error or the maximum number of function evaluations being reached).
   
   The following are returned or are directly accessible via access functions:
	 - the approximation to the given integral
	 - the estimated relative error in the approximation
	 - an error flag indicating how reliable the returned result is likely to be, and
	 - the number of function evaluations actually used in the computation.
 
   Although integrals of type (2) above include those of type (1), it is strongly
   recommended that the calling sequence designed for type (1) integrals be used
   whenever the integration domain is cuboid.  The routine has inbuilt efficiencies
   for this case.
 
   Details of the underlying algorithm can be found in Li T. and Robinson I., 
   "elrint3d: A Three-Dimensional Non-Adaptive Automatic Cubature Routine Using a 
   Sequence of Embedded Lattice Rules", ACM Transactions on Mathematical Software,
   to appear.
 
   Further details on how to use the software, including sample programs, are
   provided in the document associated with this package.
 
   Notes
	1. The name 'elrint3d' is an acronym for [E]mbedded [L]attice [R]ule
	   [Int]egrator - [3] [D]imensions.
	2. Outside this file, 'elrint3d' can only be accessed through a defined access
	   function.
	3. The algorithm is designed specifically for 64-bit double precision arithmetic
	   and includes several 'hard-wired' design decisions (for example, the lattice
	   sequence used, the maximum number of function evaluations permitted, and so on).

 ****************************************************************************************/

// Standard includes
#include    <algorithm>
#include    <limits>
#include    <cfloat>

// Application includes
#include    <infiniteMap.h>
#include    <elrint3dSequence.h>
#include    <embeddedCubRule.h>

// Constants
#define DEFAULT_ABS_TOL    0.0     // Default requested absolute accuracy - maximum accuracy
#define DEFAULT_REL_TOL    0.0     // Default requested relative accuracy - maximum accuracy
#define DEFAULT_MAX_EVALS  421376  // Default maximum number of function evaluations permitted
#define MAX_SEQ            10      // Maximum number of lattices in the augmentation sequence

// Define infinity 
#ifndef INFINITY
#define INFINITY  std::numeric_limits<double>::infinity()  // Double precision infinity value
#endif

// Function pointer definitions
typedef double (*BOUND_FUNC_1) (double); // Single-variable boundary function - used for y coordinate
typedef double (*BOUND_FUNC_2) (double, double); // Boundary function of two variables - used for z coordinate
typedef double (*FUNC_3) (double, double, double); // Function pointer for the integrand
typedef double (*INFINITE_MAP) (double, double, double&); // Infinite map

class Elrint3d : public EmbeddedCubRule<double>
{
public:
	// Define error flags indicating the cubature status

	enum
	{
		EF_LOW = -1, // Invalid status, used as initial status
		NORMAL, // Requested accuracy is met
		NOT_SUITABLE, // The algorithm is not suited to this integration - terminated early
		NOT_IMPROVEMENT, // Rounding error prevents further improvement - terminated early
		REACH_MAX_NUMBER // Maximum number of function evaluations reached
	};

private:
	// Wrapper acts as a converter from a function pointer to an integrand object.  Thus, internally, only the
	// class definition is used for the integrand even if it was defined using a function pointer externally.

	class Integrand3DWrapper : public Integrand<double>
	{
	private:
		FUNC_3 pf; // Function pointer
	public:

		double fun(const double x[]) const
		{
			return pf(x[0], x[1], x[2]);
		}

		// Constructors

		Integrand3DWrapper(FUNC_3 f) : pf(f)
		// Passing in the function pointer
		{
		}

		Integrand3DWrapper(const Integrand3DWrapper& iwa) : pf(iwa.pf)
		// Complying with the compiler option -Weffc++
		{
		}
	};

	const Integrand3DWrapper intWrapper;
	// Used to wrap a function pointer in a class definition when
	// a function pointer is used to define the integrand

	const Integrand<double>& rIntegrand;
	// Points to 'intWrapper' when a funtion pointer is used to define the integrand

	// Integration domain definitions:
	//   1. Cuboid
	//   2. Variable boundary
	//   3. Semi-infinite or infinite

	class BoundaryDef
	{
	private:
		// Integration domain for integrals of type (1) above
		const double cx1, cx2; // x interval
		const double cy1, cy2; // y interval
		const double cz1, cz2; // z interval

		// Integration domain for type (2) integrals
		// 'cx1', 'cx2' used for x interval
		const BOUND_FUNC_1 fy1, fy2; // y interval
		const BOUND_FUNC_2 fz1, fz2; // z interval

		// Internally, we need to remember the type of integration domain
		const bool byCons; // True if the y interval is specified by constants
		const bool bzCons; // True if the z interval is specified by constants

		// Internally, we need to remember if the integration domain is infinite
		bool bxInfty; // True if the x interval is semi-infinite or infinite
		bool byInfty; // True if the y interval is semi-infinite or infinite
		bool bzInfty; // True if the z interval is semi-infinite or infinite

		INFINITE_MAP inftyMap; // Only used when there is an infinite boundary

	public:
		double trans(double* px, double* py, double* pz) const;
		// Transforms the coordinates from the unit cube to the original domain
		// and returns the weight associated with the point (px, py, pz)

		bool bInfinite() const
		// Checks whether the boundary is infinite
		{
			return bxInfty || byInfty || bzInfty;
		}

		void setInfMapping(INFINITE_MAP map)
		// Changes the infinite mapping
		{
			inftyMap = map;
		}

		// Constructors for different integration domain types
		BoundaryDef(double x1, double x2, double y1, double y2, double z1, double z2, INFINITE_MAP infMap = rationalInfMap<double>);
		// Boundary specified by constants

		BoundaryDef(double x1, double x2, BOUND_FUNC_1 y1, BOUND_FUNC_1 y2, double z1, double z2, INFINITE_MAP infMap = rationalInfMap<double>);
		// Variable y boundary and other boundaries constant


		BoundaryDef(double x1, double x2, double y1, double y2, BOUND_FUNC_2 z1, BOUND_FUNC_2 z2, INFINITE_MAP infMap = rationalInfMap<double>);
		// Variable z boundary and other boundaries constant


		BoundaryDef(double x1, double x2, BOUND_FUNC_1 y1, BOUND_FUNC_1 y2, BOUND_FUNC_2 z1, BOUND_FUNC_2 z2, INFINITE_MAP infMap = rationalInfMap<double>);
		// Variable y and z boundaries
	};

	BoundaryDef bound; // Boundary definition

	const double absReqTol; // Requested absolulte tolerance
	const double relReqTol; // Requested relative tolerance

	const int maxEvals; // Maximum number of function evaluations permitted

	// Cubature information
	double estError; // Estimated error - set to maximum value initially, then updated after each sequence advancement
	double evalSum; // Running sum of function values
	double interCubs[MAX_SEQ]; // Intermediate cubatures
	int eflag; // Error flag indicating cubature status
	bool bStop; // True if the algorithm needs to be terminated
	int nEvals; // Number of function evaluations - used only when the integration domain is infinite
	bool bFirstNoImprovement; // True if little improvement can be achieved by continuing the algorithm

public:
	double evaluate();
	// Generates the sequence of cubatures and returns the final approximation

	int errFlag() const
	// Returns a value indicating the integration status
	// (This function overrides the virtual function in 'EmbeddedCubRule')

	// Possible values:
	// 0 - Normal termination.  The requested accuracy is assumed to have been met.
	// 1 - The algorithm is not suited to this integral.  For example, the integrand
	//     is not smooth even after application of the periodizing transformation.
	// 2 - Rounding error prohibits further improvement in the approximation.
	// 3 - The maximum number of function evaluations has been reached.

	// For each case in which 'errFlag' returns a non-zero value, it is believed
	// that the requested accuracy has not been met.
	{
		return eflag;
	}

	double estErr() const
	// Returns the estimated relative error
	{
		return estError;
	}

	int evals() const
	// Returns the number of function values used
	{
		if (bound.bInfinite())
			// The 'evaluate' function may have been invoked twice, in which case the
			// total number of function evaluations is stored internally in 'nEvals'
		{
			return nEvals;
		}
		else
			// Finite domain; access 'evals' directly
		{
			return EmbeddedCubRule<double>::evals();
		}
	}

protected:
	void postAdvancement(int iLat);
	// Computes the error estimate and tests for termination of the algorithm (either
	// because the requested accuracy has been achieved or because it cannot be achieved)

	bool bTerminated() const;
	// Checks if the algorithm should be terminated

	double funcEval(const double* pa) const;
	// Evaluates th eintegrand at the point 'pa'

	double trans(double * const pa) const;
	// Transforms the coordinates from the unit cube to the original domain
	// and returns the weight associated with the point (px, py, pz)
	// (This function overrides the virtual function in 'EmbeddedCubRule')

	double sequenceAdvance(int iLat);
	// Stores information used in 'postAdvancement(iLat)'
	// (This function overrides the virtual function in 'EmbeddedCubRule')

	double sum(const double* pa, int n) const;

	double errPerCycle(const double e[]) const;
	// Computes an error estimate based on the difference between successive cubatures in a
	// single augmentation cycle - only called immediatley after the first augmentation cycle

public:
	void reset();
	// Resets the routine so it calls 'evaluate' multiple times without creating a new instance each time

	// Constructors

	template<typename BD1, typename BD2>
	Elrint3d(const Integrand<double>& integrand,
			double x1, double x2, // x interval
			BD1 y1, BD1 y2, // y interval
			BD2 z1, BD2 z2, // z interval
			double rel_tol = DEFAULT_REL_TOL, // Relative tolerance
			double abs_tol = DEFAULT_ABS_TOL, // Absolute tolerance
			int max_evals = DEFAULT_MAX_EVALS) // Maximum number of function evaluations
	: EmbeddedCubRule<double>(singleton_elrint3dSequence()),
	intWrapper(NULL), rIntegrand(integrand),
	bound(x1, x2, y1, y2, z1, z2),
	absReqTol(abs_tol), relReqTol(rel_tol), maxEvals(max_evals),
	estError(DBL_MAX), evalSum(0.0), eflag(EF_LOW), bStop(false), nEvals(0), bFirstNoImprovement(false)
	// Class definition for integrand - BD1 and BD2 can be constant or BOUND_FUNC_1 or BOUND_FUNC_2
	{
	}

	template<typename BD1, typename BD2>
	Elrint3d(FUNC_3 f,
			double x1, double x2, // x interval
			BD1 y1, BD1 y2, // y interval
			BD2 z1, BD2 z2, // z interval
			double rel_tol = DEFAULT_REL_TOL, // Relative tolerance
			double abs_tol = DEFAULT_ABS_TOL, // Absolute tolerance
			int max_evals = DEFAULT_MAX_EVALS) // Maximum number of functional evaluations
	: EmbeddedCubRule<double>(singleton_elrint3dSequence()),
	intWrapper(f), rIntegrand(intWrapper),
	bound(x1, x2, y1, y2, z1, z2),
	absReqTol(abs_tol), relReqTol(rel_tol), maxEvals(max_evals),
	estError(DBL_MAX), evalSum(0.0), eflag(EF_LOW), bStop(false), nEvals(0), bFirstNoImprovement(false)
	// Function pointer definition for integrand - BD1 and BD2 can be constant or BOUND_FUNC_1 or BOUND_FUNC_2
	{
	}

};

#endif  // _ELRINT3D_H_

