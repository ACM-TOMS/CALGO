// Author: Yong-Kang Zhu (yongkang.zhu@gmail.com)
// Code can be used only for academic purpose

// ExactSum.h and ExactSum.cpp provide algorithms for
// correctly rounded floating-point summation

// GCC/G++ Compilation Options:
// -DDOUBLE: set the rounding mode to Double precision
// -DREV: reverse the floating-point representation structure for SPARC and PowerPC machines

// How to compile?
// (1) x86 Linux with GCC/G++:  -O1 -DDOUBLE
// (2) x86 Windows with Virtual C++: Win32, Release
// (3) Cygwin Windows with GCC/G++: -O1 -DDOUBLE
// (4) SPARC Solaris with GCC/G++: -O1 -DREV
// (5) PowerPC with GCC/G++: -O1 -DREV
// (6) x86 with Mac OS X: -O1

// the number of exponents for IEEE754 double
#define N_EXPONENT 2048

// the length of the accumulators, i.e., 2 X N_EXPONENT
#define N2_EXPONENT (2*N_EXPONENT)

#define HALF_MANTISSA 26 // number of bits in a split mantissa

// Max number of split mantissas that can be summed without error
#define MAX_N (1 << HALF_MANTISSA) // 2^HALF_MANTISSA
#define MAX_N_AFTER_SWAP (MAX_N - N2_EXPONENT)

// a structure for IEEE754 double precision
struct str_double
{
#ifdef REV
// floating-point number representation on SPARC and PowerPC ...
	unsigned sign : 1;
	unsigned exponent : 11;
	unsigned mantissa_high  : 20;
	unsigned mantissa_low : 32;
#else
	unsigned mantissa_low : 32;
	unsigned mantissa_high  : 20;
	unsigned exponent : 11;
	unsigned sign : 1;
#endif
};

// Computes a sum which is guaranteed-accurate
// Note: not synchronized
class ExactSum
{
private:
	int r_c;      // a global used by iFastSum

	// Part B
	int c_num;    // number of the summands which have been accumulated

	// accumulators of OnlineExactSum used by Part B
	double *t_s, *t_s2;

	// Return 1 if not correctly rounded; 0 otherwise
	int Round3(double s0, double s1, double s2);

	// set the rounding mode to the double precision
  void set_fpu (unsigned int mode);

public:
	ExactSum();
	~ExactSum();

	// Dekker's algorithm: Exact Addition for 2 numbers
	void AddTwo(double &a, double &b);

	// Part A: Sum the given array (more efficient than Part B)

	// Returns a correctly rounded sum
	// Note: a. the array starts from [1]
	//       b. after execution, num_list is destroyed
	double iFastSum(double *num_list, int n);

	// Returns a correctly rounded sum
	// Note: a. the array starts from [1]
	//       b. after execution, num_list does not change
	//       c. to satisfy b., it does not call iFastSum if n < 2000
	// Comment: since iFastSum is empirically faster when n < 2000,
	//          users can call iFastSum if n < 2000, and OnlineExactSum otherwise
	double OnlineExactSum(double *num_list, int n);


	// Part B: Online Summation if summands are not given at once, i.e.,
	//         users can feed a number or an array, and get a sum at any time
	// Note: a. use OnlineExactSum to add either a number or an array, and
	//       b. return the result if GetSum() is called

	// Adds one number, and is used with GetSum()
	void AddNumber(double x);

	// Adds an array, and is used with GetSum()
	void AddArray(double *num_list, int n);

	// Returns the current sum
	// Also see AddNumber()
	double GetSum();

	// Resets sum to zero
	void Reset();
};
