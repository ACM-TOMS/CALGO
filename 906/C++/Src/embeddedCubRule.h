#ifndef _EMBEDDED_CUB_RULE_H_
#define _EMBEDDED_CUB_RULE_H_

/*******************************************************************************
                                                                   
   'EmbeddedCubRule' is a generic class that encpsulates the application of
   a cubature rule to a sequence of embedded lattices.  Functions associated 
   with termination of the sequence are included.  To develop a specific
   application, users must extend the class by overriding pure virtual
   methods and non-pure virtual methods, as appropriate.     
                                                                                
   The class is designed in a way that makes it possible for the user to
   experiment with different dimensions, different precision requirements (e.g.,
   float, double, long double) and paralellism.                                              
                                                                                
   Implementation                                                        
   To extend this class for a specific cubature rule based on a constructed   
   embedded sequence, the following pure virtual functions must be implemented
   in the inherited class:                                                                       
   
   1. T trans(const T*) const
	  Maps a point in the integration domain (including infinite domains) to
	  [0,1]^s.   A periodizing transformation may be applied or may have already
		  been applied.  Returns the weight associated with the point.
   2. T funcEval(T*) const
	  Evaluates the integrand function at points in the original integration
	  domain, i.e., after application of the 'trans' function.
   3. T sequenceAdvance(int iLat) const
	  Returns the weighted sum of function values at the (transformed) points
	  in the lattice indexed by 'iLat', excluding previously evaluated points.
   4. void postAdvancement(int iLat);
	  Post-processing following 'sequenceAdvance', in which the estimated
	  error is calculated and a termination flag is set (including 'continue').
   5. T evaluate()
	  Generates the sequence of cubatures until it is believed that either the
	  requested accuracy has been achieved or it cannot be achieved using the
	  given sequence of lattices.
   6. bool bTerminated() const
	  True if the algorithm should be terminated.
   7. int errFlag() const
	  Returns an error flag indicating the outcome of the algorithm.
   8. T estErr() const
	  Returns the estimated relative error for the current cubature.
   9. int evals() const
	  Returns the number of function evaluations to date.
  10. T sum(const T*, int) const
	  Summation method - uses recursive (accumulative) summation, but users
	  may wish to re-implement the function using a different summation method.
  11. void reset()
	  Resets the algorithm so it can restart.  It is recommended that users
	  re-implement this function to reset the internal status of their own
	  algorithm.

   Notes
   1. This is an abstract class; users must re-implement the class to suit   
	  their own algorithm.
   2. The re-implementation requires an associated embedded lattice sequence.      
   3. The class provides some well-defined implementations for the 'evaluate'
	  and 'sequencAdvance' methods.  Users can override these implementations
	  or the implementations can be re-used by calling explicitly
	  EmbeddedCubRule<T>::sequenceAdvance() and EmbeddedCubRule<T>::evaluate().

 ********************************************************************************/

// Application includes
#include	<embeddedSequence.h>
#include	<integrand.h>

template <typename T>
class EmbeddedCubRule
// T is the data type for the design cubature rule 
// corresponding to the arithmetic precision desired
{
private:
	const EmbeddedSequence<T> * const pSequence; // The provided embedded lattice sequence
	int iLat; // Index of current embedded lattice

public:
	virtual T evaluate()
	// Generates the sequence of cubatures and returns the final approximation.

	// 1. While the termination flag is false
	//    1.1 Call 'sequenceAdvance' to compute the sum of function values
	//        needed for the current cubature
	//    1.2 Compute the sum and thence the cubature
	//    1.3 Call 'postAdvancement' to calculate the error estimate and set
	//        the termination flag
	// 2. Return the approximation to the integral
	{
		T sum = 0;
		T value = 0;

		while (true)
		{
			sum += sequenceAdvance(iLat);
			value = sum / pSequence->count(iLat);
			postAdvancement(iLat);
			if (bTerminated())
				// Either the requested accuracy has been achieved
				// or an abnormal termination flag has been set
			{
				break;
			}
			iLat++;
		}
		return value;
	}

	const EmbeddedSequence<T>* sequence() const
	// Returns the sequence used for the cubature rule
	{
		return pSequence;
	}

	virtual int errFlag() const = 0;
	// Returns a flag indicating the outcome of the integration.
	// The value 0 corresponds to normal termination (i.e., the requested accuracy
	// is believed to have been met).
	// A non-zero value is returned when an abnormal termination occurs. The
	// particular values used, and the condition they indicate, depend on the
	// specific implementation in the inherited class.

	virtual T estErr() const = 0;
	// Returns the estimated relative error

	virtual int evals() const
	// Returns the total number of function evaluations to date.  This function must
	// be overridden if more than one attempt is made to evaluate a given integral,
	// for instance by invoking different mappings of semi-infinite or infinite
	// domains onto [0,1]^s.
	{
		return pSequence->count(iCurLat());
	}

	int iCurLat() const
	// Returns the index of the current lattice
	{
		return iLat;
	}

	virtual void reset()
	// Resets the algorithm so that evaluation can start from the first lattice
	{
		iLat = 0;
	}

protected:

	virtual void postAdvancement(int iLat) = 0;
	// Computes the error estimate and sets the termination flag

	virtual bool bTerminated() const = 0;
	// True if the algorithm should be terminated

	virtual T funcEval(const T * pa) const = 0;
	// Evaluates the integrand function at points in the original integration
	// domain, i.e., after application of the 'trans' function

	virtual T trans(T * const pa) const = 0;
	// Maps a point in the integration domain to [0,1]^s.  A periodizing transformation
	// may be applied within the method or have already been applied.  Returns the weight
	// associated with the point.

	virtual T sum(const T* pa, int n) const
	// Straightforward summation algorithm.  Re-implement using an alternative method if
	// rounding error is or may be an issue.
	{
		T total = 0;
		for (int i = 0; i < n; i++)
		{
			total += pa[i];
		}
		return total;
	}

	virtual T sequenceAdvance(int iLat)
	// Returns the weighted sum of function values at the (transformed) points in the
	// lattice indexed by 'iLat', excluding previously evaluated points.
	// Users can override this function by calling 'EmbeddedCubRule::evaluationSum'
	// if they wish to save the intermediate cubature following the advancement of
	// the sequence.
	// The function could also be re-implemented if a parallel version is desired,
	// running multiple threads on a multi-core machine.
	{
		int cp = pSequence->count(iLat);

		if (iLat > 0)
		{
			cp -= pSequence->count(iLat - 1);
		}

		T values[cp]; // Storage of the weighted sum of function values

		typename EmbeddedSequence<T>::Iterator its = pSequence->iStartAtSeq(iLat);

		for (int i = 0; i < cp; i++)
		{
			T* coord = *its;
			T weight = trans(coord); // Map to original domain
			values[i] = funcEval(coord) * weight; //  Add to weighted sum
			its++;
		}
		return sum(values, cp);
	}

public:
	// Constructors
	EmbeddedCubRule(const EmbeddedSequence<T>* pLat)
	: pSequence(pLat), iLat(0)
	// Constructor with embedded sequence and integrand specified
	{
	}

	EmbeddedCubRule(const EmbeddedCubRule<T>& ecr)
	: pSequence(ecr.pSequence), iLat(ecr.iLat)
	// Same as default copy constructor, but complies with compiler option -Weffc++
	{
	}

	virtual ~EmbeddedCubRule()
	// Virtual destructor required
	{
	}

	// Assignment operator
	EmbeddedCubRule & operator=(const EmbeddedCubRule<T> & ecr)
	// Same as default assignment operator, but complies with compiler option -Weffc++
	{
		pSequence = ecr.pSequence;
		iLat = ecr.iLat;
	}
};

#endif  // _EMBEDDED_CUB_RULE_H_
