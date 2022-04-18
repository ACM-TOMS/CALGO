#ifndef _INTEGRAND_H_
#define _INTEGRAND_H_

/*****************************************************************************
                                                               
   Generic template class for defining the integrand. The template type
   can be used to handle different arithmetic precisions.                          
                                                                                
   Defining the integrand in this way provides more flexibility than the 
   standard function definiton, especially when the integrand cannot be
   expressed easily in analytic form.                                           
                                                                                
   Implementation                                                       
   To define the intended integrand, users must implement 'fun(const T[])'
   in the subclass. 

 ******************************************************************************/

template <typename T>
class Integrand
{
public:
	virtual T fun(const T x[]) const = 0;
	// Returns the function value based on the given coordinates 'x'
	// regardless of dimensionality.

	virtual ~Integrand()
	// Virtual destructor - required for an abstract class
	{
	}
};

#endif  // _INTEGRAND_H_
