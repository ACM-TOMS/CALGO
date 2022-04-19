#ifndef _INFINITE_UNIT_MAP_H_
#define	_INFINITE_UNIT_MAP_H_

/*******************************************************************************
   
   This module provides two types of transformation from a semi-infinite or 
   infinite interval to the interval [0,1], namely a rational transformation 
   and a logarithmic transformation.  Each is generic, to allow for different 
   arithmetic precisions.
   
   The effectiveness of these transformations depends on the particular integral
   being evaluated.  No general advice is provided regarding preference for one 
   transformation over the other.

 *******************************************************************************/

// Standard includes
#include    <algorithm>
#include    <cmath>

template<typename T>
T rationalInfMap(const T p, const T q, T& x)
// Implements the rational transformation.  
// (p,q) is the integration interval, with one or both of p and q infinite.
// x is the point to be mapped.  On entry, it is a value between p and q.  
// On exit, it has the corresponding mapped value between 0 and 1.
// The function returns the weight associated with x.

{
#define INFTY std::numeric_limits<T>::infinity()

	T coef = 1.0, u, v;
	v = x;

	if (-INFTY < p && p < INFTY && INFTY <= q)
		// Semi-infinite
	{
		if (v != 0.0)
		{
			coef *= 1.0 / (v * v);
			u = p - 1.0 + 1.0 / v;
		}
		else
		{
			u = 0.0;
			coef = 0.0;
		}
	}
	else
	{
		if (p <= -INFTY && -INFTY < q && q < INFTY)
			// Semi-infinite
		{
			if (v != 0.0)
			{
				coef *= 1.0 / (v * v);
				u = q + 1.0 - 1.0 / v;
			}
			else
			{
				u = 0.0;
				coef = 0.0;
			}
		}
		else
			// Infinite
		{
			if (v == 0.0 || v == 1.0)
			{
				coef = 0.0;
				u = 0.0;
			}
			else
			{
				coef *= (1.0 / ((1.0 - v)*(1.0 - v)) + 1.0 / (v * v));
				u = 1.0 / (1.0 - v) - 1.0 / v;
			}
		}
	}
	x = u;
	return coef;
#undef INFTY
}

template<typename T>
T logarithmicInfMap(const T p, const T q, T& x)
// Implements the logarithmic transformation.
// (p,q) is the integration interval, with one or both of p and q infinite.
// x is the point to be mapped.  On entry, it is a value between p and q.  
// On exit, it has the corresponding mapped value between 0 and 1.
// The function returns the weight associated with x.
{
#define INFTY std::numeric_limits<T>::infinity()

	T u, coef = 1.0;
	T v = x;
	if (-INFTY < p && p < INFTY && INFTY <= q)
		// Semi-infinite
	{
		if (v != 1.0)
		{
			coef *= 1.0 / (1.0 - v);
			u = p + std::log(1.0 / (1.0 - v));
		}
		else
		{
			u = 0.0;
			coef = 0.0;
		}
	}
	else
	{
		if (p <= -INFTY && -INFTY < q && q < INFTY)
			// Semi-infinite
		{
			if (v != 1.0)
			{
				coef *= 1.0 / (1.0 - v);
				u = q - std::log(1.0 / (1.0 - v));
			}
			else
			{
				u = 0.0;
				coef = 0.0;
			}
		}
		else
			// Infinite
		{
			if (v == 0.0 || v == 1.0)
			{
				coef = 0.0;
				u = 0.0;
			}
			else
			{
				coef *= 1.0 / (v * (1.0 - v));
				u = log(v / ((1.0 - v)));
			}
		}
	}
	x = u;
	return coef;
#undef INFTY
}

#endif  // _INFINITE_UNIT_MAP_H_

