#ifndef	_CUB_UTIL_H_
#define _CUB_UTIL_H_

/**********************************************************

  This module contains frequently-used global functions

 ***********************************************************/

// Standard includes
#include	<algorithm>
#include	<cmath>
#include    <cfloat>

using namespace std;

template<typename T>
T max3(T a, T b, T c)
// Returns the maximum among three
{
	return max(max(a, b), c); // Use stl max
}

template<typename T>
T min3(T a, T b, T c)
// Returns the minimum among three
{
	return min(min(a, b), c); // Use stl min
}

template <typename T>
T abs_err(T act, T eval)
// Returns the absolute error estimate
// 'act' is the actual value
// 'eval' is the approximation
{
	return abs(eval - act);
}

template<typename T>
T rel_err(T act, T eval)
// Returns the relative error estimate
// 'act' is the actual value
// 'eval' is the approximation
{
	return act == 0.0 ? abs_err(act, eval) : abs_err(act, eval) / abs(act);
}

template<typename T>
T kahan_sum(const T* data, int n)
// Kahan summation algorithm, also known as compensated summation
{
	T q = data[0];
	T r = 0;
	T s;
	T t;
	for (int i = 1; i < n; i++)
	{
		s = data[i] - r;
		t = q + s;
		r = (t - q) - s;
		q = t;
	}
	return q;
}

template<typename T>
T epsilon_extrapolation(T r1, T r2, T r3)
// Applies the epsilon algorithm to compute an extrapolated value based
// on the three most recent values in a 2's-copy augmentation sequence
{
	T temp, t1, t2;
	temp = r2 - r1;
	if (temp == 0.0)
	{
		t1 = DBL_MAX;
	}
	else
	{
		t1 = 1.0 / temp;
	}
	temp = r3 - r2;
	if (temp == 0.0)
	{
		t2 = DBL_MAX;
	}
	else
	{
		t2 = 1.0 / temp;
	}
	temp = t2 - t1;
	if (temp == 0.0)
	{
		return r3;
	}
	else
	{
		return r2 + 1.0 / temp;
	}
}

#endif  // _CUB_UTIL_H_
