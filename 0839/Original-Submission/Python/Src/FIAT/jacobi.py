# does several operations related to the jacobi polynomials

import math

def make_jacobi(a,b,n):
	if n == 0:
		return func.one
	else:
		return JacobiP(a,b,n)

# evaluates the nth jacobi polynomial with weight parameters a,b at a
# point x 
def eval_jacobi(a,b,n,x):
    if 0 == n:
        return 1.0;
    elif 1 == n:
        return 0.5 * ( a - b + ( a + b + 2.0 ) * x )
    else: # 2 <= n
        apb = a + b
        pn2 = 1.0
        pn1 = 0.5 * ( a - b + ( apb + 2.0 ) * x )
        p = 0
        for k in range(2,n+1):
            a1 = 2.0 * k * ( k + apb ) * ( 2.0 * k + apb - 2.0 )
            a2 = ( 2.0 * k + apb - 1.0 ) * ( a * a - b * b )
            a3 = ( 2.0 * k + apb - 2.0 ) * ( 2.0 * k + apb - 1.0 ) * ( 2.0 * k + apb )
            a4 = 2.0 * ( k + a - 1.0 ) * ( k + b - 1.0 ) * ( 2.0 * k + apb )
            a2 = a2 / a1
            a3 = a3 / a1
            a4 = a4 / a1
            p = ( a2 + a3 * x ) * pn1 - a4 * pn2
            pn2 = pn1
            pn1 = p
        return p

# evaluates the derivative of the jacobi polynomial
def eval_jacobi_deriv(a,b,n,x):
	if n == 0:
		return 0.0
	else:
		return 0.5 * ( a + b + n + 1 ) * eval_jacobi(a+1,b+1,n-1,x)

# uses Newtons method with Chebyshev points as a starting guess.  This
# algorithm is implemented from the appendix of Karniadakis and Sherwin.
def compute_gauss_jacobi_points( a , b , m ):
	x = []
	eps = 1.e-8
	max_iter = 100
	for k in range(0,m):
		r = -math.cos(( 2.0*k + 1.0) * math.pi / ( 2.0 * m ) )
		if k > 0:
			r = 0.5 * ( r + x[k-1] ) 
		j = 0
		delta = 2 * eps
		while j < max_iter:
			s = 0
			for i in range(0,k):
				s = s + 1.0 / ( r - x[i] )
			f = eval_jacobi(a,b,m,r)
			fp = eval_jacobi_deriv(a,b,m,r)
			delta = f / (fp - f * s)

			r = r - delta

			if math.fabs(delta) < eps:
				break
			else:
				j = j + 1

		x.append(r)
	return x

# again taken from Karniadakis and sherwin - the gauss-lobatto points
# include the endpoints.
def compute_gauss_lobatto_jacobi_points( a , b , n ):
	return [ -1.0 ] + compute_gauss_jacobi_points( a + 1 , b + 1 , n - 2 ) + [ 1.0 ]
