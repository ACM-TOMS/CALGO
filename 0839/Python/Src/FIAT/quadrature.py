import jacobi, math, gamma, points
from Numeric import *
import operator

# Gauss-Jacobi quadrature rules in one dimension and
# for collapsed coordinates on triangles.

# constructor takes a 2-tuple consisting of points and weights.
# notice that loose typing of Python allows us to inherit quadrature
# rules for arbitrary number of dimensions from it.
# Provides methods for obtaining the points, weights, and a member
# function for integrating that takes a callable object and returns
# its integral approximated by the quadrature rule.
# integrate can work either on a callable object that works on lists
# or else an array of already tabulated data of the correct length.
# The array type must support componentwise multiplication

class QuadratureRule:
	def __init__(self,rule):
		self.x = rule[0]
		self.w = array( rule[1] )
	def get_points(self):
		return self.x
	def get_weights(self):
		return self.w
	def __call__( self , f ):
		return self.integrate( f )
	def integrate(self,f):
		fs = array( [ f( x ) for x in self.x ] )
		return sum( self.w * fs )

# takes a 1-d integration rule and an edge_number (0,1,2) and returns the rule that integrates
# a function of two variables along that edge (counterclockwise)
def to_edge_quad(ed_no,line_rule):
	if ed_no < 0 or ed_no > 2:
		raise RuntimeError, "Illegal edge number in quadrature.to_edge_quad"
	if ed_no == 0: # hypotenuse is longer than [-1,1] , so include Jacobian of mapping
		ws = math.sqrt(2.0) * array( line_rule.get_weights() )
	else:
		ws = array( line_rule.get_weights() )
	xs = tuple( points.interval_to_edge_points( ed_no , line_rule.get_points() ) )
	return QuadratureRule( ( xs, ws ) )


# one-d gauss-jacobi quadrature.  User inputs parameters
# (a,b) which describe the weights for the Jacobi polynomials
# and m, the number of points in the quadrature rule
# the points are the zeros of P_m^{a,b} on [-1,1], and the
# weights are as given in standard numerical analysis texts.
# We took them from the appendix in Karniadakis & Sherwin
# _Spectral/hp_Element_Methods_for_CFD_, appendix B

class JacobiQuadrature(QuadratureRule):
	def __init__(self,a,b,m):
		QuadratureRule.__init__(self,make_jacobi_quadrature(a,b,m))

class JacobiQuadrature2D(QuadratureRule):
	def __init__(self,m):
		QuadratureRule.__init__(self,make_jacobi_quadrature_collapsed(m))

# returns (x,w), where x is the list of quadrature points and w is
# the list of quadrature weights.  Algorithm taken from 
# Karniadakis & Sherwin, Appendix B

def make_jacobi_quadrature(a,b,m):
	x = jacobi.compute_gauss_jacobi_points( a , b , m )

	w = []
	a1 = math.pow(2,a+b+1);
	a2 = gamma.gamma(a + m + 1);
	a3 = gamma.gamma(b + m + 1);
	a4 = gamma.gamma(a + b + m + 1);
	a5 = fact(m);
	a6 = a1 * a2 * a3 / a4 / a5;

	for k in range(0,m):
		fp = jacobi.eval_jacobi_deriv(a,b,m,x[k])
		w.append( a6 / ( 1.0  - x[k]**2.0) / fp ** 2.0 )

	return (tuple(x),array(w))

# computes tensor product of two Jacobi quadrature rules on reference
# square, then maps to the reference triangle via the collapsed coordinate
# mapping.
# returns the pair containing the list of points and the list of weights.
def make_jacobi_quadrature_collapsed(m):
	qx = JacobiQuadrature(0,0,m)
	qy = JacobiQuadrature(1,0,m)
	w = array( [ 0.5 * w1 * w2 for w1 in qx.w for w2 in qy.w ] )
	x = tuple( map( collapse_coords, \
				     [ (x,y) for x in qx.x \
				             for y in qy.x ] ) )
	return (x,w)


def fact(m):
	result = 1
	for i in range(2,m+1):
		result = result * i
	return result
	
def collapse_coords(c):
	(x,y) = c
	return (0.5 * ( 1.0 + x ) * (1.0 - y) - 1.0 , y )

