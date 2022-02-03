import func, jacobi, math
from curry import curry

from func import MemoizedFunction


# Dubiner functions are D_{i,j}(x,y) =
# c L_i(eta_1) (1-eta-2)^i P_j^{2i+1,0},
# where eta_1 and eta_2 are functions of x,y
# these are orthogonal and hierarchical functions
# spanning P_k on the reference triangle


def dub_dx( i , j , x ):
    """Evaluates the x partial derivative of D_{i,j} at a point x"""	
    if i == 0:
        return 0.0
    else:
        alpha = math.sqrt( (i+0.5)*(i+j+1)) \
                * math.pow(0.5,i) \
                * ( 1.0 + i )
        ( a , b ) = coord_transform( x )
        f1 = jacobi.eval_jacobi( 1 , 1 , i - 1 , a )
        f2 = jacobi.eval_jacobi( 2*i+1 , 0 , j , b ) \
             * math.pow(1-b , i - 1 )
        return alpha * f1 * f2

def dub_dy( i , j , x ):
    """Evaluates the y partial derivative of D_{i,j} at a point x"""	
    (a,b) = coord_transform(x)
    fa = jacobi.eval_jacobi(0,0,i,a)
    dfa = jacobi.eval_jacobi_deriv(0,0,i,a)
    gb = jacobi.eval_jacobi(2*i+1,0,j,b)
    dgb = jacobi.eval_jacobi_deriv(2*i+1,0,j,b)
    dds = dfa * gb * 0.5 * (1 + a)
    if i>0:
        dds *= math.pow(0.5*(1-b),i-1)
    tmp = dgb * math.pow(0.5 * (1-b), i)
    if i>0:
        tmp -= 0.5 * i * gb * math.pow(0.5 * (1-b), i-1)
    dds += fa * tmp
    dds *= math.sqrt((i+0.5)*(i+j+1))
    return dds

# curry dub_dx and dub_dy to get a function of just location for any
# pair i,j
def dubdx( i , j ):
    return curry( dub_dx , i , j )

def dubdy( i , j ):
    return curry( dub_dy , i , j )

class Dubiner:
    """Instances implement a Dubiner polynomial D_{i,j}.
They are callable and supply differentiation methods"""
    def __init__( self , i , j ):
        self.i = i
        self.j = j
        self.alpha = math.sqrt((i+0.5)*(i+j+1)) * (0.5 ** i)
    def __call__( self , x ):
        ( a , b ) = coord_transform( x )
        return self.alpha \
               * jacobi.eval_jacobi( 0 , 0 , self.i , a ) \
               * jacobi.eval_jacobi( 2*self.i + 1 , 0 , self.j , b ) \
               * ( ( 1.0 - b ) ** self.i )
    def __str__( self ):
        return "D_{" + str( self.i ) + "," + str( self.j ) + "}"
    def deriv( self , dir ):
	if dir == 0:
		return dubdx( self.i , self.j )
	elif dir == 1:
		return dubdy( self.i , self.j )
	else:
		raise RuntimeError, "illegal direction"


# don't recreate multiple instances of the dubiner functions.
# I keep a list of all dubiner instances I've created.
# the make_dubiner_funcs function returns this list, a subset thereof,
# or adds new terms to the list and returns it.  That way, if
# I create dubiner functions for the function space, the constraints
# on the function space, and the nodes for the finite element, I get
# to reuse the cache if possible in each instance.
current_dubiners = [ MemoizedFunction( Dubiner( 0 , 0 ) ) ]
current_dubiner_degree = 0
def make_dubiner_funcs( n ):
    global current_dubiners
    global current_dubiner_degree
    if n <= current_dubiner_degree:
        dimpn = (n+1)*(n+2)/2
        return current_dubiners[:dimpn]
    else:
        current_dubiners = current_dubiners + \
               [ MemoizedFunction( Dubiner( k - i , i ) ) \
                     for k in range(current_dubiner_degree+1,n+1) \
                     for i in range( 0 , k + 1 ) ] 
        current_dubiner_degree = n
        return current_dubiners

# returns a list of Dubiner polynomials.  called in
# functionspace.PolynomialTriangle, among other places
#def make_dubiner_funcs( n ):
#	return [ MemoizedFunction( Dubiner( k - i , i ) ) \
#		for k in range( 0 , n + 1 ) for i in range( 0 , k + 1 ) ]
  
# the mapping from the reference triangle to the reference square. 
def coord_transform( c ):
    ( x , y ) = c
    if ( y != 1.0 ):
        a = 2.0 * (1.0 + x) / (1.0 - y) - 1.0
    else:
        a = -1.0
    return ( a , y )


