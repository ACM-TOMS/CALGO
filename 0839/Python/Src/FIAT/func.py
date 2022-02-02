# func.py
# defines memoized scalar, vector, and tensor-valued functions,
# plus differentiation and algebraic operations.

import operator
from Numeric import *
from curry import curry

class MemoizedFunction:
    """MemoizedFunction wraps functions that map a tuple to a real value
       with a memoizing layer.
       The constructor takes the point --> real map
       The model is that function evaluation is expensive, so we
       memoized the results of call
       MemoizedFunctions support some operator overloading --
          scalar multiplication( __rmul__ )
          multiplication onto functions, vectors, or scalars
          addition
          negation
    """
    def __init__( self , f ):
        self.f = f
        self.cache = {}
	self.deriv_cache = {}
    def __call__( self , x ):
        if not self.cache.has_key( x ):
            self.cache[ x ] = self.f( x )
        return self.cache[ x ]
    def __str__( self ):
        return "Memoized " + str( self.f )
    def deriv( self , ii ):
	# deriv remembers when it's been called before and returns
	# the same instance.  This way, all differentiation of the
	# same object returns the same instance.
	if not self.deriv_cache.has_key( ii ):
	    self.deriv_cache[ ii ] = MemoizedFunction( self.f.deriv(ii) )
        return self.deriv_cache[ ii ]
    def __mul__( self , other ):
        if isinstance( other , MemoizedFunction ):
            return MemoizedFunction( Product( self , other ) )
        elif isinstance( other , MemoizedVectorFunction ):
            return MemoizedVectorFunction( [ self * u for u in other ] )
        elif isinstance( other , MemoizedTensorFunction ):
            return MemoizedTensorFunction( [ [ self * u for u in other_rows ] \
                                             for other_rows in other ] )
	elif isinstance( other , MemoizedSymmetricTensorFunction ):
	    return MemoizedSymmetricTensorFunction( \
		self * other.a11 , self * other.a12 , self * other.a22 )
        else:
            raise RuntimeError, "Illegal argument to MemoizedFunction.__mul__"
    def __rmul__( self , other ):
        return MemoizedFunction( Scaling( other , self ) )
    def __add__( self , other ):
        if not isinstance( other , MemoizedFunction ):
            raise RuntimeError, "Illegal operand to MemoizedFunction.__add__"
        else:
            return MemoizedFunction( Sum( self, other ) )
    def __sub__( self , other ):
        if not isinstance( other , MemoizedFunction ):
            raise RuntimeError, "Illegal operand to MemoizedFunction.__add__"
        else:
            return MemoizedFunction( Diff( self, other ) )       
    def __neg__( self ):
        return MemoizedFunction( Scaling( -1.0 , self ) )
    def __pow__( self , n ):
        return MemoizedFunction( Power( self , n ) )

class Zero( MemoizedFunction ):
    def __init__( self ):
        self.f = None
    def __call__( self , x ):
        return 0.0
    def deriv( self , dir ):
        return self
    def __str__( self ):
	return "0(func)"

zero = Zero()

class MemoizedVectorFunction:
    """Converts a list of MemoizedFunctions into a memoized
vector-valued function"""
    def __init__( self , components ):
        self.components = components
        self.cache = {}
    def __call__( self , x ):
        if not self.cache.has_key( x ):
            self.cache[ x ] = array( [ u( x ) for u in self ] )
        return self.cache[ x ]
    def __str__( self ):
        return str( map( str , self ) )
    def __len__( self ):
        return len( self.components )
    def __mul__( self , other ):
	# vector * vector ==> dot product
        if isinstance( other , MemoizedVectorFunction ):
            return reduce( operator.add , [ u * v for ( u , v ) in zip( self , other ) ] )
        # vector * function --> scaled vector
        elif isinstance( other, MemoizedFunction ):
            return MemoizedVectorFunction( [ u * other for u in self ] )
        # vector * tensor --> vector
        elif isinstance( other , MemoizedTensorFunction ):
            return MemoizedVectorFunction( \
                [ reduce( operator.add , \
                          [ u * v for ( u , v ) in zip( self , othercol ) ] ) \
                  for othercol in transpose_tensor( other ) ] )
        else:
            raise RuntimeError, "multiplication error in MemoizedVectorFunction"
    def __rmul__( self , other ):
        return MemoizedVectorFunction( [ other * u for u in self.components ] )
    def __add__( self , other ):
        return MemoizedVectorFunction( [ u + v for ( u , v ) in zip( self , other ) ] ) 
    def __neg__( self ):
        return MemoizedVectorFunction( [ -u for u in self ] )
    def __getitem__( self , i ):
        return self.components[ i ]

class MemoizedTensorFunction:
    """Converts list of lists of MemoizedFunction or a list of
    MemoizedVectorFunction into a tensor"""
    # components stored as a list of vectors
    def __init__( self , components ):
        self.cache = {}
        if isinstance( components[0] , MemoizedVectorFunction ):
            self.components = components
        elif isinstance( components[0][0] , MemoizedFunction ):
            self.components = [ MemoizedVectorFunction( u ) \
                                for u in components ]
    def __call__( self , x ):
        if not self.cache.has_key( x ):
            self.cache[ x ] = array( [ u( x ) for u in self ] )
        return self.cache[ x ]
    def __len__( self ):
        return len( self.components )
    def __mul__( self , other ):
        #Dyadic product - sum over componentwise product
        if isinstance( other , MemoizedTensorFunction ):
            return reduce( operator.add , \
                           [ reduce( operator.add , \
                                     [ u * v \
                                       for ( u , v ) in zip( selfrow, otherrow ) ] ) \
                             for ( selfrow , otherrow ) in zip( self, other ) ] )
        # tensor-vector product
        elif isinstance( other , MemoizedVectorFunction ):
            return MemoizedVectorFunction( \
                [ reduce( operator.add , \
                          [ u * v \
                            for (u,v) in zip( selfrow , other ) ] ) \
                          for selfrow in self ] )
        # scale each component of the tensor
        elif isinstance( other , MemoizedFunction ):
            return MemoizedTensorFunction( \
                [ u * other for u in self ] )
    def __rmul__( self , other ):
        return MemoizedTensorFunction( [ other * u for u in self ] )
    def __add__( self, other ):
        return MemoizedTensorFunction( \
            [ u + v for ( u , v ) in zip( self , other ) ] )
    def __neg__( self ):
        return MemoizedTensorFunction( \
            [ -u for u in self ] )
    def __getitem__( self , i ):
        return self.components[ i ]

class MemoizedSymmetricTensorFunction:
    """Only works for 2x2 tensors -- you put in the three unique components"""
    def __init__( self , a11 , a12 , a22 ):
        self.cache = {}
	self.as = [ a11 , a12 , a22 ]
	self.components = \
		[ MemoizedVectorFunction( [ self.as[0] , self.as[1] ] ) , \
		  MemoizedVectorFunction( [ self.as[1] , self.as[2] ] ) ]
    def __call__( self , x ):
        if not self.cache.has_key( x ):
	    asatx = [ a( x ) for a in self.as ]
            self.cache[ x ] \
		= array( [ [ asatx[0] , asatx[1] ] , \
			   [ asatx[1] , asatx[2] ] ] )
        return self.cache[ x ]
    def __len__( self ):
        return len( self.components )
    def __mul__( self , other ):
	# symmetric tensor * symmetric tensor ==> dyadic product
	if isinstance( other , MemoizedSymmetricTensorFunction ):
	    return self.as[0] * other.as[0] \
		+ 2 * self.as[1] * other.as[1] \
		+ self.as[2] * other.as[2]
	elif isinstance( other , MemoizedTensorFunction ):
	    if len( other ) != 2:
		raise RuntimeError, "Can't multiply in MSTF"
	    return self.a11 * other[0][0] \
		+ self.a12 * ( other[1][0] + other[0][1] ) \
		+ self.a22 * other[1][1]
	elif isinstance( other , MemoizedVectorFunction ):
	    if len( other ) != 2:
		raise RuntimeError, "Can't multiply in MSTF"
	    return MemoizedVectorFunction( \
			[ self.as[0] * other[0] + self.as[1] * other[1] ,\
			  self.as[1] * other[0] + self.as[2] * other[1] ] )
	elif isinstance( other , MemoizedFunction ):
	    return MemoizedSymmetricTensorFunction( \
		self.as[0] * other , self.as[1] * other , self.as[2] * other )
    def __rmul__( self , other ):
        return MemoizedSymmetricTensorFunction(  \
		other * self.as[0] , other * self.as[1] , other * self.as[2] )
    def __add__( self, other ):
	if isinstance( other , MemoizedSymmetricTensorFunction ):
	    return MemoizedSymmetricTensorFunction( \
		self.as[0] + other.as[0] , self.as[1] + other.as[1] , \
		self.as[2] + other.as[2] )
	elif isinstance( other , MemoizedTensorFunction ):
	    if len( other ) != 2:
		raise RuntimeError, "Can't add in MSTF"
	    return MemoizedTensorFunction( \
		[ u + v for (u,v) in zip( self , other ) ] )
    def __neg__( self ):
        return MemoizedTensorFunction( \
	    -self.as[0] , -self.as[1] , - self.as[2] )
    def __getitem__( self , i ):
	return self.components[ i ]

class Sum:
    """Models f + g for MemoizedFunctions f and g."""
    def __init__( self , f , g ):
        self.f = f
        self.g = g
    def __call__( self , x ):
        return self.f( x ) + self.g( x )
    def str( self ):
        return " ( " + str( self.f ) + " + " + str( self.g ) + " ) "
    def deriv( self , ii ):
	return self.f.deriv( ii ) + self.g.deriv( ii )

class Diff:
    """Models f - g for MemoizedFunctions f and g."""
    def __init__( self , f , g ):
        self.f = f
        self.g = g
    def __call__( self , x ):
        return self.f( x ) - self.g( x )
    def str( self ):
        return " ( " + str( self.f ) + " - " + str( self.g ) + " ) "
    def deriv( self , ii ):
        return self.f.deriv( ii ) - self.g.deriv( ii )

class Product:
    """Models f * g for MemoizedFunctions f and g."""
    def __init__( self , f , g ):
	self.f = f
	self.g = g
    def __call__( self , x ):
        return self.f( x ) * self.g( x )
    def str( self ):
        return " ( " + str( self.f ) + " * " + str( self.g ) + " ) "
    def deriv( self , ii ):
        return self.f.deriv( ii ) * self.g + self.f * self.g.deriv( ii )

class Scaling:
    """Models alpha * f for numeric alpha and MemoizedFunction f."""
    def __init__( self , alpha , f ):
        self.alpha = alpha
        self.f = f
    def __call__( self , x ):
        return self.alpha * self.f( x )
    def str( self ):
        return " ( " + str( self.alpha ) + " * " + str( self.f ) + " ) "
    def deriv( self , ii ):
        return self.alpha * self.f.deriv( ii )

class Power:
    """Models f ** n for MemoizedFunction f and numeric n."""
    def __init__( self , f, n ):
        self.f = f
        self.n = n
    def __call__( self , x ):
        return ( self.f( x ) ) ** self.n

def derivative( i , f ):
    return f.deriv( i )

def D( i ):
    """Given i, returns the function that differentiates an object
in the ith coordinate direction"""
    # since scoping for lambda functions is messed up in Python, I use
    # curry
    return curry( derivative , i )

def gradient( f ):
    """Computes vector-valued gradient of scalar functions, and
tensor-valued gradient of vector functions"""
    if isinstance( f , MemoizedFunction ):
        return MemoizedVectorFunction( [ D(0)( f ) , D(1)( f ) ] )
    elif isinstance( f , MemoizedVectorFunction ):
        return MemoizedTensorFunction( [ gradient( u ) for u in f ] )
    else:
        raise RuntimeError, "Gradient not implemented for that type"

def divergence( f ):
    """Computes scalar-valued divergence of vector-valued functions,
and vector-valued divergence of tensor functions"""
    if isinstance( f , MemoizedVectorFunction ):
        return reduce( operator.add , \
                       [ f[i].deriv( i ) \
                         for i in range(0,len( f )) ] )
    elif isinstance( f , MemoizedTensorFunction ) \
		or isinstance( f , MemoizedSymmetricTensorFunction ):
        return MemoizedVectorFunction( \
            [ divergence( v ) for v in f ]  )
    else:
        raise RuntimeError, "Divergence not implemented for that type"

def curl( f ):
    return MemoizedVectorFunction( [ D(1)( f ) , -D(0)( f ) ] )

def epsilon( f ):
	"""Symmetric part of the gradient of a vector-valued function"""
	if not isinstance( f , MemoizedVectorFunction ):
		raise RuntimeError, "Can't do that -- epsilon"
	return MemoizedSymmetricTensorFunction( \
		D(0)( f[0] ) , 0.5 * ( D(0)( f[1] ) + D(1)( f[0] ) ), \
		D(1)( f[1] ) )

# Operator mapping potentials to Airy stress
def J( f ):
	"""Airy stress operator, maps scalars to symmetric tensors"""
	d00 = lambda f: D(0)( D(0)( f ) )
	d11 = lambda f: D(1)( D(1)( f ) )
	d01 = lambda f: D(0)( D(1)( f ) )
	return MemoizedSymmetricTensorFunction( \
		d11( f ) , -d01( f ) , d00( f ) )

def transpose_tensor( A ):
    if not isinstance( A , MemoizedTensorFunction ):
        raise RuntimeError, "I barfed -- transpose_tensor"
    else:
        return MemoizedTensorFunction( apply( zip, tuple( A ) ) )

# these functions model the barycentric coordinates.
lam0base = lambda x: -0.5 * ( x[0] + x[1] )
def lam0deriv( i ):
	return lambda x : -0.5
lam0base.deriv = lam0deriv
lam0 = MemoizedFunction( lam0base )

lam1base = lambda x: 0.5 * ( 1. + x[0] )
def lam1deriv( i ):
	if i == 0:
		return lambda x : 0.5
	elif i == 1:
		return lambda x : 0.0
	else:
		raise RuntimeError, "illegal input -- lam1deriv"
lam1base.deriv = lam1deriv
lam1 = MemoizedFunction( lam1base )

lam2base = lambda x: 0.5 * ( 1. + x[1] )
def lam2deriv( i ):
	if i == 0:
		return lambda x : 0.0
	elif i == 1:
		return lambda x : 0.5
	else:
		raise RuntimeError, "illegal input -- lam1deriv"
lam2base.deriv = lam2deriv
lam2 = MemoizedFunction( lam2base )

# cubic bubble function, vanishes on the boundary
bubble = lam0 * lam1 * lam2

# the two basic monomials
xbase = lambda x: x[0]
def xbasederiv( i ):
	if i == 0: return lambda x: 1.0
	else: return lambda x: 0.0
xbase.deriv = xbasederiv
x = MemoizedFunction( xbase )

ybase = lambda x: x[1]
def ybasederiv( i ):
	if i == 0: return lambda x: 0.0
	else: return lambda x: 1.0
ybase.deriv = ybasederiv
y = MemoizedFunction( ybase )

# the position vector in R2, consisting of monomials.
r = MemoizedVectorFunction( [ x , y ] )




