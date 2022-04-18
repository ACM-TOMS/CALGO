import dualbasis, functional, points, functionspace , operator
import func, quadrature
from func import gradient, curl, bubble, r
from math import sqrt

from Numeric import *

class VectorLagrangeDualBasis( dualbasis.DualBasis ):
    def __init__( self , n ):
        ls = reduce( operator.add , [ [ functional.ComponentPointEvaluation( i , points.make_vertex_point( v ) ) \
               for v in range( 0 , 3 ) ] \
           + [ functional.ComponentPointEvaluation( i , pt ) \
               for pt in points.make_edge_points( 0 , n ) ] \
           + [ functional.ComponentPointEvaluation( i , pt ) \
               for pt in points.make_edge_points( 1 , n ) ] \
           + [ functional.ComponentPointEvaluation( i , pt ) \
               for pt in points.make_edge_points( 2 , n ) ] \
           + [ functional.ComponentPointEvaluation( i , pt ) \
               for pt in points.make_interior_points( n ) ] for i in (0,1) ] )
	dimPn = (n+1)*(n+2)/2
        vertex_nodes = [ [ 0 , dimPn ] , \
			[ 1 , 1+dimPn ] , \
			[ 2 , 2+dimPn ] ]
        edge_nodes = [ range( 3 , 3 + ( n - 1 ) ) +\
			range( 3 + dimPn , dimPn + 3 + ( n - 1 ) ) , \
                       range( 3 + ( n - 1 ) , 3 + 2 * ( n - 1 ) ) + \
			range(3+(n-1)+dimPn,dimPn+3+2*(n-1)) , \
                       range( 3 + 2 * ( n - 1 ) , 3 + 3 * ( n - 1 ) ) + \
                       range(dimPn+3+2*( n - 1 ) ,dimPn+3+3 * ( n - 1 ) ) ]
        interior_nodes = range( 3+3*(n-1),(n+1)*(n+2)/2) + \
			range( dimPn+3+3*(n-1),dimPn+(n+1)*(n+2)/2) 
        dualbasis.DualBasis.__init__( self , ls , \
                                      vertex_nodes , edge_nodes , \
                                      interior_nodes )

class Lagrange( functionspace.ListTypeFiniteElementSpace ):
    def __init__( self , n ):
        U = functionspace.VectorPolynomialTriangle( n )
        UDual = VectorLagrangeDualBasis( n )
        functionspace.ListTypeFiniteElementSpace.__init__( self , U , UDual )

rt2 = sqrt(2.0)
rt2inv = 1.0 / rt2
normals = [ ( rt2inv , rt2inv ) , ( -1.0 , 0.0 ) , ( 0.0 , -1.0 ) ]

# makes the nodes corresponding to point evaluation at n+1 points on the interior
# of an edge.  This function can be used in all three H(div) elements
def make_edge_nodes( ed_no , n ):
	return [ functional.PointNormal( x , normals[ ed_no ] ) \
			for x in points.make_edge_points( ed_no , n + 2 ) ]


# interior nodes for BDM and BDFM are the same, so we write one function that we will
# use in both constructors
def make_BD_interior_nodes( n ):
	if n == 1:
		return []
	else:
		q = quadrature.JacobiQuadrature2D( 2 * n )
		P = functionspace.PolynomialTriangle( n - 1 )
		ps = P.getBasis()
		return [ functional.FunctionMoment( gradient( p ) , q ) \
				for p in ps[1:] ] \
			+ [ functional.FunctionMoment( curl( bubble * p ) , q ) \
				for p in ps[:((n-1)*n/2)] ]

# Raviart-Thomas interior nodes are moments against vector polynomials of degree n - 1
def make_RT_interior_nodes( n ):
	if n == 0:
		return []
	else:
		q = quadrature.JacobiQuadrature2D( 2 * n )
		V = functionspace.VectorPolynomialTriangle( n - 1 )
		vs = V.getBasis()
		return [ functional.FunctionMoment( v , q ) \
						for v in vs ] 


# dual spaces -- lists of nodes.  There are no vertex nodes, edge nodes are just point evaluation,
# and interior nodes are as given above. 
class BDMDualSpace( dualbasis.DualBasis ):
	def __init__( self , n ):
		ls = reduce( operator.add , [ make_edge_nodes( ed , n ) for ed in [ 0 , 1 , 2 ] ] ) \
		  + make_BD_interior_nodes( n )
		vertex_nodes = [ [] , [] , [] ]
		edge_nodes = [ range(0,n+1) , range( n+1,2*(n+1)) , range( 2*(n+1) , 3*(n+1) ) ]
		interior_nodes = range( 3*( n + 1) , (n+1)*(n+2) ) 
		dualbasis.DualBasis.__init__( self , ls , vertex_nodes , edge_nodes , interior_nodes )

class BDFMDualSpace( dualbasis.DualBasis ):
	def __init__( self , n ):
		ls = reduce( operator.add , [ make_edge_nodes( ed , n-1 ) for ed in [ 0 , 1 , 2 ] ] ) \
		  + make_BD_interior_nodes( n )
		vertex_nodes = [ [] , [] , [] ]
		edge_nodes = [ range(0,n) , range( n,2*(n)) , range( 2*(n) , 3*(n) ) ]
		interior_nodes = range( 3*( n) , (n+1)*(n+2)-3 ) 
		dualbasis.DualBasis.__init__( self , ls , vertex_nodes , edge_nodes , interior_nodes )

class RTDualSpace( dualbasis.DualBasis ):
	def __init__( self , n ):
		ls = reduce( operator.add , [ make_edge_nodes( ed , n ) for ed in [ 0 , 1 , 2 ] ] ) \
		  + make_RT_interior_nodes( n )
		vertex_nodes = [ [] , [] , [] ]
		edge_nodes = [ range(0,n+1) , range( n+1,2*(n+1)) , range( 2*(n+1) , 3*(n+1) ) ]
		interior_nodes = range( 3*( n + 1) , (n+1)*(n+3) ) 
		dualbasis.DualBasis.__init__( self , ls , vertex_nodes , edge_nodes , interior_nodes )
# function spaces for the H(div) elements.  We don't include BDM_k since it is trivially
# ( P_k )^2.  These two spaces serve as examples of spaces that are:
#	- constructed as the null space of a collection of linear functionals (BDFM)
#	- constructed by orthogonalizing a given set of functions (RT)

class BDFMFunctionSpace( functionspace.ListTypeConstrainedSpace ):
	def __init__( self , n ):
		line_rule = quadrature.JacobiQuadrature( 0 , 0 , 2 * n )
		Ls = [ functional.NormalEdgeMoment( ed , n , normals[ed] , line_rule ) \
				for ed in [ 0 , 1 , 2 ] ]
		U = functionspace.VectorPolynomialTriangle( n )
		functionspace.ListTypeConstrainedSpace.__init__( self , U , Ls )

class RTFunctionSpace( functionspace.ListTypeOrthonormalizedFunctionSpace ):
	def __init__( self , n ):
		q = quadrature.JacobiQuadrature2D( 2 * ( n + 1 ) )
		Vn = functionspace.VectorPolynomialTriangle( n )
		vns = Vn.getBasis()
		Vnp1 = functionspace.VectorPolynomialTriangle( n + 1 )
		P = functionspace.PolynomialTriangle( n )
		ps = P.getBasis()
		pnonlys = ps[(n*(n+1)/2):]
		us = vns + [ p * r for p in pnonlys ]
		functionspace.ListTypeOrthonormalizedFunctionSpace.__init__( self , Vnp1 , us , q )

# These three classes define the H(div) finite elements.
# They just instantiate the base function space, the dual space, and
# call the base class constructor for VectorFiniteElementSpace

class BDMElement( functionspace.ListTypeFiniteElementSpace ):
	def __init__( self , n ):
		V = functionspace.VectorPolynomialTriangle( n )
		Vdual = BDMDualSpace( n )
		functionspace.ListTypeFiniteElementSpace.__init__( self , V , Vdual )

class BDFMElement( functionspace.ListTypeFiniteElementSpace ):
	def __init__( self , n ):
		V = BDFMFunctionSpace( n )
		Vdual = BDFMDualSpace( n )
		functionspace.ListTypeFiniteElementSpace.__init__( self , V , Vdual )
		
class RTElement( functionspace.ListTypeFiniteElementSpace ):
	def __init__( self, n ):
		V = RTFunctionSpace( n )
		Vdual = RTDualSpace( n )
		functionspace.ListTypeFiniteElementSpace.__init__( self , V , Vdual )

        
