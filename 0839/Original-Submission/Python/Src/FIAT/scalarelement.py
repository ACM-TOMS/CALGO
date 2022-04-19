import dualbasis, functional, points, functionspace, quadrature, jacobi, math
from Numeric import *
from func import lam0,lam1,lam2


class LagrangeDualBasis( dualbasis.DualBasis ):
    def __init__( self , n ):
        ls = [ functional.PointEvaluation( points.make_vertex_point( v ) ) \
               for v in range( 0 , 3 ) ] \
           + [ functional.PointEvaluation( pt ) \
               for pt in points.make_edge_points( 0 , n ) ] \
           + [ functional.PointEvaluation( pt ) \
               for pt in points.make_edge_points( 1 , n ) ] \
           + [ functional.PointEvaluation( pt ) \
               for pt in points.make_edge_points( 2 , n ) ] \
           + [ functional.PointEvaluation( pt ) \
               for pt in points.make_interior_points( n ) ]
        vertex_nodes = [ [ 0 ] , [ 1 ] , [ 2 ] ]
        edge_nodes = [ range( 3 , 3 + ( n - 1 ) ) , \
                       range( 3 + ( n - 1 ) , 3 + 2 * ( n - 1 ) ) , \
                       range( 3 + 2 * ( n - 1 ) , 3 + 3 * ( n - 1 ) ) ]
        interior_nodes = range( 3 + 3 * ( n - 1 ) , ( n + 1 ) * ( n + 2 ) / 2 )
        dualbasis.DualBasis.__init__( self , ls , \
                                      vertex_nodes , edge_nodes , \
                                      interior_nodes )

class NonconformingDualBasis( dualbasis.DualBasis ):
    def __init__( self , n ):
        gjs = jacobi.compute_gauss_jacobi_points( 0 , 0 , n  )
        ls = [ functional.PointEvaluation( pt ) \
               for pt in points.interval_to_edge_points( 0 , gjs ) ] \
           + [ functional.PointEvaluation( pt ) \
               for pt in points.interval_to_edge_points( 1 , gjs ) ] \
           + [ functional.PointEvaluation( pt ) \
               for pt in points.interval_to_edge_points( 2 , gjs ) ] \
           + [ functional.PointEvaluation( pt ) \
               for pt in points.make_interior_points( n ) ]
        if n % 2 == 0:
            q = quadrature.JacobiQuadrature2D( 3 * n )
            ls.append( lambda f : q.integrate ( f ) )
        vertex_nodes = [ [] , [] , [] ]
        edge_nodes = [ range( 0 , n ) , \
                       range( n , 2*n ) , \
                       range( 2*n , 3*n ) ]
        interior_nodes = range( (3*n) , len( ls ) )
        dualbasis.DualBasis.__init__( self , ls ,\
                                      vertex_nodes , \
                                      edge_nodes , \
                                      interior_nodes ) 



def make_nonconforming_function_space( n ):
	Pn = functionspace.PolynomialTriangle( n )
	if n % 2 != 0:
		return Pn
	else:
		vs = Pn.getBasis() + [ nc_extra( n ) ]
		Pnp1 = functionspace.PolynomialTriangle( n + 1 )	
		q = quadrature.JacobiQuadrature2D( 2 * n + 3 )
		return functionspace.OrthonormalizedFunctionSpace( Pnp1 , vs , q )

class Lagrange( functionspace.FiniteElementSpace ):
    def __init__( self , n ):
        U = functionspace.PolynomialTriangle( n )
        UDual = LagrangeDualBasis( n )
        functionspace.FiniteElementSpace.__init__( self , U , UDual )
        
class Nonconforming( functionspace.FiniteElementSpace ):
    def __init__( self , n ):
        U = make_nonconforming_function_space( n )
        Udual = NonconformingDualBasis( n )
        functionspace.FiniteElementSpace.__init__( self , U , Udual )

def nc_extra( n ):
    return ( lam1 - lam0 ) * ( lam2 - lam1 ) * ( lam0 - lam2 )  \
           * ( ( lam0 * lam1 ) ** ( ( n - 2 ) / 2 )  \
               + ( lam1 * lam2 ) ** ( ( n - 2 ) / 2 ) \
               + ( lam0 * lam2 ) ** ( ( n - 2 ) / 2 ) )
