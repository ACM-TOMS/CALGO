import functionspace, quadrature, dualbasis, points, functional, jacobi
from func import divergence, curl, D, J, bubble, epsilon
from Numeric import *

from math import sqrt

rt2 = sqrt(2.0)
rt2inv = 1.0 / rt2
normals = [ ( rt2inv , rt2inv ) , ( -1.0 , 0.0 ) , ( 0.0 , -1.0 ) ]


# Tensors Sigma^k = { tau P_{k+2}(S) : div tau in P_{k}(R^2)
class ArnoldWintherSigmaSpace( functionspace.ListTypeConstrainedSpace ):
	def __init__( self , k ):
		if k < 1:
			raise RuntimeError, "Illegal parameter: ArnoldWinther"
		# base function space
		U = functionspace.SymmetricTensorPolynomialTriangle( k + 2 )

		# but divergences of each component must only be in
		# P_k, so I need to get homogeneous polynomials
		# of degree k+1 to specify my constraints

		P = functionspace.PolynomialTriangle( k + 1 )
		ps = P.getBasis()
		a = (k+1)*(k+2)/2
		pkp1_homog = ps[a:]

		# I need a quadrature rule to integrate:
		# div( P_{k+2 ) against P_k+1 -- order 2*(k+1) 
		integrate = quadrature.JacobiQuadrature2D( 2 * ( k + 1 ) )

		ls = [ lambda( f ): integrate( p * divergence( f[0] )) \
			for p in pkp1_homog ] \
		   + [ lambda( f ): integrate( p * divergence( f[1] )) \
			for p in pkp1_homog ]

		# and now initialize the base class -- creates the basis
		functionspace.ListTypeConstrainedSpace.__init__( self , U , ls )

def make_ArnoldWinther_Mks( k ):
	if k < 2:
		return []
	else:
		Pkm2 = functionspace.PolynomialTriangle( k - 2 )
		ps = Pkm2.getBasis()
		b2 = bubble * bubble
		Pkp4 = functionspace.PolynomialTriangle( k + 4 )
		mkpotentials = [ b2 * p for p in ps ]
		q = quadrature.JacobiQuadrature2D( 2 * ( k + 4 ) )
		MkpotentialSpace = \
			functionspace.OrthonormalizedFunctionSpace( \
				Pkp4 , mkpotentials , q ) 
		return [ J( m ) for m in MkpotentialSpace.getBasis() ]

def make_epsilon_pks( k ):
	if k < 1:
		return []
	else:
		PkVec = functionspace.VectorPolynomialTriangle( k )
		pks = PkVec.getBasis()
		# throw away ones that are constant, plus a linear
		dimpk = (k+1)*(k+2)/2
		pkvs_touse = pks[1:dimpk] + pks[dimpk+2:2*dimpk]
		# return symmetric parts of gradients
		return map( epsilon , pkvs_touse )

def make_ArnoldWinther_Nks( k ):
	return make_ArnoldWinther_Mks( k ) + make_epsilon_pks( k )
	

class ArnoldWintherSigmaDualBasis( dualbasis.DualBasis ):
	def __init__( self , k ):
		gjs = jacobi.compute_gauss_jacobi_points( 0 , 0 , k + 1 )

		integrate = quadrature.JacobiQuadrature2D( 2 * k + 3 )

		vertex_ls = \
	[ functional.SymmetricTensorComponentPointEvaluation( i , x ) \
		for x in [ points.make_vertex_point( j ) \
				for j in [0,1,2] ] \
		for i in [0,1,2] ]

		edge_ls = \
	[ functional.ComponentPointNormal( i , x , normals[ed] ) \
		for ed in [ 0 , 1 , 2 ] \
#		for x in points.interval_to_edge_points( ed , gjs ) \
		for x in points.make_edge_points( ed , k + 2 ) \
		for i in [ 0 , 1 ] ]

		interior_ls = \
	[ functional.FunctionMoment( tau , integrate ) \
		for tau in make_ArnoldWinther_Nks( k ) ]

		ls = vertex_ls + edge_ls + interior_ls

		vertex_nodes = [ [ 0 , 1 , 2 ] , \
				 [ 3 , 4 , 5 ] , \
				 [ 6 , 7 , 8 ] ]
		edge_nodes = [ range( 9 , 9 + 2 * ( k + 1 ) ) , \
			       range( 9 + 2*( k + 1 ) , 9 + 4*( k + 1 ) ) , \
			       range( 9+4*(k+1) , 9+6*(k+1) ) ]

		interior_nodes = range( 9+6*(k+1) , len(ls) )

		dualbasis.DualBasis.__init__( self , ls , \
						vertex_nodes , \
						edge_nodes , \
						interior_nodes )

class ArnoldWintherStressElement( functionspace.ListTypeFiniteElementSpace ):
	def __init__( self , k ):
		functionspace.ListTypeFiniteElementSpace.__init__( self , \
			ArnoldWintherSigmaSpace( k ) , \
			ArnoldWintherSigmaDualBasis( k ) )
