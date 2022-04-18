# This module gives a smattering of classes implementing various
# linear functionals.

import points, quadrature, jacobi , operator
from Numeric import *

class PointEvaluation:
	def __init__( self , pt ):
		self.pt = pt
	def __call__( self , f ):
		return f( self.pt )
	def __str__( self ):
		return "Delta_{" + str( self.pt ) + "}"

class ComponentPointEvaluation:
	def __init__( self , comp , pt ):
		self.comp = comp
		self.pt = pt
	def __call__( self , v ):
		return v[self.comp]( self.pt )

class SymmetricTensorComponentPointEvaluation:
	def __init__( self , comp , pt ):
		self.comp = comp
		self.pt = pt
	def __call__( self , tau ):
		return tau.as[self.comp]( self.pt )

class FunctionMoment:
	def __init__( self , g , integrate ):
		self.g = g
		self.integrate = integrate
	def __call__( self , f ):
		return self.integrate( f * self.g )

class EdgeMoment:
        def __init__( self , ed_no , degree , line_rule ):
                self.edge_rule = quadrature.to_edge_quad( ed_no , line_rule )
                self.vals = array( [ jacobi.eval_jacobi( 0 , 0 , degree , x ) \
                                for x in line_rule.get_points( ) ] )
        def __call__( self , f ):
                fvals = array( f( self.edge_rule.get_points() ) )
		ws = self.edge_rule.get_weights()
                return sum( ws * fvals * self.vals )

class PointNormal:
        def __init__(self,x,n):
                self.x = x
                self.n = n
        def __call__(self,v):
                vx = v(self.x)
                return vx[0] * self.n[0] + vx[1] * self.n[1]

class ComponentPointNormal:
	def __init__( self, component , x , n ):
		self.component = component
		self.x = x
		self.n = n
	def __call__( self , w ):
		return reduce( operator.add , \
				[ c * n \
				for (c,n) in zip( w[self.component](self.x) , \
						  self.n ) ] )

class NormalEdgeMoment:
        def __init__( self , ed_no , degree , normal , line_rule ):
                self.edge_rule = quadrature.to_edge_quad( ed_no , line_rule )
                self.vals = array( [ jacobi.eval_jacobi( 0 , 0 , degree , x ) \
                                        for x in line_rule.get_points( ) ] )
                self.normal = normal                    
        def __call__( self , f ):
                fvals = array( [ (self.normal[0] * f[0] + self.normal[1] * f[1])( x ) for x in self.edge_rule.get_points() ] )
                return sum( self.edge_rule.get_weights() * fvals * self.vals )

