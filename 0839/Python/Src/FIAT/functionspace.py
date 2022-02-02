# functionspace.py

import types, operator, dubiner, points, func
from Numeric import *
import LinearAlgebra

pinv = LinearAlgebra.generalized_inverse
svd = LinearAlgebra.singular_value_decomposition

class FunctionSpace:
    def __init__( self , fs , dmats ):
        self.cache = {}
        self.fs = fs
        self.dmats = dmats
    def getBasis(self):
        return self.fs
    def getDmat( self , i ):
        return self.dmats[ i ]
    def tabulateBasis( self , xs ):
        "Array A_{i,j} = phi_j ( x_i )"
        if not self.cache.has_key( xs ):
            self.cache[ xs ] = [ f( xs ) for f in self.fs ] 
        return self.cache[ xs ]

class LinCombFunctionSpace( FunctionSpace ):
    def __init__( self , U , backward ):
        if not isinstance( U , LinCombFunctionSpace ):
            self.U = U
            self.backward = backward
            self.forward = pinv( backward )
        else:
            self.U = U.U
            self.backward = dot( U.backward , backward )
            self.forward = dot( pinv( backward ) , U.forward )
        fs = [ make_member( self.U , dofs ) \
               for dofs in transpose( self.backward ) ]
        dmats = [ dot( self.forward , dot( dmat , self.backward ) ) \
                  for dmat in self.U.dmats ]
        FunctionSpace.__init__( self , fs , dmats )

class ConstrainedFunctionSpace( LinCombFunctionSpace ):
    def __init__( self , U , ls ):
        a = array( [ [ l( u ) for u in U.getBasis() ] for l in ls ] )
        (U_a , Sig_a , Vt_a) = svd( a , 1 )
        vcal = Vt_a[ len( ls ): ]
        LinCombFunctionSpace.__init__( self , transpose( vcal ) )

class OrthonormalizedFunctionSpace( LinCombFunctionSpace ):
    def __init__( self , U , vs , integrate ):
        a = array( [ [ integrate( u * v ) \
                       for v in vs ] for u in U.getBasis() ] )
        (U_a , Sig_a , Vt_a) = svd( a )
        LinCombFunctionSpace.__init__( self , U , U_a )

class FiniteElementSpace( LinCombFunctionSpace ):
    def __init__( self , U , dual ):
        v = array( [ [ node( u ) \
                       for u in U.getBasis() ] \
                     for node in dual.getBasis() ] )
        vinv = LinearAlgebra.inverse( v )
        self.dual = dual
        LinCombFunctionSpace.__init__( self , U , vinv )

class ListTypeFunctionSpace:
    def __init__( self , vs ):
        self.cache = {}
        self.vs = vs
    def getBasis( self ):
        return self.vs
    def tabulateBasis( self , xs ):
        if not self.cache.has_key( xs ):
            self.cache[ xs ] = [ v( xs ) for v in self.vs ] 
        return self.cache[ xs ]

class ListTypeLinCombFunctionSpace( ListTypeFunctionSpace ):
    def __init__( self , U , backward ):
        if not isinstance( U , ListTypeLinCombFunctionSpace ):
            self.U = U
            self.backward = backward
        else:
            self.U = U.U
            self.backward = dot( U.backward , backward )	
        bs = [ make_member( self.U , dof ) \
               for dof in transpose( self.backward ) ]
        ListTypeFunctionSpace.__init__( self , bs )

class ListTypeConstrainedSpace( ListTypeLinCombFunctionSpace ):
     def __init__( self , V , Ls ):
         vs = V.getBasis()
         Lmat = array( [ [ l(v) for v in vs ] for l in Ls ] )
         (U_L,Sig_L,Vt_L) = svd(Lmat,1)
         Vcal = Vt_L[ len(Ls): ]
         ListTypeLinCombFunctionSpace.__init__( self , V , transpose( Vcal ) )

class ListTypeOrthonormalizedFunctionSpace( ListTypeLinCombFunctionSpace ):
    def __init__( self, U , vs , q ):
        a = array( [ [ q.integrate( u * v ) \
                       for v in vs ]\
                     for u in U.getBasis() ] )
        (U_a , Sig_a , Vt_a) = svd( a )
        ListTypeLinCombFunctionSpace.__init__( self , U , U_a )

class ListTypeFiniteElementSpace( ListTypeLinCombFunctionSpace ):
    def __init__( self , U , dual ):
	us = U.getBasis()
	nodes = dual.getBasis()
        v = array( [ [ node( u ) \
                       for u in us ] \
                     for node in nodes ] )
        vinv = LinearAlgebra.inverse( v )
        self.dual = dual
        ListTypeLinCombFunctionSpace.__init__( self , U , vinv )

class PolynomialTriangle( FunctionSpace ):
    def __init__( self, n ):
        bs = dubiner.make_dubiner_funcs( n )
        if n == 0:
            dmats = [ array( [ [ 0.0 ] ] , "d" ) ] * 2
        else:       
            pts = points.make_lattice( n )
            v = array( [ [ b( x ) for b in bs ] for x in pts ] )
            vinv = LinearAlgebra.inverse( v )
            dtildes = [ transpose( [ map( func.D(i)( b ) , pts ) \
                                     for b in bs ] ) \
                        for i in (0,1) ]
            dmats = [ dot( vinv , dtilde ) for dtilde in dtildes ]
        FunctionSpace.__init__( self , bs , dmats )

class VectorPolynomialTriangle( ListTypeFunctionSpace ):
    def __init__( self , n ):
        ps = dubiner.make_dubiner_funcs( n )
        z = func.zero
        vs = [ func.MemoizedVectorFunction( [ p , z ] ) \
               for p in ps ] \
               + [ func.MemoizedVectorFunction( [ z , p ] ) \
                   for p in ps ]
        ListTypeFunctionSpace.__init__( self , vs )

class TensorPolynomialTriangle( ListTypeFunctionSpace ):
    def __init__( self , n ):
        ps = dubiner.make_dubiner_funcs( n )
        z = func.zero
        ts = [ func.MemoizedTensorFunction( \
            [ [ p , z ] , [ z , z ] ] ) \
               for p in ps ] \
           + [ func.MemoizedTensorFunction( \
            [ [ z , p ] , [ z , z ] ] ) \
               for p in ps ] \
           + [ func.MemoizedTensorFunction( \
            [ [ z , z ] , [ p , z ] ] ) \
               for p in ps ] \
           + [ func.MemoizedTensorFunction( \
            [ [ z , z ] , [ z , p ] ] ) \
               for p in ps ]
        ListTypeFunctionSpace.__init__( self , ts)

class SymmetricTensorPolynomialTriangle( ListTypeFunctionSpace ):
    def __init__( self , n ):
	ps = dubiner.make_dubiner_funcs( n )
	z = func.zero
	ts = [ func.MemoizedSymmetricTensorFunction( p , z , z ) \
		for p in ps ] \
	   + [ func.MemoizedSymmetricTensorFunction( z , p , z ) \
		for p in ps ] \
	   + [ func.MemoizedSymmetricTensorFunction( z , z , p ) \
		for p in ps ]
	ListTypeFunctionSpace.__init__( self , ts )

def make_member( U , dof ):
	if isinstance( U , FunctionSpace ):
		return func.MemoizedFunction( FunctionSpaceMember( U , dof ) )
	else:
		return reduce( operator.add , \
			[ a * u \
			for ( a , u ) in zip( dof , U.getBasis() ) ] )

class FunctionSpaceMember:
	def __init__( self, U , dof ):
		self.U = U
		self.dof = dof
	def __call__( self , pts ):
		us = self.U.tabulateBasis( pts )
		return reduce( operator.add , \
			[ a * u for ( a, u ) in zip( self.dof , us ) ] )
	def deriv( self, i ):
		return FunctionSpaceMember( \
			self.U , dot( self.U.getDmat( i ) , self.dof ) )