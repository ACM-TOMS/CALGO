# Curries a function
# Given f( x1 , x2 , x3 ), for example,
# g = curry( f , a , b ) is a function of one argument, x3,
# so that g( c ) has the same value as f( a , b , c )
# This class is useful in bypassing the scoping bugs for higher
# order functions in Python
class curry:
    def __init__( self , f , *first_args ):
        self.f = f
        self.first_args = first_args[:]    # copy
    def __call__( self , *xs ):
        return self.f( *(self.first_args + xs) )

