# example1.py
# tabulates the Lagrange basis functions at a set of quadrature points
#

# get the modules we need
from FIAT import quadrature, scalarelement, func
from FIAT.func import D

# set the degree of the space
degree = 3

# create a quadrature rule
integrate = quadrature.JacobiQuadrature2D( 2 * degree - 1 )

# create the finite element space
U = scalarelement.Lagrange( degree )

# output in the format we talk about in the manual
xs = integrate.get_points()
ws = integrate.get_weights()
nqp = len( xs )

# print the quadrature rule
print nqp
for ( x , w ) in zip( xs , ws ):
	print x[0], x[1], w
print

# print the basis functions
us = U.getBasis()
print len( us )
for u in us:
	u_at_qps = map( u , xs )
	ux_at_qps = map( D(0)( u ) , xs )
	uy_at_qps = map( D(1)( u ) , xs )
	print
	for i in range( 0 , nqp ):
		print u_at_qps[ i ] , ux_at_qps[ i ] , uy_at_qps[ i ]
