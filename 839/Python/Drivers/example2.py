# example2.py
# tabulates the Raviart-Thomas basis functions at quadrature points
#

# get the modules we need
from FIAT import quadrature, vectorelement, func
from FIAT.func import divergence

# set the degree of the space
degree = 2

# create a quadrature rule
integrate = quadrature.JacobiQuadrature2D( 2 * degree - 1 )

# create the finite element space
U = vectorelement.RTElement( degree )

# output in the format we talk about in the manual
xs = integrate.get_points()
ws = integrate.get_weights()
nqp = len( xs )

# print the quadrature rule
print nqp
for ( x , w ) in zip( xs , ws ):
	print x[0], x[1], w
print

# print the basis functions and their divergences
# We can just print the divergence since the Piola transform typically
# used to map the reference element to a physical element renders
# (div(psi),w) invariant.
us = U.getBasis()
print len( us )
for u in us:
	u_at_qps = map( u , xs )
	div_u_at_qps = map( divergence( u ) , xs )
	print
	for i in range( 0 , nqp ):
		print u_at_qps[i][0] , u_at_qps[i][1] , div_u_at_qps[i]
