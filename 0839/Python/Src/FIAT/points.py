# points.py
# computes a lattice, plus gives lists of points
# for vertices, edges, and interior of lattice if requested
# precursor to Lagrange-type elements.
# all done with functional programming / list comprehension
# provides other services for manipulating points, such
# as taking the coordinates on [-1,1] and interpreting them
# as points along an edge of the reference triangle

from curry import curry

def make_lattice(n):
	return tuple( [ ( -1.0 + (2.0 * jj) / n  , -1.0 + (2.0 * ii) / n ) \
			for ii in xrange(0,n+1) \
			for jj in xrange(0,n-ii+1)] )

def make_vertex_point(v_no):
	if v_no == 0:
		return (-1.0,-1.0)
	elif v_no == 1:
		return (1.0,-1.0)
	elif v_no == 2:
		return (-1.0,1.0)
	else:
		raise RuntimeError, "Illegal argument to make_vertex_point"

# equispaced edge points interior to a particular edge.  Think of n as
# the order of polynomial defined using those points plus the
# vertices.  so then, n == 2 corresponds to only the edge midpoint, n
# == 3 gives two points in the interior
def make_edge_points(ed_no,n):
	if n == 1: return ()
	elif ed_no == 0:  # hypotenuse
		return tuple( [ ( -x , x ) \
				for x in edge_point_range(n) ] )
	elif ed_no == 1:  # vertical edge
		return tuple( [ ( -1.0 , -x ) \
				for x in edge_point_range(n) ] )
	elif ed_no == 2:
		return tuple( [ ( x , -1.0 ) \
				for x in edge_point_range(n) ] )
	else:
		raise RuntimeError, "Illegal argument to make_edge_points"

# given order n, returns the points strictly inside the triangle that
# lie on the lattice of order k.  Returns the empty list for n < 3
def make_interior_points(n):
	if n == 1 or n == 2: return ()
	else:
		return tuple( [ ( -1.0 + (2.0 * jj) / n  , -1.0 + (2.0 * ii) / n ) \
					for ii in range(1,n-1) \
					for jj in range(1,n-ii)] )

def interior_point_range(n):
	# hack to get around scoping issues in python
	f = curry( lambda m,i : -1.0 + ( 2.0 * i ) / m , n ) 
	return tuple( map( f , range(1,n) ) )

# useful for Gaussian quadrature along an edge, for example
def interval_to_edge_points(ed_no,pts):
	if ed_no == 0:   # hypotenuse
		return tuple( [ (-x,x) for x in pts ] )
	elif ed_no == 1: # vertical edge
		return tuple( [ (-1.0,-x) for x in pts ] )
	elif ed_no == 2: # horizontal edge
		return tuple( [ (x,-1.0) for x in pts ] )
	else:
		raise RuntimeError, "Illegal argument to interval_points_to_edge_points"
		
def edge_point_range(n):
	f = curry( lambda m,i : -1.0 + (2.0*i)/m , n )
        return tuple( map( f , range(1,n) ) )

