/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  P2MESH : DRIVER FOR TESTING THE PACKAGE ON USER DEFINED MARKERS         |
 |                                                                          |
 |  date         : 1999, 11 October                                         |
 |  version      : 1.0                                                      |
 |  file         : testall_marker.cc                                        |
 |  authors      : Enrico Bertolazzi (1) & Gianmarco Manzini (2)            |
 |  affiliations :                                                          |
 |                                                                          |
 |  (1) Department of Mechanics and Structures Engineering                  |
 |      University of Trento                                                |
 |      via Mesiano 77, I -- 38050 Trento, Italy                            |
 |      email : enrico.bertolazzi@ing.unitn.it                              |
 |                                                                          |
 |  (2) Institute of Numerical Analysis -- CNR                              |
 |      via Ferrata 1, I -- 27100 Pavia, Italy                              |
 |      email: gianmarco.manzini@ian.pv.cnr.it                              |
 |                                                                          |
 |  purpose:                                                                |
 |    test P2MESH package on user defined markers.                          |
 |                                                                          |
 |  output file:                                                            |
 |    testall_marker.out                                                    |
 |                                                                          |
 |  documentation:                                                          |
 |    userman.ps -- user manual for P2MESH programmers                      |
 |                                                                          |
\*--------------------------------------------------------------------------*/

# ifndef FILE_MMESH
# define FILE_MMESH "meshes/mmesh"
# endif

# include "p2mesh.hh"
# include <string>
# include <math.h>

// declare the name of the user defined classes
class Vertex ;
class Edge   ;
class Quad   ;
class Mesh   ;

struct Marker {
  string name ;
  double value ;
} ;

inline
ostream &
operator << (ostream & s, Marker const & m)
{ s << "( " << m.name << " = " << m.value << " )" ; return s ; }

inline
istream &
operator >> (istream & s, Marker & m)
{ s >> m.name >> m.value ; return s ; }

// setup common code & data
class Common : public p2_common<Vertex,Edge,Quad,Mesh,
                                4,true,
                                double, int, unsigned,
                                Marker, Marker, Marker>
{ } ; 

// define class Vertex
class Vertex : public p2_vertex<Common> {
  typedef Common::Vmark Vmark ;
  Vmark v_marker ;
public:
  static void Set_BC(Vertex & V, Vmark const & m)
    { V.v_marker = m ; }
  Vmark const & marker(void) const { return v_marker ; }
} ;

// define class Edge
class Edge : public p2_edge<Common> {
  typedef Common::Emark Emark ;
  Emark e_marker ;
public:
  static void Set_BC(Edge & E, Emark const & m)
    { E.e_marker = m ; }
  Emark const & marker(void) const { return e_marker ; }
} ;

// define class Quad
class Quad : public p2_poly<Common> {
  typedef Common::Pmark Pmark ;
  Pmark q_marker ;
public:
  static void Set_BC(Quad & Q, Pmark const & m)
    { Q.q_marker = m ; }
  Pmark const & marker(void) const { return q_marker ; }
} ;

// define class Mesh
class Mesh : public p2_mesh<Common> {
	     
  typedef Common::Unsigned Unsigned ;
  typedef Common::Integer  Integer ;
  typedef Common::Real     Real ;

public:
  void print_marker(ostream & outs) const ;
} ;

static
void
test_build_mesh(ostream & file_out) {

  typedef Common::Unsigned Unsigned ;
  typedef Common::Integer  Integer ;
  typedef Common::Real     Real ;
  typedef Common::Vmark    Vmark ;
  typedef Common::Emark    Emark ;
  typedef Common::Pmark    Pmark ;

  cout << "testing build_mesh..." ; cout.flush() ;
  
  Unsigned const nv = 7 ;
  Real     const XY[2*nv] = { 0,0,  1,0,  2,0,  2,2,  0,2,  0,1,  1,1  } ;
  Vmark    mv[nv] ;
  for ( Unsigned i = 0 ; i < nv ; ++i )
    { mv[i].name = "Vertex" ; mv[i].value = i ; }
  
  Unsigned const ne = 9 ;  
  Unsigned const E[2*ne]  = { 0,1,  1,2,  2,3,  3,4,  4,5,  5,0,  1,6,  6,5,  6,3 } ;    
  Emark    me[ne] ;
  for ( Unsigned i = 0 ; i < ne ; ++i )
    { me[i].name = "Edge" ; me[i].value = i ; }
   
  Unsigned const np = 3 ;  
  Unsigned const P[4*np]  = { 0,1,6,5, 1,2,3,6,  5,6,3,4 } ;
  Pmark    mp[np] ; 
  for ( Unsigned i = 0 ; i < np ; ++i )
    { mp[i].name = "Poly" ; mp[i].value = i ; }

  Mesh mesh ; // instantiate an "empty" mesh object

  // method 82
  mesh.build_mesh( nv, XY, mv, Vertex :: Set_BC, // vertex arrays
		   ne,  E, me, Edge   :: Set_BC, // edge   arrays
		   np,  P, mp, Quad   :: Set_BC, // quad   arrays
		   0 ) ;                         // C MESH offset

  file_out << "** build_mesh **" << endl ;
  mesh.print_marker( file_out );
  cout << "done" << endl ;

}

static
void
test_read_mesh(ostream & file_out) {
  cout << "testing read_mesh..." ; cout.flush() ;
  Mesh mesh ; // instantiate an "empty" mesh object
  mesh.read_mesh( FILE_MMESH,       // file name (TRIANGLE format)
                  Vertex :: Set_BC, // set markers for vertices
                  Edge   :: Set_BC, // set markers for edges
                  Quad   :: Set_BC, // set markers for quads
                  1) ;              // FORTRAN MESH

  file_out << "** read_mesh(\"" << FILE_MMESH << "\",...) **"
           << endl ;
  mesh.print_marker( file_out ) ;
  cout << "done" << endl ;
}

// test built-in iterators and print markers
void Mesh::print_marker(ostream & outs) const {

  CIterator<Vertex> v( *this ) ;
  CIterator<Edge>   e( *this ) ;
  CIterator<Quad>   q( *this ) ;
  
  // print markers on all (boundary+internal) vertices
  outs << endl << "Print vertex markers" << endl ;
  foreach( v )
    outs << "vertex : "
         << local_number( *v )
    	 << " -- marker : " << v -> marker() << endl ;
  
  // print markers on all (boundary+internal) edge
  outs << endl << "Print edge markers" << endl ;
  foreach( e ) {
    outs << "edge : "
	 << local_number( *e )
	 << " -- marker : " << e -> marker() << endl ;
  }
  
  // print markers on all (boundary+internal) vertices
  outs << endl << "Print quad markers" << endl ;
  foreach ( q ) {
    outs << "quad : "
	 << local_number( *q )
	 << " -- marker : " << q -> marker()
         << endl ;
  }
  
  outs << endl ;
}
 
int
main() {

  ofstream file_out("testall_marker.out") ;
  
  file_out << "****************************************************" << endl ;
  test_build_mesh(file_out) ;
  test_read_mesh (file_out) ;
  file_out << "****************************************************" << endl ;
  file_out.close() ;
  cout << "all tests done" << endl ;
}

// end of file testall_marker.hh
// by Enrico Bertolazzi & Gianmarco Manzini

