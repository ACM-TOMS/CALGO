/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  P2MESH : DRIVER FOR TESTING THE PACKAGE ON QUADRILATERAL-BASED MESHES   |
 |                                                                          |
 |  date         : 1999, 11 October                                         |
 |  version      : 1.0                                                      |
 |  file         : testall_quad.cc                                          |
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
 |    test all the methods of P2MESH package on quadrilateral-based meshes. |
 |                                                                          |
 |  output file:                                                            |
 |    testall_quad.out   : this ASCII file is generated if the code         |
 |                         is compiled with -DLIST=false                    | 
 |    testall_quad_l.out : this ASCII file is generated if the code         |
 |                        is compiled with -DLIST=true                      |
 |                                                                          |
 | ``LIST = true'' forces the construction of vertex topology lists         |
 |                                                                          |
 |  documentation:                                                          |
 |    userman.ps -- user manual for P2MESH programmers                      |
 |                                                                          |
\*--------------------------------------------------------------------------*/

# ifndef FILE_QMESH
# define FILE_QMESH "meshes/qmesh"
# endif

# ifndef FILE_SMESH
# define FILE_SMESH "meshes/mesh.grd"
# endif

# ifndef LIST
  # define LIST true
# endif

# include "p2mesh.hh"
# include "p2print.hh"

# include <stdio.h>
# include <math.h>

// declare the name of the user defined classes
class Vertex ;
class Edge   ;
class Quad   ;
class Mesh   ;

//setup common code & data

class Common : public p2_common<Vertex,Edge,Quad,Mesh,4,LIST>
{ } ; 

// define class Vertex
class Vertex : public p2_vertex<Common> {
  unsigned v_marker ;
public:
  static void Set_Marker(Vertex & V, unsigned const & marker ) ; 
  unsigned marker(void) const { return v_marker ; }
} ;

// define class Edge
class Edge : public p2_edge<Common> {
  unsigned e_marker ;
public:
  static void Set_Marker(Edge & E, unsigned const & marker ) ; 
  unsigned marker(void) const { return e_marker ; }
} ;

// define class Quad
class Quad : public p2_poly<Common> {
  unsigned q_marker ;
public:
  static void Set_Marker(Quad & Q, unsigned const & marker ) ; 
  unsigned marker(void) const { return q_marker ; }
} ;

// define class Mesh
class Mesh : public p2_mesh<Common>,
	     public p2_print_mesh<Common>
{ public: void print_marker(ostream & outs) const ; } ;

// ************************
// * Vertex CLASS METHODS *
// ************************
void
Vertex::Set_Marker(Vertex & V, unsigned const & m ) 
{ V . v_marker = m ; }

// **********************
// * Edge CLASS METHODS *
// **********************

void
Edge::Set_Marker(Edge & E, unsigned const & m ) 
{ E . e_marker = m ; }

// **************************
// * Quad CLASS METHODS *
// **************************

void
Quad::Set_Marker(Quad & Q, unsigned const & m ) 
{ Q . q_marker = m ; }

// **********************
// * Mesh CLASS METHODS *
// **********************

static
void
test_tensor_mesh(ostream & file_out) {

  cout << "testing tensor_mesh..." ; cout . flush() ;
  
  Mesh mesh ; // instantiate an "empty" mesh object

  unsigned const nx = 2 ; 
  unsigned const ny = 2 ;

  double xmin = 0.5 ;
  double xmax = 1.5 ;

  double ymin = 0.5 ;
  double ymax = 1.5 ;
  
  // method 78
  mesh . tensor_mesh( xmin, xmax, ymin, ymax, // bounding box 
                      nx, ny,                 // box partitioning
                      Vertex :: Set_Marker,   // set markers for vertices
                      Edge   :: Set_Marker,   // set markers for edges
                      Quad   :: Set_Marker) ; // set markers for quads

  file_out << endl
           << "*****************" << endl
           << "** tensor_mesh **" << endl
           << "*****************" << endl
           << "NX = " << nx << " NY = " << ny << endl
           << endl ;

  mesh . print_all( file_out ) ;
  mesh . print_marker( file_out ) ;
  cout << "done" << endl ;
}

static
void
test_std_tensor_mesh(ostream & file_out) {

  cout << "testing std_tensor_mesh..." ; cout . flush() ;
  
  Mesh mesh ; // instantiate an "empty" mesh object

  unsigned const nx = 2 ; 
  unsigned const ny = 2 ;
  
  // method 79
  mesh . std_tensor_mesh( nx, ny,                 // box partitioning
                          Vertex :: Set_Marker,   // set markers for vertices
                          Edge   :: Set_Marker,   // set markers for edges
                          Quad   :: Set_Marker) ; // set markers for quads

  
  file_out << endl
           << "*********************" << endl
           << "** std_tensor_mesh **" << endl
           << "*********************" << endl
           << "NX = " << nx << " NY = " << ny << endl
           << endl ;

  mesh . print_all( file_out ) ;
  mesh . print_marker( file_out ) ;
  cout << "done" << endl ;
}

static
void
shape(double const & s, double const & t, 
      double       & x, double       & y ) {
  static double const PI2 = 2 * atan( 1.0 ) ; // pi/2
  x = ( 1 + t ) * sin( s * PI2 ) ;
  y = ( 1 + t ) * cos( s * PI2 ) ;
}

static
void
test_map_mesh(ostream & file_out) {

  cout << "testing map_mesh..." ; cout . flush() ;
  
  Mesh mesh ; // instantiate an "empty" mesh object
  unsigned const nx = 2 ; 
  unsigned const ny = 2 ;

  // method 80
  mesh . map_mesh( shape,                  // bounding box 
                   nx, ny,                 // box partitioning
                   Vertex :: Set_Marker,   // set markers for vertices
                   Edge   :: Set_Marker,   // set markers for edges
                   Quad   :: Set_Marker) ; // set markers for quads

  file_out << endl
           << "**************" << endl
           << "** map_mesh **" << endl
           << "**************" << endl
           << "NX = " << nx << " NY = " << ny << endl
           << endl ;

  mesh . print_all( file_out ) ;
  mesh . print_marker( file_out ) ;
  cout << "done" << endl ;
}

static
void
test_read_map_mesh(ostream & file_out) {

  cout << "testing read_map_mesh..." ; cout . flush() ;
  
  Mesh mesh ; // instantiate an "empty" mesh object

  // method 81
  mesh . read_map_mesh( FILE_SMESH,             // file name (TRIANGLE format)
                        Vertex :: Set_Marker,   // set markers for vertices
                        Edge   :: Set_Marker,   // set markers for edges
                        Quad   :: Set_Marker) ; // set markers for quads
      
  file_out << endl
           << "*******************" << endl
           << "** read_map_mesh **" << endl
           << "*******************" << endl
           << " FILE = " << FILE_SMESH << endl
           << endl ;
  mesh . print_all( file_out ) ;
  mesh . print_marker( file_out ); 
  cout << "done" << endl ;
}

static
void
test_build_mesh(ostream & file_out) {

  cout << "testing build_mesh..." ; cout . flush() ;
  
  unsigned const nv = 7 ;
  double  XY[2*nv]  = { 0.0, 0.0,    1.0, 0.0,    2.0, 0.0,    2.0, 2.0,    
			0.0, 2.0,    0.0, 1.0,    1.0, 1.0 } ;
  unsigned mv[nv]   = { 1, 2, 3, 4, 5, 6, 7} ;


  unsigned const ne = 9 ;  
  unsigned E[2*ne]  = { 0, 1,    1, 2,    2, 3,    3, 4,    4, 5,    
			5, 0,    1, 6,    6, 5,    6, 3 } ;    
  unsigned me[ne]   = { 1, 2, 3, 4, 5, 6, 7, 8, 9} ;

  unsigned const np = 3 ;  
  unsigned P[4*np]  = { 0, 1, 6, 5,    1, 2, 3, 6,    5, 6, 3, 4 } ;
  unsigned mp[np]   = { 1, 2, 3 } ; 
  
  Mesh mesh ; // instantiate an "empty" mesh object
  
  // method 82
  mesh . build_mesh( nv, XY, mv, Vertex :: Set_Marker, // vertex arrays
		     ne,  E, me, Edge   :: Set_Marker, // edge   arrays
		     np,  P, mp, Quad   :: Set_Marker, // quad   arrays
		     0 ) ;                             // C MESH offset 
    
  file_out << endl
           << "****************" << endl
           << "** build_mesh **" << endl
           << "****************" << endl
           << endl ;
  mesh . print_all( file_out ) ;
  mesh . print_marker( file_out );
  cout << "done" << endl ;
}

static
void
test_read_mesh(ostream & file_out) {

  cout << "testing read_mesh..." ; cout . flush() ;

  Mesh mesh ; // instantiate an "empty" mesh object

  // method 83
  mesh . read_mesh( FILE_QMESH,             // file name (TRIANGLE format)
                    Vertex :: Set_Marker,   // set markers for vertices
                    Edge   :: Set_Marker,   // set markers for edges
                    Quad   :: Set_Marker,   // set markers for quads
                    1) ;                    // FORTRAN MESH
  file_out << endl
           << "***************" << endl
           << "** read_mesh **" << endl
           << "***************" << endl
           << " FILE = " << FILE_QMESH << endl
           << endl ;
  mesh . print_all( file_out ) ;
  mesh . print_marker( file_out ) ;
  cout << "done" << endl ;
}

// test built-in iterators and print markers
void Mesh::print_marker(ostream & outs) const {
 
  unsigned const All      = 0 ;
  unsigned const Boundary = 1 ;
  unsigned const Internal = 2 ;

  // print markers on all (boundary+internal) vertices
  outs << endl << "Print vertex markers" << endl ;
  CIterator<Vertex> v( *this, All ) ;                         // method 109
  foreach( v ) {
    Vertex const & current_vertex = *v ;                      // method 116
    outs << "vertex : "                                      
	 << this -> local_number( current_vertex )            // method 74
	 << " -- marker : " << current_vertex . marker() << endl ;
  }
  
  // print markers on all (boundary+internal) edge
  outs << endl << "Print edge markers" << endl ;
  CIterator<Edge> e( *this, Boundary ) ;                       // method 109
  foreach( e ) {
    Edge const & current_edge = *e ;                           // method 116
    outs << "edge (bnd) : "                                      
	 << this -> local_number( current_edge )               // method 74
	 << " -- marker : " << current_edge . marker() << endl ;
  }
  e . set_loop( *this, Internal ) ;                            // method 112
  foreach( e ) {
    Edge const & current_edge = *e ;                           // method 116
    outs << "edge (int) : "                                      
	 << this -> local_number( current_edge )               // method 74
	 << " -- marker : " << current_edge . marker() << endl ;
  }
  
  // print markers on all (boundary+internal) vertices
  outs << endl << "Print quad markers" << endl ;
  CIterator<Quad> q ;                                         // method 110
  q . set_loop( *this ) ;                                     // method 111
  for ( q . begin() ;                                         // method 113
	! q . end_of_loop() ;                                 // method 114
	++q ) {                                               // method 115
    Quad const & current_quad = *q ;                          // method 116
    outs << "quad : "                                      
	 << this -> local_number( current_quad )              // method 74
	 << " -- marker : " << current_quad . marker() << endl ;
  }
  outs << endl ;
}
 
int
main() {
  ofstream file_out( LIST ? "testall_quad_l.out" : "testall_quad.out") ;
  
  test_std_tensor_mesh(file_out) ;
  test_tensor_mesh    (file_out) ;
  test_map_mesh       (file_out) ;
  test_build_mesh     (file_out) ;
  test_read_mesh      (file_out) ;
  test_read_map_mesh  (file_out) ;

  file_out . close() ;
  cout << "all tests done" << endl ;
}

// end of file testall_quad.cc
// by Enrico Bertolazzi & Gianmarco Manzini
