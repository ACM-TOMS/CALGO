/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  P2MESH : DRIVER FOR TESTING THE PACKAGE ON TRIANGLE-BASED MESHES        |
 |                                                                          |
 |  date         : 1999, 11 October                                         |
 |  version      : 1.0                                                      |
 |  file         : testall_triangle.cc                                      |
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
 |    test all the methods of P2MESH package on triangle-based meshes.      |
 |                                                                          |
 |  output file:                                                            |
 |    testall_triangle.out   : this ASCII file is generated if the code     |
 |                             is compiled with -DLIST=false                | 
 |    testall_triangle_l.out : this ASCII file is generated if the code     |
 |                             is compiled with -DLIST=true                 |
 |                                                                          |
 | ``LIST = true'' forces the construction of vertex topology lists         |
 |                                                                          |
 |  documentation:                                                          |
 |    userman.ps -- user manual for P2MESH programmers                      |
 |                                                                          |
\*--------------------------------------------------------------------------*/

# ifndef FILE_TMESH
# define FILE_TMESH "meshes/mesh"
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
class Triangle ;
class Mesh   ;

//setup common code & data
class Common : public p2_common<Vertex,Edge,Triangle,Mesh,3,LIST>
{ } ; 

// define class Vertex
class Vertex : public p2_vertex<Common> {
  unsigned v_marker ;
public:
  static void Set_BC(Vertex & V, unsigned const & marker ) ; 
  unsigned marker(void) const { return v_marker ; }
} ;

// define class Edge
class Edge : public p2_edge<Common> {
  unsigned e_marker ;
public:
  static void Set_BC(Edge & E, unsigned const & marker ) ; 
  unsigned marker(void) const { return e_marker ; }
} ;

// define class Triangle
class Triangle : public p2_poly<Common> {
  unsigned t_marker ;
public:
  static void Set_BC(Triangle & T, unsigned const & marker ) ; 
  unsigned marker(void) const { return t_marker ; }
} ;

// define class Mesh
class Mesh : public p2_mesh<Common>,
	     public p2_print_mesh<Common>
{ public: void print_marker(ostream & file_out) const ; } ;

// ************************
// * Vertex CLASS METHODS *
// ************************
void
Vertex::Set_BC(Vertex & V, unsigned const & m ) 
{ V . v_marker = m ; }

// **********************
// * Edge CLASS METHODS *
// **********************

void
Edge::Set_BC(Edge & E, unsigned const & m ) 
{ E . e_marker = m ; }

// **************************
// * Triangle CLASS METHODS *
// **************************

void
Triangle::Set_BC(Triangle & T, unsigned const & m ) 
{ T . t_marker = m ; }

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

  double xmin = 1.0 ;
  double xmax = 2.0 ;

  double ymin = 1.0 ;
  double ymax = 2.0 ;

  for ( unsigned kind = 0 ; kind < 3 ; ++kind ) {

    // method 78
    mesh . tensor_mesh( xmin, xmax, ymin, ymax, // bounding box 
                        nx, ny,                 // box partitioning
                        Vertex   :: Set_BC,     // set markers for vertices
                        Edge     :: Set_BC,     // set markers for edges
                        Triangle :: Set_BC,     // set markers for triangles
                        kind ) ;                // mesh kind

    file_out << endl
             << "*****************" << endl
             << "** tensor_mesh **" << endl
             << "*****************" << endl
             << "KIND = " << kind << endl
             << "NX = " << nx << " NY = " << ny << endl
             << endl ;

    mesh . print_all( file_out ) ;
    mesh . print_marker( file_out ) ;

  }
  cout << "done" << endl ;
}

static
void
test_std_tensor_mesh(ostream & file_out) {

  cout << "testing std_tensor_mesh..." ; cout . flush() ;
  
  Mesh mesh ; // instantiate an "empty" mesh object

  unsigned const nx = 2 ; 
  unsigned const ny = 2 ;

  for ( unsigned kind = 0 ; kind < 3 ; ++kind ) {
  
    // method 79
    mesh . std_tensor_mesh( nx, ny,             // box partitioning
                            Vertex   :: Set_BC, // set markers for vertices
                            Edge     :: Set_BC, // set markers for edges
                            Triangle :: Set_BC, // set markers for triangles
                            kind ) ;            // mesh kind

    file_out << endl
             << "*********************" << endl
             << "** std_tensor_mesh **" << endl
             << "*********************" << endl
             << "KIND = " << kind << endl
             << "NX = " << nx << " NY = " << ny << endl
             << endl ;

    mesh . print_all( file_out ) ;
    mesh . print_marker( file_out ) ;
    
  }
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

  for ( unsigned kind = 0 ; kind < 3 ; ++kind ) {
  
    // method 80
    mesh . map_mesh( shape,              // bounding box 
                     nx, ny,             // box partitioning
                     Vertex   :: Set_BC, // set markers for vertices
                     Edge     :: Set_BC, // set markers for edges
                     Triangle :: Set_BC, // set markers for triangles
                     kind ) ;            // mesh kind

    file_out << endl
             << "**************" << endl
             << "** map_mesh **" << endl
             << "**************" << endl
             << "KIND = " << kind << endl
             << "NX = " << nx << " NY = " << ny << endl
             << endl ;

    mesh . print_all( file_out ) ;
    mesh . print_marker( file_out ) ;
  }
  cout << "done" << endl ;
}

static
void
test_read_map_mesh(ostream & file_out) {

  cout << "testing read_map_mesh..." ; cout . flush() ;
  
  Mesh mesh ; // instantiate an "empty" mesh object

  for ( unsigned kind = 0 ; kind < 3 ; ++kind ) {

    // method 81
    mesh . read_map_mesh( FILE_SMESH,          // file name
                          Vertex   :: Set_BC,  // set markers for vertices
                          Edge     :: Set_BC,  // set markers for edges
                          Triangle :: Set_BC,  // set markers for triangles
                          kind ) ;             // mesh kind
      
    file_out << endl
             << "**************" << endl
             << "** map_mesh **" << endl
             << "**************" << endl
             << "KIND = " << kind << endl
             << "FILE = " << FILE_SMESH << endl
             << endl ;

    mesh . print_all( file_out ) ;
    mesh . print_marker( file_out ) ;
  }
  cout << "done" << endl ;
}

static
void
test_build_mesh(ostream & file_out) {

  cout << "testing build_mesh..." ; cout . flush() ;

  unsigned const nv = 9 ;
  double XY[2*nv]   = {  0.0, 0.0,   0.5,  0.0,   0.5, 0.5,   
			 1.0, 0.5,   1.0,  1.0,   0.5, 1.0,   
			 0.0, 1.0,   0.25, 0.5,   0.5, 0.75 } ;
  unsigned mv[nv]   = { 1, 2, 3, 4, 5, 9, 6, 7, 8 } ;
 
  unsigned const ne = 16 ;
  unsigned E[2*ne]  = {  0, 1,   1, 2,   2, 3,   3, 4,    
                         4, 5,   5, 6,   6, 7,   7, 0,    
			 1, 7,   7, 8,   8, 6,   8, 4,    
			 5, 8,   2, 7,   3, 8,   8, 2 } ; 
  unsigned me[ne]   = { 1, 10, 12, 16,  8, 15,  4,  3, 
			2,  5,  6,  7,  9, 11, 13, 14 } ;
  

  unsigned const np = 8 ;
  unsigned P[3*np]  = {  0, 1, 7,   1, 2, 7,   2, 3, 8,   4, 8, 3,
			 8, 4, 5,   8, 5, 6,   6, 7, 8,   2, 8, 7 } ; 
  unsigned mp[np]   = { 1, 4, 5, 7, 3, 6, 1, 8 } ;

  Mesh mesh ; // instantiate an "empty" mesh object
  
  // method 82
  mesh . build_mesh( nv, XY, mv, Vertex   :: Set_BC, // vertex   arrays
		     ne,  E, me, Edge     :: Set_BC, // edge     arrays
		     np,  P, mp, Triangle :: Set_BC, // triangle arrays
		     0 ) ;                           // C MESH offset 
    
  file_out << endl
           << "**************" << endl
           << "** build_mesh **" << endl
           << "**************" << endl
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
  mesh . read_mesh( FILE_TMESH,         // file name (TRIANGLE format)
                    Vertex   :: Set_BC, // set markers for vertices
                    Edge     :: Set_BC, // set markers for edges
                    Triangle :: Set_BC, // set markers for triangles
                    1) ;                // FORTRAN MESH offset
  file_out << endl
           << "***************" << endl
           << "** read_mesh **" << endl
           << "***************" << endl
           << "FILE = " << FILE_TMESH << endl
           << endl ;

  mesh . print_all( file_out ) ;
  mesh . print_marker( file_out );
 
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
  CIterator<Triangle> t ;                                     // method 110
  t . set_loop( *this ) ;                                     // method 111
  for ( t . begin() ;                                         // method 113
	! t . end_of_loop() ;                                 // method 114
	++t ) {                                               // method 115
    Triangle const & current_triangle = *t ;                  // method 116
    outs << "quad : "                                      
	 << this -> local_number( current_triangle )          // method 74
	 << " -- marker : " << current_triangle . marker() << endl ;
  }
  outs << endl ;
}

int
main() {
  ofstream file_out( LIST ? "testall_triangle_l.out" : "testall_triangle.out" ) ;

  test_std_tensor_mesh(file_out) ;
  test_tensor_mesh    (file_out) ;
  test_map_mesh       (file_out) ;
  test_build_mesh     (file_out) ;
  test_read_mesh      (file_out) ;
  test_read_map_mesh  (file_out) ;

  file_out . close() ;

  cout << "all tests done" << endl ;
}

// end of file testall_triangle.cc
// by Enrico Bertolazzi & Gianmarco Manzini

