/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  P2MESH : MODULE FOR PRINTING INTERNAL MESH STRUCTURE                    |
 |                                                                          |
 |  date         : 1999, 11 October                                         |
 |  version      : 1.0                                                      |
 |  file         : p2print.hh                                               |
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
 |   module for printing the mesh structure (only for testing usage)        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

# ifndef P2PRINT_HH
# define P2PRINT_HH

# ifndef P2MESH_HH
  # include "p2mesh.hh"
# endif

template <typename COMMON>
class p2_print_mesh {
  
  typedef typename COMMON::P2V      P2V ;
  typedef typename COMMON::P2E      P2E ;
  typedef typename COMMON::P2P      P2P ;
  typedef typename COMMON::P2M      P2M ;

  typedef typename COMMON::Unsigned Unsigned ;
  typedef typename COMMON::Integer  Integer ;
  typedef typename COMMON::Real     Real ;

  P2M const & the_mesh(void) const 
    { return static_cast<P2M const &>(*this) ; }
   
public:

  static inline Real round_zero(Real const & a) { 
    if ( a < 1e-10 && a > -1e-10 ) return 0 ;
    else                           return a ;
  }

  void print_all( ostream & outs ) const ;
  void print_all( ostream & outs, P2V const & v ) const ;
  void print_all( ostream & outs, P2E const & e ) const ;
  void print_all( ostream & outs, P2P const & p ) const ;
} ;

// *****************
// * CLASS METHODS *
// *****************

template <typename COMMON> inline
void
p2_print_mesh<COMMON>::print_all(ostream & outs, P2V const & v) const {

  outs << endl
       << "Vertex:\t" << the_mesh() . local_number(v)           // method 74
       << " -- ( " << v . x()                                   // method 10
       << " , "    << v . y()                                   // method 11
       << " )"     << endl ;
  
  if ( COMMON::List ) {
    outs << endl
         << "\t" << v . n_vertex() << " adjacent vertices"      // method  1
	 << endl ;
    for ( Unsigned nv = 0 ; nv < v . n_vertex() ; ++nv ) {      // method  1
      P2V const & adjv = v . vertex(nv) ;                       // method  4 
      outs << "\tvertex(" << nv << ")"
	   << " loc = "  << v . local_number( adjv )            // method  7
           << " glob = " << the_mesh() . local_number( adjv )   // method 74
           << endl ;
    }

    outs << endl 
	 << "\t" << v . n_edge() << " adjacent edge(s)"         // method  2
         << endl ;
    for ( Unsigned ne = 0 ; ne < v . n_edge() ; ++ne ) {        // method  2
      P2E const & adje = v . edge(ne) ;                         // method  5
      outs << "\tedge(" << ne << ")"
	   << " loc = "  << v . local_number( adje )            // method  8
           << " glob = " << the_mesh() . local_number( adje )   // method 75
	   << endl ;
    }

    outs << endl
         << "\t" << v . n_poly() << " adjacent polygon(s)"      // method  3
	 << endl ;
    for ( Unsigned np = 0 ; np < v . n_poly() ; ++np ) {        // method  3
      P2P const & adjp = v . poly(np) ;                         // method  6
      outs << "\tpoly(" << np << ")"
	   << " loc = "  << v . local_number( adjp )            // method  9
           << " glob = " << the_mesh() . local_number( adjp )   // method 76
	   << endl ;
    }
  }
}

template <typename COMMON> inline
void
p2_print_mesh<COMMON>::print_all( ostream & outs, P2E const & e ) const {
  
  Real const t0 = 0.25 ;
  Real const t1 = 0.75 ;
  
  outs << endl
       << "Edge:\t" << the_mesh() . local_number(e)             // method 75
       << endl << endl ;
  
  for ( Unsigned nv = 0 ; nv < e . n_vertex() ; ++nv ) {        // method 12
    P2V const & adjv = e . vertex( nv ) ;                       // method 15
    outs << "\tVertex(" << nv << ")" 
	 << " loc = "   << e . local_number( adjv )             // method 18
	 << " glob = "  << the_mesh() . local_number( adjv )    // method 74
	 << " ( " << round_zero( e . x(nv) )                    // method 22
	 << " , " << round_zero( e . y(nv) )                    // method 23
	 << " )" << endl ;
  }

  if ( ! e . ok_poly(1) )                                       // method 21
    outs << "\tedge is on the boundary" << endl ; 

  for ( Unsigned ip = 0 ; ip < e . n_poly() ; ++ip) {           // method 14
    P2P const & adjp = e . poly(ip) ;                           // method 17
    switch (ip) {
    case 0 : outs << endl << "\tLeft polygon " ; break ;
    case 1 : outs << endl << "\tRigth polygon" ; break ; }
    outs << " loc = "  << e . local_number( adjp )              // method 20
         << " glob = " << the_mesh() . local_number( adjp )     // method 76
	 << endl ;
  }

  for ( Unsigned ne = 0 ; ne < e . n_edge() ; ++ne) {           // method 13
    P2E const & adje = e . edge( ne ) ;                         // method 16 
    outs << "\tadj edge(" << ne << ")  " 
         << " loc = "  << e . local_number( adje )              // method 19
         << " glob = " 
         << the_mesh() . local_number( adje )                   // method 72
         << endl ;
  }

  outs << endl
       << "\tMidpoint:            "
       << "( "  << round_zero( e . xm() )                       // method 24
       << " , " << round_zero( e . ym() )                       // method 25  
       << " )"
       << endl
       << "\t(1/4 and 3/4) point: " 
       << "( "  << round_zero( e . xt( t0 ) )                   // method 26
       << " , " << round_zero( e . yt( t0 ) )                   // method 27
       << " )"  << " , "
       << "( "  << round_zero( e . xt( t1 ) )                   // method 26
       << " , " << round_zero( e . yt( t1 ) )                   // method 27
       << " )"
       << endl
       << "\tOrthogonal vector:   " 
       << "[ "  << round_zero( e . nx() )                       // method 28
       << " , " << round_zero( e . ny() )                       // method 29
       << " ]"
       << endl
       << "\tTangent vector:      " 
       << "[ "  << round_zero( e . tx() )                       // method 30
       << " , " << round_zero( e . ty() )                       // method 31
       << " ]"
       << endl
       << "\tLength:              "
       << e . length()                                          // method 32
       << endl ;
}

template <typename COMMON> inline
void
p2_print_mesh<COMMON>::print_all( ostream & outs, P2P const & p ) const {

  outs << "Poly: " << the_mesh() . local_number( p )            // method 76
       << endl
       << "\tCentroid: ( " << round_zero( p . xc() )            // method 55
       << " , "            << round_zero( p . yc() )            // method 56
       << " )"
       << endl 
       << "\tArea:     "   << round_zero( p . area() )          // method 57
       << endl << endl ;
  
  for ( Unsigned nv = 0 ; nv < p . n_vertex() ; ++nv ) {        // method 33
    P2V const & adjv = p . vertex( nv ) ;                       // method 36
    outs << "\tvertex(" << nv << ")"
	 << " loc = "   << p . local_number( adjv )             // method 39
	 << " glob = "  << the_mesh() . local_number( adjv )    // method 74
	 << " coor ( "  << p . x( nv )                          // method 44
	 << " , "       << p . y( nv )                          // method 45
	 << " )" << endl ;
  }

  Real const t0 = 0.25 ;
  Real const t1 = 0.75 ;
  for ( Unsigned ne = 0 ; ne < p . n_edge() ; ++ne ) {          // method 34
    P2E const & adje = p . edge(ne) ;                           // method 37
    outs << endl
         << "\tedge("  << ne << ")"
         << " loc = "  << p . local_number( adje )              // method 40
	 << " glob = " << the_mesh() . local_number( adje )     // method 75
	 << " orientation " << p . ok_oriented(ne)              // method 43
	 << endl 
	 << "\tMidpoint:          ( "
	 << round_zero( p . xm( ne ) ) << " , "                 // method 46
	 << round_zero( p . ym( ne ) ) << " )"                  // method 47
	 << endl
         << "\t1/4 & 3/4 nodes:   ( "
	 << round_zero( p . xt( ne, t0 ) ) << " , "             // method 48
	 << round_zero( p . yt( ne, t0 ) ) << " ),  ( "         // method 49
	 << round_zero( p . xt( ne, t1 ) ) << " , "             // method 48
	 << round_zero( p . yt( ne, t1 ) ) << " )"              // method 49
	 << endl
         << "\tOrthogonal vector: [ "
	 << round_zero( p . nx( ne ) ) << " , "                 // method 50
	 << round_zero( p . ny( ne ) ) << " ]"                  // method 51
	 << endl
         << "\tTangent vector:    [ "
	 << round_zero( p . tx( ne ) ) << " , "                 // method 52
	 << round_zero( p . ty( ne ) ) << " ]"                  // method 53
	 << endl
         << "\tLength:              "
	 << round_zero( p . length( ne ) )                      // method 54
	 << endl ;
  }

  outs << endl ;

  for ( Unsigned np = 0 ; np < p . n_poly() ; ++np ) {          // method 35
    outs << "\tpoly(" << np << ")" ;
    if ( p . ok_poly( np ) ) {                                  // method 42
      P2P const & adjp = p . poly(np) ;                         // method 38
      outs << " loc = "  << p . local_number( adjp )            // method 41
	   << " glob = " << the_mesh() . local_number( adjp )   // method 76
	   << endl ;
    } else {
      outs << " not existing" << endl ;
    }
  }
 
  outs << endl
       << "\tMapping polygon: actual ==> reference ==> actual "
       << endl ;
 
  Unsigned const _size = COMMON::Size ;
  for ( Unsigned i = 0 ; i < _size ; ++i ) {
    Real _s, _t, _x1, _y1 ;
    Real _x = p . xm( i ) ;                                     // method 46
    Real _y = p . ym( i ) ;                                     // method 47

    p . xy_to_st( _x, _y, _s,  _t) ;                            // method 59
    p . st_to_xy( _s, _t, _x1, _y1) ;                           // method 58

    outs << "\t(x,y) = (" << round_zero( _x )
         << "," << round_zero( _y ) << ") ==> " 
	 <<   "(s,t) = (" << round_zero( _s )
	 << "," << round_zero( _t ) << ") ==> " 
         <<   "(x,y) = (" << round_zero( _x1 )
         << "," << round_zero( _y1 ) << ")"
	 << endl ;
  }
  
  Real J[2][2], iJ[2][2] ;
  p . jacobian(0,0,J) ;                                         // method 60
  p . inverse_jacobian(0,0,iJ) ;                                // method 61

  outs << endl
       << "\tJacobian:\t" << setw(8) << round_zero(J[0][0]) << "\t"
                          << setw(8) << round_zero(J[0][1]) << endl
       << "\t         \t" << setw(8) << round_zero(J[1][0]) << "\t"
                          << setw(8) << round_zero(J[1][1]) << endl
       << endl
       << "\tInverse: \t" << setw(8) << round_zero(iJ[0][0]) << "\t"
                          << setw(8) << round_zero(iJ[0][1]) << endl
       << "\t         \t" << setw(8) << round_zero(iJ[1][0]) << "\t"
                          << setw(8) << round_zero(iJ[1][1]) << endl
       << endl ;
}

template <typename COMMON> inline
void
p2_print_mesh<COMMON>::print_all( ostream & outs ) const {
  Unsigned i ;
  
  the_mesh() . report( outs ) ;                                 // method 84

  outs << endl << "Vertices: " 
       << the_mesh() . n_vertex()  << " (total) "               // method 62
       << the_mesh() . n_bvertex() << " (boundary) "            // method 63
       << the_mesh() . n_ivertex() << " (internal) "            // method 64
       << endl << "Edges:    " 
       << the_mesh() . n_edge()    << " (total) "               // method 65
       << the_mesh() . n_bedge()   << " (boundary) "            // method 66
       << the_mesh() . n_iedge()   << " (internal) "            // method 67
       << endl << "Polygons: " 
       << the_mesh() . n_poly()    << " (total) "               // method 68
       << the_mesh() . n_bpoly()   << " (boundary) "            // method 69
       << the_mesh() . n_ipoly()   << " (internal) "            // method 70
       << endl << endl ;

  for ( i = 0 ; i < the_mesh() . n_vertex() ; ++i ) {           // method 62
                                                                // method 74 
                                                                // method 71
    if ( i != the_mesh() . local_number( the_mesh() . vertex(i) ) ) {
      cerr << "Vertex numbering inconsistency detected !\n" ;
    }
  }
  outs << "Vertex  numbering is consistent" << endl ;

  for ( i = 0 ; i < the_mesh() . n_edge() ; ++i )  {            // method 65
                                                                // method 75 
                                                                // method 72
    if ( i != the_mesh() . local_number( the_mesh() . edge(i) ) )
      cerr << "Edge numbering inconsistency detected !\n" ;
  }
  outs << "Edge    numbering is consistent" << endl ;
    
  for ( i = 0 ; i < the_mesh() . n_poly() ; ++i )  {            // method 68
                                                                // method 76
                                                                // method 73
    if ( i != the_mesh() . local_number( the_mesh() . poly(i) ) )
      cerr << "Polygon numbering inconsistency detected !\n" ;
  }
  outs << "Polygon numbering is consistent" << endl ;

  bool ok_mesh = the_mesh() . test_mesh() ;                     // method 85
  if ( ok_mesh ) 
    outs << "BUILT-IN CONSISTENCY CHECK : MESH IS OK" << endl ;
  else
    outs << "BUILT-IN CONSISTENCY CHECK : INCONSISTENCY DETECTED" << endl ;

  Real xmin,ymin, xmax,ymax ;
  the_mesh() . bbox( xmin,ymin, xmax,ymax ) ;                   // method 77
  outs << endl
       << "Bounding box: " 
       << setw(8) << round_zero(xmin) << "  "
       << setw(8) << round_zero(ymin) << "  " 
       << setw(8) << round_zero(xmax) << "  "
       << setw(8) << round_zero(ymax) << endl ;

  // print info about vertices
  CIterator<P2V> iv( the_mesh() ) ;                             // method 108
  foreach( iv ) print_all( outs, *iv ) ;

  // print info about edges
  CIterator<P2E> ie( the_mesh() ) ;                             // method 108
  foreach( ie ) print_all( outs, *ie ) ;

  // print info about triangle
  CIterator<P2P> ip( the_mesh() ) ;                             // method 108
  foreach( ip ) print_all( outs, *ip ) ;

  // built-in info proc
  outs << endl << "Built-in info for meshes: " << endl ;
  the_mesh() . print( outs, 0 ) ;                               // method 86
  outs << endl << "Built-in info method for vertices: " << endl ;
  foreach( iv ) the_mesh() . print( outs, *iv, 0 ) ;            // method 87
  outs << endl << endl << "Built-in info method for edges: " << endl ;
  foreach( ie ) the_mesh() . print( outs, *ie, 0 ) ;            // method 88
  outs << endl << endl << "Built-in info method for polygons: " << endl ;
  foreach( ip ) the_mesh() . print( outs, *ip, 0 ) ;            // method 89

}

# endif

// end of file p2print.hh
// by Enrico Bertolazzi & Gianmarco Manzini
