/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  P2MESH : a polygon-based mesh manager                                   |
 |                                                                          |
 |  date         : 2001, November 28                                        |
 |  version      : 1.2                                                      |
 |               : 1.1 (1999/12/06)                                         |
 |  file         : p2mesh.hh                                                |
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
 |  documentation: (in the doc directory of the package)                    |
 |                                                                          |
 |    userman.ps -- user manual for P2MESH programmers                      |
 |    kernel.ps  -- internal implementation description                     |
 |    primer.ps  -- beginner introduction with commented examples           |
 |                                                                          |
\*--------------------------------------------------------------------------*/

// DEFINITION PART

# ifndef  __cplusplus
# error You must use C++ for p2mesh.hh
# endif

# ifndef P2MESH_HH
# define P2MESH_HH

// standard C lib (for sqrt)
# include <math.h>

// standard include IO
# include <memory>
# include <iostream>
# include <fstream>
# include <iomanip>

// standard string class
# include <string>

// STL lib
# include <algorithm>
# include <vector>

# ifdef P2MESH_VERBOSE
  # define P2MESH_MSG(A) do { cout << A << endl ; } while(0)
# else
  # define P2MESH_MSG(A)
# endif

# ifndef P2MESH_NOT_USE_STD
  using namespace std ;
# endif

namespace p2_mesh_namespace {
  template <typename ForwardIterator, typename T, typename Compare>
  inline
  bool
  binary_search(ForwardIterator first,
                ForwardIterator last,
                const T& value,
                Compare comp,
                ForwardIterator & i)
  {
    i = lower_bound(first, last, value, comp);
    return i != last && !comp(value, *i);
  }

  template <typename T>
  class cmp_pedge {
  public:
    bool operator() (const T & a, const T & b) {
      if ( &a -> vertex(0) == &b -> vertex(0)) {
        if ( &a -> vertex(1) > &b -> vertex(1) ) return true ;
      } else {
        if ( &a -> vertex(0) > &b -> vertex(0) ) return true ;
      }
      return false ;
    }
  } ;
  
  double const epsi     = 1e-10 ;
  double const big_epsi = 1e-6 ;
}

/*
     #####
    #     #   ####   #    #  #    #   ####   #    #
    #        #    #  ##  ##  ##  ##  #    #  ##   #
    #        #    #  # ## #  # ## #  #    #  # #  #
    #        #    #  #    #  #    #  #    #  #  # #
    #     #  #    #  #    #  #    #  #    #  #   ##
     #####    ####   #    #  #    #   ####   #    #
*/

template <
  typename P2V_type,
  typename P2E_type,
  typename P2P_type,
  typename P2M_type,
  unsigned SIZE_value    = 3,
  bool     LIST_value    = false,
  typename REAL_type     = double,
  typename INTEGER_type  = int,
  typename UNSIGNED_type = unsigned,
  typename VMARK_type    = unsigned,
  typename EMARK_type    = unsigned,
  typename PMARK_type    = unsigned>
class p2_common {
public:

  typedef P2V_type      P2V ;
  typedef P2E_type      P2E ;
  typedef P2P_type      P2P ;
  typedef P2M_type      P2M ;

  typedef VMARK_type    Vmark ;
  typedef EMARK_type    Emark ;
  typedef PMARK_type    Pmark ;

  typedef REAL_type     Real ;
  typedef INTEGER_type  Integer ;
  typedef UNSIGNED_type Unsigned ;

  static unsigned const Size = SIZE_value ;
  static bool     const List = LIST_value ;

  static istream & eatline(istream & s) {
    while ( s.get() != '\n' && s.good() ) {}
    return s ;
  }

  static istream & eatchar(istream & s) { s.get() ; return s ; }
  static istream & eatcomments(istream & s) {
    char c = s.peek() ;
    while ( ( c == '!' || c == '%' || c == '#' || c == ';' || c == '$')
            && s.good() ) { s >> eatline ; c = s.peek() ; }
    return s ;
  }
  
  static inline Real abs(Real const & a) { return a > 0 ? a : -a ; }
  
  static inline void msg_error(char const method[], char const msg[]) { 
    cerr << endl << endl
         << "P2MESH in method: ``" << method << "''" << endl
         << "Fatal error: " << msg << endl << endl ;
    exit(0) ;
  }

# ifdef P2MESH_DEBUG

  static void test_ok(bool ok, char const method[],
                      Unsigned const i, char const msg[]) {
    if ( !ok ) {
      cerr << endl << endl
           << "P2MESH in method: ``" << method << "(" << i << ")''" << endl
           << "Fatal error: " << msg << endl << endl ;
      exit(0) ;
    }
  }

  static void check_range(Unsigned const i,
                          Unsigned const imax,
                          char const method[]) {
    if ( i >= imax ) {
      cerr << endl << endl
           << "P2MESH in method: ``" << method << "(" << i << ")''" << endl 
           << "Fatal error: local index out of range." << endl << endl ;
      exit(0) ;
    }
  }

  static void check_grange(Unsigned const i,
                           Unsigned const imax,
                           char const method[]) {
    if ( i >= imax ) {
      cerr << endl << endl
           << "P2MESH in method: ``" << method << "(" << i << ")''" << endl 
           << "Fatal error: global index out of range." << endl << endl ;
      exit(0) ;
    }
  }
# else
  static inline void
  check_range(Unsigned const, Unsigned const, char const []) {} ;
  static inline void
  check_grange(Unsigned const, Unsigned const, char const []) {} ;
  static inline void
  test_ok(bool, char const [], Unsigned const, char const []) {} ;
# endif
} ;

template <unsigned P2_LIST>   class p2_vertex_variant ;
template <typename P2_COMMON> class p2_vertex ;
template <typename P2_COMMON> class p2_edge ;
template <typename P2_COMMON> class p2_poly ;
template <typename P2_COMMON> class p2_mesh ;

template <>
class p2_vertex_variant<1> {
protected:
  void     **psV_, **psE_, **psP_ ;
  unsigned p2_nv, p2_ne, p2_np ;

  void Reset(void) { p2_nv = p2_ne = p2_np = 0 ; }

  unsigned Assign(void ** ptr) {
    unsigned ntot = p2_nv + p2_ne + p2_np ;
    psV_ = ptr ;
    psE_ = psV_ + p2_nv ;
    psP_ = psE_ + p2_ne ;
    Reset() ;
    return ntot ;
  } ;

  void IncVertex(void) { ++p2_nv ; } ;
  void IncEdge(void)   { ++p2_ne ; } ;
  void IncPoly(void)   { ++p2_np ; } ;

  unsigned nsV(void) const { return p2_nv ; }
  unsigned nsE(void) const { return p2_ne ; }
  unsigned nsP(void) const { return p2_np ; }

  void * psV(unsigned const i) const {
    if ( i >= p2_nv ) {
      cerr << endl << endl
           << "P2MESH in method: ``p2_vertex::vertex(" << i << ")''" << endl 
           << "Fatal error: local index out of range." << endl << endl ;
      exit(0) ;
    }
    return psV_[i] ;
  }
  
  void * psE(unsigned const i) const {
      if ( i >= p2_ne ) {
      cerr << endl << endl
           << "P2MESH in method: ``p2_vertex::edge(" << i << ")''" << endl 
           << "Fatal error: local index out of range." << endl << endl ;
      exit(0) ;
    }
    return psE_[i] ;
  }
  
  void * psP(unsigned const i) const {
      if ( i >= p2_np ) {
      cerr << endl << endl
           << "P2MESH in method: ``p2_vertex::poly(" << i << ")''" << endl 
           << "Fatal error: local index out of range." << endl << endl ;
      exit(0) ;
    }
    return psP_[i] ;
  }

  void IVertex(void * pv) { psV_[p2_nv++] = pv ; }
  void IEdge  (void * pe) { psE_[p2_ne++] = pe ; }
  void IPoly  (void * pp) { psP_[p2_np++] = pp ; }

public:
  p2_vertex_variant(void) : psV_(NULL), psE_(NULL), psP_(NULL),
                            p2_nv(0), p2_ne(0), p2_np(0) {}
} ;

template <>
class p2_vertex_variant<0> {
protected:

  void Reset(void) {}
  unsigned Assign(void **) { return 0 ; } ;

  void IncVertex(void) {} ;
  void IncEdge(void)   {} ;
  void IncPoly(void)   {} ;

  unsigned nsV(void) const { return 0 ; }
  unsigned nsE(void) const { return 0 ; }
  unsigned nsP(void) const { return 0 ; }

  void * psV(unsigned const i) const {
   cerr << endl << endl
        << "P2MESH in method: ``p2_vertex::vertex(" << i << ")''" << endl
        << "Fatal error: vertex list not defined." << endl << endl ;
   exit(0) ;
   return NULL ;
  }
  
  void * psE(unsigned const i) const {
    cerr << endl << endl
         << "P2MESH in method: ``p2_vertex::edge(" << i << ")''" << endl
         << "Fatal error: edge list not defined." << endl << endl ;
    exit(0) ;
    return NULL ;
  }

  void * psP(unsigned const i) const {
    cerr << endl << endl
         << "P2MESH in method: ``p2_vertex::poly(" << i << ")''" << endl
         << "Fatal error: polygon list not defined." << endl << endl ;
    exit(0) ;
    return NULL ;
  }

  void IVertex(void *) {}
  void IEdge  (void *) {}
  void IPoly  (void *) {}
} ;

/*
    #     #
    #     #  ######  #####    #####  ######  #    #
    #     #  #       #    #     #    #        #  #
    #     #  #####   #    #     #    #####     ##
     #   #   #       #####      #    #         ##
      # #    #       #   #      #    #        #  #
       #     ######  #    #     #    ######  #    #
*/

template <typename P2_COMMON>
class p2_vertex : public P2_COMMON,
                  private p2_vertex_variant<P2_COMMON::List> {           
public:
  typedef typename P2_COMMON::Real     Real ;
  typedef typename P2_COMMON::Integer  Integer ;
  typedef typename P2_COMMON::Unsigned Unsigned ;
  typedef typename P2_COMMON::P2V      P2V ;
  typedef typename P2_COMMON::P2E      P2E ;
  typedef typename P2_COMMON::P2P      P2P ;
  typedef typename P2_COMMON::P2M      P2M ;

  friend class p2_edge<P2_COMMON> ;
  friend class p2_poly<P2_COMMON> ;
  friend class p2_mesh<P2_COMMON> ;
  //friend class P2M ; commented to avoid a gcc 2.95 crash

private:
  Real p2_x, p2_y ;

  void InsertVertex(P2V * pv) { IVertex(static_cast<void*>(pv)) ; }
  void InsertEdge  (P2E * pe) { IEdge  (static_cast<void*>(pe)) ; }
  void InsertPoly  (P2P * pp) { IPoly  (static_cast<void*>(pp)) ; }

public:

  // propagate constant  
  static unsigned const Size = P2_COMMON::Size ;
  static bool     const List = P2_COMMON::List ;

  // this line issues an error if you try to use a copy constructor
  p2_vertex<P2_COMMON> const & operator = (p2_vertex<P2_COMMON> const &) {
    msg_error("p2_vertex::operator = ", "attempt to use copy constructor") ;
    return *this ;
  }

  // access vertex
  Real const & x(void) const { return p2_x ; }
  Real const & y(void) const { return p2_y ; }
  Real       & x(void)       { return p2_x ; }
  Real       & y(void)       { return p2_y ; }

  Unsigned n_vertex(void) const { return static_cast<Unsigned>(nsV()) ; }
  Unsigned n_edge  (void) const { return static_cast<Unsigned>(nsE()) ; }
  Unsigned n_poly  (void) const { return static_cast<Unsigned>(nsP()) ; }

  // access surrounding entities
  P2V const & vertex(Unsigned const nv) const { return * static_cast<P2V*>(psV(nv)) ; }
  P2V       & vertex(Unsigned const nv)       { return * static_cast<P2V*>(psV(nv)) ; }
  P2E const & edge  (Unsigned const ne) const { return * static_cast<P2E*>(psE(ne)) ; }
  P2E       & edge  (Unsigned const ne)       { return * static_cast<P2E*>(psE(ne)) ; }
  P2P const & poly  (Unsigned const np) const { return * static_cast<P2P*>(psP(np)) ; }
  P2P       & poly  (Unsigned const np)       { return * static_cast<P2P*>(psP(np)) ; }

  // getting local numbering
  Unsigned local_number(P2V const & rV) const ;
  Unsigned local_number(P2E const & rE) const ;
  Unsigned local_number(P2P const & rP) const ;
} ;

/*
    #######
    #        #####    ####   ######
    #        #    #  #    #  #
    #####    #    #  #       #####
    #        #    #  #  ###  #
    #        #    #  #    #  #
    #######  #####    ####   ######
 */

template <typename P2_COMMON>
class p2_edge : public P2_COMMON {
public:
  typedef typename P2_COMMON::Real     Real ;
  typedef typename P2_COMMON::Integer  Integer ;
  typedef typename P2_COMMON::Unsigned Unsigned ;
  typedef typename P2_COMMON::P2V      P2V ;
  typedef typename P2_COMMON::P2E      P2E ;
  typedef typename P2_COMMON::P2P      P2P ;
  typedef typename P2_COMMON::P2M      P2M ;

  friend class p2_vertex<P2_COMMON> ;
  friend class p2_poly  <P2_COMMON> ;
  friend class p2_mesh  <P2_COMMON> ;
  //friend class P2M ; commented to avoid a gcc 2.95 crash

private:
  P2P * p2_p[2] ; // pointer to the left(0) and right(1) poly
  P2V * p2_v[2] ; // pointer to the starting (0) and terminal (1) vertex

public:

  // propagate constant  
  static unsigned const Size = P2_COMMON::Size ;
  static bool     const List = P2_COMMON::List ;

  // this line issues an error if you try to use a copy constructor
  p2_edge<P2_COMMON> const & operator = (p2_edge<P2_COMMON> const &) {
    msg_error("p2_edge::operator = ", "attempt to use copy constructor");
    return *this ;
  }

  Unsigned n_vertex(void) const { return 2 ; }
  Unsigned n_edge  (void) const { return p2_p[1] == NULL ? Size - 1 : 2*(Size-1) ; }
  Unsigned n_poly  (void) const { return p2_p[1] == NULL ? 1 : 2 ; }

  // accessing surrounding entities
  P2V const & vertex(Unsigned const nv) const
  { check_range( nv, 2, "p2_edge::vertex" ) ; return *p2_v[nv] ; }
  P2V & vertex(Unsigned const nv)
  { check_range( nv, 2, "p2_edge::vertex" ) ; return *p2_v[nv] ; }
  
  P2E const & edge(Unsigned const ne) const ;
  P2E       & edge(Unsigned const ne) ;
  
  P2P const & poly(Unsigned const np) const
  { check_range( np, 2, "p2_edge::poly" ) ; return *p2_p[np] ; }
  P2P & poly(Unsigned const np)
  { check_range( np, 2, "p2_edge::poly" ) ; return *p2_p[np] ; }

  // return ``true'' if the polygon np exists
  bool        ok_poly(Unsigned const np) const { return p2_p[np] != NULL ; }

  // getting local numbering
  Unsigned local_number(P2V const & rV) const ;
  Unsigned local_number(P2E const & rE) const ;
  Unsigned local_number(P2P const & rP) const ;

  // accessing vertex coodinate
  Real const & x(Unsigned const nv) const
  { check_range( nv, 2, "p2_edge::x") ; return p2_v[nv] -> p2_x ; }
  Real const & y(Unsigned const nv) const
  { check_range( nv, 2, "p2_edge::y") ; return p2_v[nv] -> p2_y ; }

  // accessing midpoint
  Real xm(void) const { return 0.5 * (p2_v[0] -> p2_x + p2_v[1] -> p2_x) ; }
  Real ym(void) const { return 0.5 * (p2_v[0] -> p2_y + p2_v[1] -> p2_y) ; }

  // accessing point on edge at parametric coordinate
  Real xt(Real const & t) const
  { return p2_v[1] -> p2_x * t + p2_v[0] -> p2_x * (1-t) ; }
  Real yt(Real const & t) const
  { return p2_v[1] -> p2_y * t + p2_v[0] -> p2_y * (1-t) ; }

  // accessing normal, tangential and length
  Real nx(void) const { return p2_v[1] -> p2_y - p2_v[0] -> p2_y ; }
  Real ny(void) const { return p2_v[0] -> p2_x - p2_v[1] -> p2_x ; }
  Real tx(void) const { return p2_v[1] -> p2_x - p2_v[0] -> p2_x ; }
  Real ty(void) const { return p2_v[1] -> p2_y - p2_v[0] -> p2_y ; }
  Real length(void) const {
    Real a = nx() ;
    Real b = ny() ;
    return sqrt(a*a+b*b) ;
  }
} ;

/*
    ######
    #     #   ####   #      #   #   ####    ####   #    #
    #     #  #    #  #       # #   #    #  #    #  ##   #
    ######   #    #  #        #    #       #    #  # #  #
    #        #    #  #        #    #  ###  #    #  #  # #
    #        #    #  #        #    #    #  #    #  #   ##
    #         ####   ######   #     ####    ####   #    #
 */

template <typename P2_COMMON>
class p2_poly : public P2_COMMON {
public:
  typedef typename P2_COMMON::Real     Real ;
  typedef typename P2_COMMON::Integer  Integer ;
  typedef typename P2_COMMON::Unsigned Unsigned ;
  typedef typename P2_COMMON::P2V      P2V ;
  typedef typename P2_COMMON::P2E      P2E ;
  typedef typename P2_COMMON::P2P      P2P ;
  typedef typename P2_COMMON::P2M      P2M ;

  friend class p2_vertex<P2_COMMON> ;
  friend class p2_edge  <P2_COMMON> ;
  friend class p2_mesh  <P2_COMMON> ;
  //friend class P2M ; commented to avoid a gcc 2.95 crash

private:
  P2E * p2_e [P2_COMMON::Size] ; // pointers to the edges of the poly
  P2V * p2_v [P2_COMMON::Size] ; // pointers to the vertices of the poly

  Real const & p2_x(Unsigned const i) const { return p2_v[ i ] -> p2_x ; }
  Real const & p2_y(Unsigned const i) const { return p2_v[ i ] -> p2_y ; }
  
  static unsigned const p2_3 = P2_COMMON::Size-1 ; // to suppress warining

public:

  // propagate constant  
  static unsigned const Size = P2_COMMON::Size ;
  static bool     const List = P2_COMMON::List ;

  // this line issues an error if you try to use a copy constructor
  p2_poly<P2_COMMON> const & operator = (p2_poly<P2_COMMON> const &) {
    msg_error("p2_poly::operator = ", "attempt to use copy constructor");
    return *this ;
  }

  Unsigned n_vertex(void) const { return Size ; }
  Unsigned n_edge  (void) const { return Size ; }
  Unsigned n_poly  (void) const { return Size ; }

  // accessing surrounding entities
  P2V const & vertex(Unsigned const nv) const
  { check_range( nv, Size, "p2_poly::vertex") ; return *p2_v[nv] ; }
  P2V  & vertex(Unsigned const nv)
  { check_range( nv, Size, "p2_poly::vertex") ; return *p2_v[nv] ; }
  
  P2E const & edge(Unsigned const ne) const
  { check_range( ne, Size, "p2_poly::edge") ; return *p2_e[ne] ; }
  P2E  & edge(Unsigned const ne)
  { check_range( ne, Size, "p2_poly::edge") ; return *p2_e[ne] ; }
  
  P2P const & poly(Unsigned const np) const ;
  P2P       & poly(Unsigned const np) ;

  // return ``true'' if the polygon np exists
  bool        ok_poly(Unsigned const np) const ;

  // getting local numbering
  Unsigned local_number(P2V const & rV) const ;
  Unsigned local_number(P2E const & rE) const ;
  Unsigned local_number(P2P const & rP) const ;

  // attibute of the poly
  bool ok_oriented(Unsigned const ne) const {
    check_range( ne, Size, "p2_poly::ok_oriented" ) ;
    return p2_e[ne] -> p2_p[0] == this ;
  }

  // the coordinate of the vertex
  Real const & x(Unsigned const nv) const
  { check_range( nv, Size, "p2_poly::x" ) ; return p2_v[nv] -> p2_x ; }
  Real const & y(Unsigned const nv) const
  { check_range( nv, Size, "p2_poly::y" ) ; return p2_v[nv] -> p2_y ; }

  // the centroid of the poly
  Real xc(void) const {
    Real XX = p2_v[0] -> p2_x + p2_v[1] -> p2_x + p2_v[2] -> p2_x ;
    for ( Unsigned i = 3 ; i < Size ; ++i ) XX += p2_v[i] -> p2_x ;
    return XX / Size ;
  }
  Real yc(void) const {
    Real YY = p2_v[0] -> p2_y + p2_v[1] -> p2_y + p2_v[2] -> p2_y ;
    for ( Unsigned i = 3 ; i < Size ; ++i ) YY += p2_v[i] -> p2_y ;
    return YY / Size ;
  }

  // the area of the poly
  Real area(void) const {
    Real A = p2_v[Size-1] -> p2_x * p2_v[0] -> p2_y
           - p2_v[0] -> p2_x * p2_v[Size-1]-> p2_y ;
    for ( Unsigned i = 1 ; i < Size ; ++i)
      A += p2_v[i-1] -> p2_x * p2_v[i] -> p2_y
         - p2_v[i] -> p2_x * p2_v[i-1] -> p2_y ;
    return 0.5 * A ;
  }

  // outward normals
  Real nx(Unsigned const ne) const
  { return p2_v[(ne+1)%Size] -> p2_y - p2_v[ne%Size] -> p2_y ; }

  Real ny(Unsigned const ne) const
  { return p2_v[ne%Size] -> p2_x - p2_v[(ne+1)%Size] -> p2_x ; }

  // counterclockwise tangential
  Real tx(Unsigned const ne) const
  { return p2_v[(ne+1)%Size] -> p2_x - p2_v[ne%Size] -> p2_x ; } ;

  Real ty(Unsigned const ne) const
  { return p2_v[(ne+1)%Size] -> p2_y - p2_v[ne%Size] -> p2_y ; } ;

  // the midpoints of the edges
  Real xm(Unsigned const ne) const
  { check_range( ne, Size, "p2_poly::xm" ) ; return p2_e[ne] -> xm() ; }
  Real ym(Unsigned const ne) const
  { check_range( ne, Size, "p2_poly::ym" ) ; return p2_e[ne] -> ym() ; }

  // the point at coordinate t on the i edge
  Real xt(Unsigned const ne, Real const & t) const
  { return p2_v[(ne+1)%Size]->p2_x * t + p2_v[ne%Size]->p2_x * (1-t) ; }

  Real yt(Unsigned const ne, Real const & t) const
  { return p2_v[(ne+1)%Size]->p2_y * t + p2_v[ne%Size]->p2_y * (1-t) ; }

  Real length(Unsigned const ne) const
  { check_range( ne, Size, "p2_poly::length" ) ; return p2_e[ne] -> length() ; }

private:
  void jacobian_triangle(Real const & s, Real const & t, Real J[2][2]) const ;
  void jacobian_quad(Real const & s, Real const & t, Real J[2][2]) const ;
  void st_to_xy_triangle(Real const &, Real const &, Real &, Real &) const ;
  void st_to_xy_quad(Real const &, Real const &, Real &, Real &) const ;
  void xy_to_st_triangle(Real const &, Real const &, Real &, Real &) const ;
  void xy_to_st_quad(Real const &, Real const &, Real &, Real &) const ;
public:

  // special function for triangles & rectangles
  void jacobian(Real const & s, Real const & t, Real J[2][2]) const {
    if ( Size == 3 ) jacobian_triangle(s,t,J) ;
    else             jacobian_quad(s,t,J)     ;
  }

  void inverse_jacobian(Real const & s, Real const & t, Real J[2][2]) const;

  // transform s,t coordinate to x,y real coordinate
  void st_to_xy(Real const & s, Real const & t, Real & xx, Real & yy) const {
    if ( Size == 3 ) st_to_xy_triangle(s,t,xx,yy) ;
    else             st_to_xy_quad(s,t,xx,yy)     ;
  }

  // transform x, y coordinate to s, t local coordinate
  void xy_to_st(Real const & xx, Real const & yy, Real & s, Real & t) const {
    if ( Size == 3 ) xy_to_st_triangle(xx,yy,s,t) ;
    else             xy_to_st_quad(xx,yy,s,t)     ;
  }

  // transform grad s,t to grad x,y
  void grad_st_to_xy(Real const & s, Real const & t, Real const gst[2], Real gxy[2]) const ;

  // transform grad s,t to grad x,y
  void grad_xy_to_st(Real const & x, Real const & y, Real const gxy[2], Real gst[2]) const ;
} ;

/*
    #     #
    ##   ##  ######   ####   #    #
    # # # #  #       #       #    #
    #  #  #  #####    ####   ######
    #     #  #            #  #    #
    #     #  #       #    #  #    #
    #     #  ######   ####   #    #
*/

template <typename P2OBJ> class Iterator ;
template <typename P2OBJ> class CIterator ;

template <typename P2_COMMON>
class p2_mesh : public P2_COMMON {
public:
  typedef typename P2_COMMON::Real     Real ;
  typedef typename P2_COMMON::Integer  Integer ;
  typedef typename P2_COMMON::Unsigned Unsigned ;

  typedef typename P2_COMMON::P2V      P2V ;
  typedef typename P2_COMMON::P2E      P2E ;
  typedef typename P2_COMMON::P2P      P2P ;
  typedef typename P2_COMMON::P2M      P2M ;

  typedef typename P2_COMMON::Vmark    Vmark ;
  typedef typename P2_COMMON::Emark    Emark ;
  typedef typename P2_COMMON::Pmark    Pmark ;

  typedef vector<P2V> P2_LIST_V ;
  typedef vector<P2E> P2_LIST_E ;
  typedef vector<P2P> P2_LIST_P ;

  typedef void (*Mark_Vertex) (P2V &, Vmark const &) ;
  typedef void (*Mark_Edge)   (P2E &, Emark const &) ;
  typedef void (*Mark_Poly)   (P2P &, Pmark const &) ;
  typedef void (*Shape_Fun)   (Real const &, Real const &, Real &, Real &) ;

  friend class p2_vertex<P2_COMMON> ;
  friend class p2_edge  <P2_COMMON> ;
  friend class p2_poly  <P2_COMMON> ;

  //friend class P2M ; commented to avoid a gcc 2.95 crash

  typedef typename P2_LIST_V::iterator       vertex_iterator ;
  typedef typename P2_LIST_V::const_iterator vertex_const_iterator ;

  typedef typename P2_LIST_E::iterator       edge_iterator ;
  typedef typename P2_LIST_E::const_iterator edge_const_iterator ;

  typedef typename P2_LIST_P::iterator       poly_iterator ;
  typedef typename P2_LIST_P::const_iterator poly_const_iterator ;
  
  static unsigned const p2_3 = P2_COMMON::Size-1 ;
  // to suppress KCC warnings on subscript range

private:
  Unsigned p2_n_vertex  ; // total numbers of vertices
  Unsigned p2_n_bvertex ; // boundary vertices
  Unsigned p2_n_ivertex ; // internal vertices
  P2_LIST_V p2_vlist ;

  vertex_iterator p2_VA, p2_VB, p2_VC ;

  Unsigned p2_n_edge    ; // total numbers of edges
  Unsigned p2_n_bedge   ; // boundary edges
  Unsigned p2_n_iedge   ; // internal edges
  P2_LIST_E p2_elist ;

  edge_iterator p2_EA, p2_EB, p2_EC ;

  Unsigned p2_n_poly    ; // total numbers of polys
  Unsigned p2_n_bpoly   ; // boundary polys
  Unsigned p2_n_ipoly   ; // internal polys
  P2_LIST_P p2_plist ;

  poly_iterator p2_PA, p2_PB, p2_PC ;

  vector<void*> p2_vertex_list ;

  static Unsigned Mark(Unsigned const i,  Unsigned const j,
                       Unsigned const nx, Unsigned const ny) ;

  void map_meshQ(Unsigned const nx, Unsigned const ny, Mark_Poly mark_poly) ;

  void map_meshT0(Unsigned const nx, Unsigned const ny,
                  Mark_Edge mark_edge, Mark_Poly mark_poly) ;

  void map_meshT1(Unsigned const nx, Unsigned const ny,
                  Mark_Edge mark_edge, Mark_Poly mark_poly) ;

  void map_meshT2(Unsigned const nx, Unsigned const ny,
                  Mark_Vertex mark_vertex,
                  Mark_Edge   mark_edge,
                  Mark_Poly   mark_poly) ;

  void map_mesh_allocate(Unsigned const nx, Unsigned const ny,
                         Unsigned const kind) ;

  void map_mesh_internal(Unsigned const nx, Unsigned const ny,
                         Mark_Vertex mark_vertex,
                         Mark_Edge   mark_edge,
                         Mark_Poly   mark_poly,
                         Unsigned const kind) ;

  void JointEdges(void) ;
  void BuildEdges(void) ;
  void Reorder(void) ;
  void ReorderList(void) ;

public:
  // this line issue an error if you try to use a copy constructor
  p2_mesh<P2_COMMON> const & operator = (p2_mesh<P2_COMMON> const &) {
     msg_error("p2_mesh::operator = ", "attempt to use copy constructor") ;
     return *this ;
  }

  void bbox(Real & xmin, Real & ymin, Real & xmax, Real & ymax) const ;

  Unsigned const & n_vertex (void) const { return p2_n_vertex ; }
  Unsigned const & n_bvertex(void) const { return p2_n_bvertex ; }
  Unsigned const & n_ivertex(void) const { return p2_n_ivertex ; }
  Unsigned const & n_edge   (void) const { return p2_n_edge ; }
  Unsigned const & n_bedge  (void) const { return p2_n_bedge ; }
  Unsigned const & n_iedge  (void) const { return p2_n_iedge ; }
  Unsigned const & n_poly   (void) const { return p2_n_poly ; }
  Unsigned const & n_bpoly  (void) const { return p2_n_bpoly ; }
  Unsigned const & n_ipoly  (void) const { return p2_n_ipoly ; }

  P2V const & vertex(Unsigned const nv) const
  { check_grange( nv, p2_n_vertex, "p2_mesh::vertex") ; return p2_vlist[nv] ; }
  P2V & vertex(Unsigned const nv)
  { check_grange( nv, p2_n_vertex, "p2_mesh::vertex") ; return p2_vlist[nv] ; }

  P2E const & edge(Unsigned const ne) const
  { check_grange( ne, p2_n_edge, "p2_mesh::edge") ; return p2_elist[ne] ; }
  P2E & edge(Unsigned const ne)
  { check_grange( ne, p2_n_edge, "p2_mesh::edge") ; return p2_elist[ne] ; }

  P2P const & poly(Unsigned const np) const
  { check_grange( np, p2_n_poly, "p2_mesh::poly") ; return p2_plist[np] ; }
  P2P & poly(Unsigned const np)
  { check_grange( np, p2_n_poly, "p2_mesh::poly") ; return p2_plist[np] ; }

  Unsigned local_number(P2V const & rV) const
  { return (&rV) - (&p2_vlist[0]) ; }
  
  Unsigned local_number(P2E const & rE) const
  { return (&rE) - (&p2_elist[0]) ; }
  
  Unsigned local_number(P2P const & rP) const
  { return (&rP) - (&p2_plist[0]) ; }

  // iterators for vertex
  vertex_iterator       vertex_begin(void)       { return p2_VA ; }
  vertex_const_iterator vertex_begin(void) const { return p2_VA ; }
  vertex_iterator       vertex_end  (void)       { return p2_VC ; }
  vertex_const_iterator vertex_end  (void) const { return p2_VC ; }

  vertex_iterator       bvertex_begin(void)       { return p2_VA ; }
  vertex_const_iterator bvertex_begin(void) const { return p2_VA ; }
  vertex_iterator       bvertex_end  (void)       { return p2_VB ; }
  vertex_const_iterator bvertex_end  (void) const { return p2_VB ; }

  vertex_iterator       ivertex_begin(void)       { return p2_VB ; }
  vertex_const_iterator ivertex_begin(void) const { return p2_VB ; }
  vertex_iterator       ivertex_end  (void)       { return p2_VC ; }
  vertex_const_iterator ivertex_end  (void) const { return p2_VC ; }

  // iterators for edge
  edge_iterator       edge_begin(void)        { return p2_EA ; }
  edge_const_iterator edge_begin(void) const  { return p2_EA ; }
  edge_iterator       edge_end  (void)        { return p2_EC ; }
  edge_const_iterator edge_end  (void) const  { return p2_EC ; }

  edge_iterator       bedge_begin(void)       { return p2_EA ; }
  edge_const_iterator bedge_begin(void) const { return p2_EA ; }
  edge_iterator       bedge_end  (void)       { return p2_EB ; }
  edge_const_iterator bedge_end  (void) const { return p2_EB ; }

  edge_iterator       iedge_begin(void)       { return p2_EB ; }
  edge_const_iterator iedge_begin(void) const { return p2_EB ; }
  edge_iterator       iedge_end  (void)       { return p2_EC ; }
  edge_const_iterator iedge_end  (void) const { return p2_EC ; }

  // iterators for poly
  poly_iterator       poly_begin(void)       { return p2_PA ; }
  poly_const_iterator poly_begin(void) const { return p2_PA ; }
  poly_iterator       poly_end  (void)       { return p2_PC ; }
  poly_const_iterator poly_end  (void) const { return p2_PC ; }

  poly_iterator       bpoly_begin(void)       { return p2_PA ; }
  poly_const_iterator bpoly_begin(void) const { return p2_PA ; }
  poly_iterator       bpoly_end  (void)       { return p2_PB ; }
  poly_const_iterator bpoly_end  (void) const { return p2_PB ; }

  poly_iterator       ipoly_begin(void)       { return p2_PB ; }
  poly_const_iterator ipoly_begin(void) const { return p2_PB ; }
  poly_iterator       ipoly_end  (void)       { return p2_PC ; }
  poly_const_iterator ipoly_end  (void) const { return p2_PC ; }

private:
  vertex_iterator Get_base (P2V const * const) { return p2_VA ; }
  vertex_iterator Get_start(P2V const * const) { return p2_VB ; }
  vertex_iterator Get_end  (P2V const * const) { return p2_VC ; }

  edge_iterator   Get_base (P2E const * const) { return p2_EA ; }
  edge_iterator   Get_start(P2E const * const) { return p2_EB ; }
  edge_iterator   Get_end  (P2E const * const) { return p2_EC ; }

  poly_iterator   Get_base (P2P const * const) { return p2_PA ; }
  poly_iterator   Get_start(P2P const * const) { return p2_PB ; }
  poly_iterator   Get_end  (P2P const * const) { return p2_PC ; }

  vertex_const_iterator Get_base (P2V const * const) const { return p2_VA ; }
  vertex_const_iterator Get_start(P2V const * const) const { return p2_VB ; }
  vertex_const_iterator Get_end  (P2V const * const) const { return p2_VC ; }

  edge_const_iterator Get_base (P2E const * const) const { return p2_EA ; }
  edge_const_iterator Get_start(P2E const * const) const { return p2_EB ; }
  edge_const_iterator Get_end  (P2E const * const) const { return p2_EC ; }

  poly_const_iterator Get_base (P2P const * const) const { return p2_PA ; }
  poly_const_iterator Get_start(P2P const * const) const { return p2_PB ; }
  poly_const_iterator Get_end  (P2P const * const) const { return p2_PC ; }

public:
  void map_mesh(Shape_Fun shape_,
                Unsigned const nx, Unsigned const ny,
                Mark_Vertex mark_vertex,
                Mark_Edge   mark_edge,
                Mark_Poly   mark_poly,
                Unsigned const kind = 0) ;

  void tensor_mesh(Real const & xmin_, Real const & xmax_,
                   Real const & ymin_, Real const & ymax_,
                   Unsigned const nx,  Unsigned const   ny,
                   Mark_Vertex mark_vertex,
                   Mark_Edge   mark_edge,
                   Mark_Poly   mark_poly,
                   Unsigned const kind = 0) ;

  void std_tensor_mesh(Unsigned const nx, Unsigned const ny,
                       Mark_Vertex mark_vertex,
                       Mark_Edge   mark_edge,
                       Mark_Poly   mark_poly,
                       Unsigned const kind = 0) ;

  void build_mesh(Unsigned const nv,
                  Real     const *XY,
                  Vmark    const *mv,
                  Mark_Vertex mark_vertex,

                  Unsigned const ne,
                  Unsigned const *E,
                  Emark    const *me,
                  Mark_Edge  mark_edge,

                  Unsigned const np,
                  Unsigned const *P,
                  Pmark    const *mp,
                  Mark_Poly  mark_poly,
                  Unsigned const base = 0) ;

  void read_mesh(char const * const file_name,
                 Mark_Vertex mark_vertex,
                 Mark_Edge   mark_edge,
                 Mark_Poly   mark_poly,
                 Unsigned const base = 0) ;

  void read_map_mesh(char const * const file_name,
                     Mark_Vertex mark_vertex,
                     Mark_Edge   mark_edge,
                     Mark_Poly   mark_poly,
                     Unsigned const kind = 0) ;

  void report(ostream & s) const ;
  bool test_mesh(void) const ;

  void print(ostream & s, Unsigned const base = 0) const ;
  void print(ostream & s, P2E const & rE, Unsigned const base = 0) const ;
  void print(ostream & s, P2V const & rV, Unsigned const base = 0) const ;
  void print(ostream & s, P2P const & rP, Unsigned const base = 0) const ;

  friend class Iterator<P2V> ;
  friend class Iterator<P2E> ;
  friend class Iterator<P2P> ;
  friend class CIterator<P2V> ;
  friend class CIterator<P2E> ;
  friend class CIterator<P2P> ;
} ;

/*
      ###
       #      #####  ######  #####     ##     #####   ####   #####
       #        #    #       #    #   #  #      #    #    #  #    #
       #        #    #####   #    #  #    #     #    #    #  #    #
       #        #    #       #####   ######     #    #    #  #####
       #        #    #       #   #   #    #     #    #    #  #   #
      ###       #    ######  #    #  #    #     #     ####   #    #
*/

template <typename P2OBJ>
class Iterator {
public:
  typedef typename P2OBJ::P2M      P2M ;
  typedef typename P2OBJ::Unsigned Unsigned ;
  typedef typename P2OBJ::Real     Real ;

private:
  P2M * pMesh ;
  typename vector<P2OBJ>::iterator i ;
  typename vector<P2OBJ>::iterator start ;
  typename vector<P2OBJ>::iterator past_end ;

public:
  Iterator(void) {}

  Iterator(P2M & Mesh_, Unsigned const p2loop_ = 0 )
    { set_loop(Mesh_, p2loop_) ; }

  void set_loop(P2M & Mesh_, Unsigned const p2loop_ = 0) {
    pMesh = & Mesh_ ;
    P2OBJ * const cast = 0 ;
    switch ( p2loop_ ) {
    case 0: start    = pMesh -> Get_base(cast) ;
            past_end = pMesh -> Get_end (cast) ;
            break ;
    case 1: start    = pMesh -> Get_base (cast) ;
            past_end = pMesh -> Get_start(cast) ;
            break ;
    case 2: start    = pMesh -> Get_start(cast) ;
            past_end = pMesh -> Get_end  (cast) ;
            break ;
    default:
            P2M::test_ok(false, "Iterator::set_loop", p2loop_,
                                "bad loop specification") ;
    }
  }

  void begin(void)       { i = start ; }
  bool end_of_loop(void) { return i == past_end ; }

  Iterator<P2OBJ> const & operator ++ (void) { ++i ; return *this ; }

  // smart pointer stuff
  P2OBJ const * operator -> () const { return &*i ; }
  P2OBJ       * operator -> ()       { return &*i ; }
  P2OBJ const & operator *  () const { return *i ; }
  P2OBJ       & operator *  ()       { return *i ; }

} ;

/*
     #####    ###
    #     #    #      #####  ######  #####     ##     #####   ####   #####
    #          #        #    #       #    #   #  #      #    #    #  #    #
    #          #        #    #####   #    #  #    #     #    #    #  #    #
    #          #        #    #       #####   ######     #    #    #  #####
    #     #    #        #    #       #   #   #    #     #    #    #  #   #
     #####    ###       #    ######  #    #  #    #     #     ####   #    #
*/

template <typename P2OBJ>
class CIterator {
public:
  typedef typename P2OBJ::P2M      P2M ;
  typedef typename P2OBJ::Unsigned Unsigned ;
  typedef typename P2OBJ::Real     Real ;

private:
  P2M const * pMesh ;
  typename vector<P2OBJ>::const_iterator i ;
  typename vector<P2OBJ>::const_iterator start ;
  typename vector<P2OBJ>::const_iterator past_end ;

public:
  CIterator(void) {}

  CIterator(P2M const & Mesh_, Unsigned const p2loop_ = 0 )
    { set_loop(Mesh_, p2loop_) ; }

  void set_loop(P2M const & Mesh_, Unsigned const p2loop_ = 0) {
    pMesh = & Mesh_ ;
    P2OBJ const * const cast = 0 ;
    switch ( p2loop_ ) {
    case 0: start    = pMesh -> Get_base(cast) ;
            past_end = pMesh -> Get_end (cast) ;
            break ;
    case 1: start    = pMesh -> Get_base (cast) ;
            past_end = pMesh -> Get_start(cast) ;
            break ;
    case 2: start    = pMesh -> Get_start(cast) ;
            past_end = pMesh -> Get_end  (cast) ;
            break ;
    default:
            P2M::test_ok(false, "CIterator::set_loop", p2loop_,
                                "bad loop specification") ;
    }
  }

  void begin(void)       { i = start ; }
  bool end_of_loop(void) { return i == past_end ; }

  CIterator<P2OBJ> const & operator ++ (void) { ++i ; return *this ; }

  // smart pointer stuff
  P2OBJ const * operator -> () const { return &*i ; }
  P2OBJ const & operator *  () const { return *i ; }
} ;

# ifndef P2MESH_NO_FOREACH
  # define foreach(X) for ( X.begin() ; ! X.end_of_loop() ; ++X )
# endif

// IMPLEMENTATION PART

/*
    #     #
    #     #  ######  #####    #####  ######  #    #
    #     #  #       #    #     #    #        #  #
    #     #  #####   #    #     #    #####     ##
     #   #   #       #####      #    #         ##
      # #    #       #   #      #    #        #  #
       #     ######  #    #     #    ######  #    #
*/

template <typename P2_COMMON> inline
typename P2_COMMON::Unsigned
p2_vertex<P2_COMMON>::local_number(P2V const & rV) const {
  for ( Unsigned nv = 0 ; nv < n_vertex() ; ++nv )
    if ( &rV == & vertex(nv) ) return nv ;
  msg_error("p2_vertex::local_number(Vertex &)", "bad reference" ) ;
  return 0 ;
}

template <typename P2_COMMON> inline
typename P2_COMMON::Unsigned
p2_vertex<P2_COMMON>::local_number(P2E const & rE) const {
  for ( Unsigned ne = 0 ; ne < n_edge() ; ++ne )
    if ( &rE == & edge(ne) ) return ne ;
  msg_error("p2_vertex::local_number(Edge &)", "bad reference" ) ;
  return 0 ;
}

template <typename P2_COMMON> inline
typename P2_COMMON::Unsigned
p2_vertex<P2_COMMON>::local_number(P2P const & rP) const {
  for ( Unsigned np = 0 ; np < n_poly() ; ++np )
    if ( &rP == & poly(np) ) return np ;
  msg_error("p2_vertex::local_number(Poly &)", "bad reference" ) ;
  return 0 ;
}

/*
    #######
    #        #####    ####   ######
    #        #    #  #    #  #
    #####    #    #  #       #####
    #        #    #  #  ###  #
    #        #    #  #    #  #
    #######  #####    ####   ######
*/

template <typename P2_COMMON> inline
typename P2_COMMON::P2E const &
p2_edge<P2_COMMON>::edge(Unsigned const ne) const {
  check_range( ne, n_edge(), "p2_edge::edge") ;
  Unsigned ipoly = ne < Size - 1 ? 0 : 1 ;
  P2P const & P = poly(ipoly) ;
  Unsigned le = (P.local_number(static_cast<P2E const&>(*this)) + ne + 1 + ipoly ) % Size ;
  return P.edge( le ) ;
}

template <typename P2_COMMON> inline
typename P2_COMMON::P2E &
p2_edge<P2_COMMON>::edge(Unsigned const ne) {
  check_range( ne, n_edge(), "p2_edge::edge") ;
  Unsigned ipoly = ne < Size - 1 ? 0 : 1 ;
  P2P & P = poly(ipoly) ;
  Unsigned le = (P.local_number(static_cast<P2E&>(*this)) + ne + 1 + ipoly) % Size ;
  return P.edge( le ) ;
}

// getting local numbering
template <typename P2_COMMON> inline
typename P2_COMMON::Unsigned
p2_edge<P2_COMMON>::local_number(P2V const & rV) const {
  if ( &rV == p2_v[0] ) return 0 ;
  if ( &rV == p2_v[1] ) return 1 ;
  msg_error("p2_edge::local_number(vertex &)", "bad reference" ) ;
  return 0 ;
}

template <typename P2_COMMON> inline
typename P2_COMMON::Unsigned
p2_edge<P2_COMMON>::local_number(P2E const & rE) const {
  for ( Unsigned i = 0 ; i < n_edge() ; ++i )
    if ( &rE == & edge(i) ) return i ;
  msg_error("p2_edge::local_number(edge &)", "bad reference") ;
  return 0 ;
}

template <typename P2_COMMON> inline
typename P2_COMMON::Unsigned
p2_edge<P2_COMMON>::local_number(P2P const & rP) const {
  if ( p2_p[0] == &rP ) return 0 ;
  if ( p2_p[1] == &rP ) return 1 ;
  msg_error("p2_edge::local_number(poly &)", "bad reference") ;
  return 0 ;
}

/*
    ######
    #     #   ####   #      #   #   ####    ####   #    #
    #     #  #    #  #       # #   #    #  #    #  ##   #
    ######   #    #  #        #    #       #    #  # #  #
    #        #    #  #        #    #  ###  #    #  #  # #
    #        #    #  #        #    #    #  #    #  #   ##
    #         ####   ######   #     ####    ####   #    #
*/

template <typename P2_COMMON> inline
typename P2_COMMON::P2P const &
p2_poly<P2_COMMON>::poly(Unsigned const np) const {
  check_range( np, Size, "p2_poly::poly" ) ;
  register P2E const * pE = p2_e[np] ;
  register P2P const * pP = this == pE->p2_p[0] ? pE->p2_p[1] : pE->p2_p[0] ;
  test_ok( pP != NULL, "p2_poly::poly", np,"invalid polygon access") ;
  return *pP ;
}

template <typename P2_COMMON> inline
typename P2_COMMON::P2P &
p2_poly<P2_COMMON>::poly(Unsigned const np) {
  check_range( np, Size, "p2_poly::poly" ) ;
  register P2E * pE = p2_e[np] ;
  register P2P * pP = this == pE -> p2_p[0] ? pE -> p2_p[1] : pE -> p2_p[0] ;
  test_ok( pP != NULL, "p2_poly::poly", np,"invalid polygon access") ;
  return *pP ;
}

template <typename P2_COMMON> inline
bool
p2_poly<P2_COMMON>::ok_poly(Unsigned const np) const {
  check_range( np, Size, "p2_poly::ok_poly" ) ;
  register P2E const * pE = p2_e[np] ;
  register P2P const * pP = this != pE->p2_p[0] ? pE->p2_p[0] : pE->p2_p[1] ;
  return pP != NULL ;
}

// getting local numbering
template <typename P2_COMMON> inline
typename P2_COMMON::Unsigned
p2_poly<P2_COMMON>::local_number(P2V const & rV) const {
  for ( Unsigned i = 0 ; i < Size ; ++i )
    if ( &rV == p2_v[i] ) return i ;
  msg_error("p2_poly::local_number(vertex &)", "bad reference" ) ;
  return 0 ;
}

template <typename P2_COMMON> inline
typename P2_COMMON::Unsigned
p2_poly<P2_COMMON>::local_number(P2E const & rE) const {
  for ( Unsigned i = 0 ; i < Size ; ++i )
    if ( &rE == p2_e[i] ) return i ;
  msg_error("p2_poly::local_number(edge &)", "bad reference") ;
  return 0 ;
}

template <typename P2_COMMON>
typename P2_COMMON::Unsigned
p2_poly<P2_COMMON>::local_number(P2P const & rP) const {
  for ( Unsigned np = 0 ; np < Size ; ++np )
    if ( ok_poly(np) )
      if ( &rP == &poly(np) ) return np ;
  msg_error("p2_poly::local_number(poly &)", "bad reference") ;
  return 0 ;
}

// transform s,t coordinate to x,y real coordinate
template <typename P2_COMMON> inline
void
p2_poly<P2_COMMON>::st_to_xy_triangle(Real const & s, Real const & t,
                                      Real & xx, Real & yy) const {
  Real r = 1-s-t ;
  xx = r * p2_x(0) + s * p2_x(1) + t * p2_x(2) ;
  yy = r * p2_y(0) + s * p2_y(1) + t * p2_y(2) ;
}

// transform s,t coordinate to x,y real coordinate
template <typename P2_COMMON> inline
void
p2_poly<P2_COMMON>::st_to_xy_quad(Real const & s, Real const & t,
                                  Real & xx, Real & yy) const {
  Real st = s*t ;
  Real a0 = 1-s-t+st ;
  Real a1 = 1+s-t-st ;
  Real a2 = 1+s+t+st ;
  Real a3 = 1-s+t-st ;
  xx = 0.25 * ( a0*p2_x(0) + a1*p2_x(1) + a2*p2_x(2) + a3*p2_x(p2_3) ) ;
  yy = 0.25 * ( a0*p2_y(0) + a1*p2_y(1) + a2*p2_y(2) + a3*p2_y(p2_3) ) ;
}

// transform x, y coordinate to s, t local coordinate
template <typename P2_COMMON> inline
void
p2_poly<P2_COMMON>::xy_to_st_triangle(Real const & xin, Real const & yin,
                                      Real & s, Real & t) const {
  Real xx  = xin - p2_x(0) ;
  Real yy  = yin - p2_y(0) ;
  Real v1x = p2_x(1) - p2_x(0) ;
  Real v1y = p2_y(1) - p2_y(0) ;
  Real v2x = p2_x(2) - p2_x(0) ;
  Real v2y = p2_y(2) - p2_y(0) ;

  Real den = v1x * v2y - v1y * v2x ;
  s = ( xx * v2y - yy * v2x ) / den ;
  t = ( v1x * yy - v1y * xx ) / den ;

}

// transform x, y coordinate to s, t local coordinate
template <typename P2_COMMON>
void
p2_poly<P2_COMMON>::xy_to_st_quad(Real const & xx, Real const & yy,
                                  Real & s, Real & t) const {
  
  Real xx2 = 2*xx ;
  Real yy2 = 2*yy ;
                                  
  Real c0x = p2_x(0) + p2_x(1) - xx2 ;
  Real c0y = p2_y(0) + p2_y(1) - yy2 ;

  Real c1x = p2_x(2) + p2_x(p2_3) - xx2 ;
  Real c1y = p2_y(2) + p2_y(p2_3) - yy2 ;

  Real c2x = p2_x(1) - p2_x(0) ;
  Real c2y = p2_y(1) - p2_y(0) ;

  Real c3x = p2_x(2) - p2_x(p2_3) ;
  Real c3y = p2_y(2) - p2_y(p2_3) ;

  Real d0x = p2_x(0) + p2_x(p2_3) - xx2 ;
  Real d0y = p2_y(0) + p2_y(p2_3) - yy2 ;

  Real d1x = p2_x(1) + p2_x(2) - xx2 ;
  Real d1y = p2_y(1) + p2_y(2) - yy2 ;

  Real d2x = p2_x(p2_3) - p2_x(0) ;
  Real d2y = p2_y(p2_3) - p2_y(0) ;

  Real d3x = p2_x(2) - p2_x(1) ;
  Real d3y = p2_y(2) - p2_y(1) ;

  Real as = c2x * c3y - c2y * c3x ;
  Real bs = c2y * c1x - c2x * c1y + c0y * c3x - c0x * c3y ;
  Real cs = c0x * c1y - c0y * c1x ;
  Real ds = sqrt( bs * bs - 4 * as * cs ) ;

  Real at = d2x * d3y - d2y * d3x ;
  Real bt = d2y * d1x - d2x * d1y +  d0y * d3x - d0x * d3y ;
  Real ct = d0x * d1y - d0y * d1x ;
  Real dt = sqrt( bt * bt - 4 * at * ct ) ;

  if ( abs(as) > p2_mesh_namespace::big_epsi * abs(bs) ) {
    Real s1 = bs - ds ;
    Real s2 = bs + ds ;
    s = ( abs(s1) > abs(s2) ? s2 : s1 ) / (2*as) ;
  } else {
    s = cs / bs ;
  }
  
  if ( abs(at) > p2_mesh_namespace::big_epsi * abs(bt) ) {
    Real t1 = bt - dt ;
    Real t2 = bt + dt ;
    t = ( abs(t1) > abs(t2) ? t2 : t1 ) / (2*at) ;
  } else {
    t = ct / bt ;
  }
}

template <typename P2_COMMON> inline
void
p2_poly<P2_COMMON>::jacobian_triangle(Real const &, Real const &,
                                      Real J[2][2]) const {
  J[0][0] = p2_v[1] -> p2_x - p2_v[0] -> p2_x ;
  J[0][1] = p2_v[2] -> p2_x - p2_v[0] -> p2_x ;
  J[1][0] = p2_v[1] -> p2_y - p2_v[0] -> p2_y ;
  J[1][1] = p2_v[2] -> p2_y - p2_v[0] -> p2_y ;
}

template <typename P2_COMMON> inline
void
p2_poly<P2_COMMON>::jacobian_quad(Real const & s, Real const & t,
                                  Real J[2][2]) const {
  Real dx = p2_x(0) - p2_x(1) + p2_x(2) - p2_x(p2_3) ;
  Real dy = p2_y(0) - p2_y(1) + p2_y(2) - p2_y(p2_3) ;
  Real ss = 0.5*(s + 1) ;
  Real tt = 0.5*(t + 1) ;
  J[0][0] = 0.5*(p2_x(1)    - p2_x(0) + tt * dx) ;
  J[0][1] = 0.5*(p2_x(p2_3) - p2_x(0) + ss * dx) ;
  J[1][0] = 0.5*(p2_y(1)    - p2_y(0) + tt * dy) ;
  J[1][1] = 0.5*(p2_y(p2_3) - p2_y(0) + ss * dy) ;
}

template <typename P2_COMMON> inline
void
p2_poly<P2_COMMON>::inverse_jacobian(Real const & s, Real const & t,
                                     Real InvJ[2][2]) const {
  Real J[2][2] ;
  jacobian(s,t,J) ;
  Real detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0] ;
  InvJ[0][0] =   J[1][1] / detJ ;
  InvJ[0][1] = - J[0][1] / detJ ;
  InvJ[1][0] = - J[1][0] / detJ ;
  InvJ[1][1] =   J[0][0] / detJ ;
}

// transform grad s,t to grad x,y
template <typename P2_COMMON> inline
void
p2_poly<P2_COMMON>::grad_st_to_xy(Real const & s,
                                  Real const & t,
                                  Real const   gst[2],
                                  Real         gxy[2]) const {
  Real InvJ[2][2] ;
  inverse_jacobian(s,t,InvJ) ;
  gxy[0] = InvJ[0][0] * gst[0] + InvJ[1][0] * gst[1] ;
  gxy[1] = InvJ[0][1] * gst[0] + InvJ[1][1] * gst[1] ;
}

// transform grad s,t to grad x,y
template <typename P2_COMMON> inline
void
p2_poly<P2_COMMON>::grad_xy_to_st(Real const & xx,
                                  Real const & yy,
                                  Real const   gxy[2],
                                  Real         gst[2]) const {

  Real s, t, J[2][2] ;
  xy_to_st(xx, yy, s, t) ;
  jacobian(s,t,J) ;
  gst[0] = J[0][0] * gxy[0] + J[1][0] * gxy[1] ;
  gst[1] = J[0][1] * gxy[0] + J[1][1] * gxy[1] ;
}

/*
        #     #
        ##   ##  ######   ####   #    #
        # # # #  #       #       #    #
        #  #  #  #####    ####   ######
        #     #  #            #  #    #
        #     #  #       #    #  #    #
        #     #  ######   ####   #    #
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//
//  8 ----- 3 ------7
//  |               |
//  4       0       2
//  |               |
//  5 ----- 1 ----- 6
//
template <typename P2_COMMON>
typename P2_COMMON::Unsigned
p2_mesh<P2_COMMON>::Mark(Unsigned const i,  Unsigned const j,
                         Unsigned const nx, Unsigned const ny) {
  static Unsigned const mark[] = { 5, 1, 6, 4, 0, 2, 8, 3, 7 } ;
  Unsigned mi = 1 ;
  if ( i == 0 )       mi = 0 ;
  else if ( i == nx ) mi = 2 ;

  Unsigned mj = 1 ;
  if ( j == 0 )       mj = 0 ;
  else if ( j == ny ) mj = 2 ;

  return mark[ mi + 3 * mj ] ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::BuildEdges() {
  P2MESH_MSG("enter p2_mesh::BuildEdges()") ;
  
  typedef p2_mesh_namespace::cmp_pedge<P2E*> cmp_pedge ;

  p2_n_edge = 0 ;
  p2_elist.resize( Size * p2_n_poly + 1 ) ;
  vector<P2E*> ptre_list ;

  edge_iterator ie  = p2_elist.begin() ;
  typename vector<P2E*>::iterator ipe ;

  // build edges
  for ( poly_iterator ip = p2_plist.begin() ;
        ip != p2_plist.end() ;
        ++ip ) {

    P2V * pb = ip -> p2_v[0] ;

    for ( Unsigned ne = 0 ; ne < Size ; ++ne ) {
      P2V * pa = pb ; pb = ip -> p2_v[(ne+1) % Size] ;

      // setup edge
      ie -> p2_v[0] = pa   ; ie -> p2_v[1] = pb ;
      ie -> p2_p[0] = NULL ; ie -> p2_p[1] = NULL ;

      // search edge in the list
      P2E * pE = &*ie ;
      bool found = p2_mesh_namespace::binary_search
         (ptre_list.begin(), ptre_list.end(), pE, cmp_pedge(), ipe) ;

      if ( found ) {
        pE = *ipe ; // found
      } else { // not found (insert vertex in the list)
        ptre_list.insert( ipe, pE ) ;
        ++ie ;
        ++p2_n_edge ;
      }

      P2P ** ppP = & pE -> p2_p[found ? 1 : 0] ;
      if ( *ppP != NULL )
         msg_error("p2_mesh::BuildEdges()", 
                   "an edge is referenced twice from the same side" ) ;
      *ppP = &*ip ;
      ip -> p2_e[ne] = pE ;
    }
  }

  p2_elist.resize( p2_n_edge ) ;

  P2MESH_MSG("exit p2_mesh::BuildEdges()") ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::JointEdges() {
  
  typedef p2_mesh_namespace::cmp_pedge<P2E*> cmp_pedge ;
  
  P2MESH_MSG("enter p2_mesh::JointEdges()") ;

  vector<P2E*> elist(p2_n_edge) ;
  typename vector<P2E*>::iterator ipe = elist.begin() ;
  for ( edge_iterator ie = p2_elist.begin() ;
        ie != p2_elist.end() ;
        ++ipe, ++ie )
    *ipe = &*ie ;

  P2MESH_MSG("p2_mesh::JointEdges() ...sort edge list") ;
  sort(elist.begin(), elist.end(), cmp_pedge() ) ;

  P2E e ;
  P2E * pe = &e;

  for ( poly_iterator ip = p2_plist.begin() ;
        ip != p2_plist.end() ;
        ++ip ) {

    // check vertex pointer
    for ( Unsigned nev = 0 ; nev < Size ; ++nev )
      if ( ip -> p2_v[nev] == NULL )
        msg_error( "p2_mesh::JointEdges()", "incomplete polygon found" ) ;

    // search edges to join
    for ( Unsigned ne = 0 ; ne < Size ; ++ne ) {
      Unsigned ne1 = (ne+1) % Size ;
      e.p2_v[0] = ip -> p2_v[ne] ;
      e.p2_v[1] = ip -> p2_v[ne1] ;

      bool found = p2_mesh_namespace::binary_search
        (elist.begin(), elist.end(), pe, cmp_pedge(), ipe) ;

      if ( !found ) {
        e.p2_v[0] = ip -> p2_v[ne1] ;
        e.p2_v[1] = ip -> p2_v[ne] ;
        bool ok = p2_mesh_namespace::binary_search
          (elist.begin(), elist.end(), pe, cmp_pedge(), ipe) ;
        if ( !ok )
          msg_error( "p2_mesh::JointEdges()",
                     "try to build a polygon with a not existing edge" ) ;
      }

      // join edge to the polygon
      ip -> p2_e[ ne ] = *ipe ;

      P2P ** ppP = & (*ipe) -> p2_p[ found ? 0 : 1 ] ;
      if ( *ppP != NULL )
        msg_error( "p2_mesh::JointEdges()",
                   "try to assign a polygon to an already assigned edge side" );
      *ppP = &*ip ;

    }
  }
  P2MESH_MSG("exit p2_mesh::JointEdges()") ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::Reorder() {

  P2MESH_MSG("enter p2_mesh::Reorder()") ;

  edge_iterator ie ;
  poly_iterator ip ;

  P2MESH_MSG("p2_mesh::Reorder() ...orient and count boundary edges") ;
  // reverse boundary edge when necessary
  // and count boundary edges
  p2_n_bedge = 0 ;
  for ( ie = p2_elist.begin() ; ie != p2_elist.end() ; ++ie ) {
    if ( ie -> p2_p[0] == NULL ) {
      swap( ie -> p2_v[0], ie -> p2_v[1] ) ;
      swap( ie -> p2_p[0], ie -> p2_p[1] ) ;
    }
    if ( ie -> p2_p[1] == NULL ) ++p2_n_bedge ;
    if ( ie -> p2_p[0] == NULL )
      msg_error( "p2_mesh::Reorder()", "isolated edge found" ) ;
    
    if ( ie -> p2_v[0] == NULL || ie -> p2_v[1] == NULL )
      msg_error( "p2_mesh::Reorder()", "incomplete edge found" ) ;
  }
    
  P2MESH_MSG("p2_mesh::Reorder() ...numbering vertices, edges, triangles") ;

  // auxiliary vectors
  vector<Integer> vlist(p2_n_vertex) ;
  vector<Integer> elist(p2_n_edge) ;
  vector<Integer> plist(p2_n_poly) ;

  vector<P2E*> belist(p2_n_bedge) ;
  typename vector<P2E*>::iterator pbe, pbe1 ;

  fill( vlist.begin(), vlist.end(), -1 ) ;
  fill( elist.begin(), elist.end(), -1 ) ;
  fill( plist.begin(), plist.end(), -1 ) ;

  pbe = belist.begin() ;
  for ( ie = p2_elist.begin() ; ie != p2_elist.end() ; ++ie )
    if ( ie -> p2_p[1] == NULL ) *pbe++ = &*ie ;

  p2_n_bvertex = p2_n_bedge = p2_n_bpoly = 0 ;
  // number the boundary edge, boundary vertex and boundary polygons
  for ( pbe = belist.begin() ; pbe != belist.end() ; ++pbe ) {
    if ( *pbe == NULL ) continue ; // visited

    P2E * pE = *pbe ;
    *pbe = NULL ; // mark as visited

    P2V * pA = pE -> p2_v[0] ;
    P2V * pV = pA ;

  ok_found:

    P2P * pP = pE -> p2_p[0] ;

    Integer & vmark = vlist[ local_number(*pV) ];
    Integer & emark = elist[ local_number(*pE) ];
    Integer & pmark = plist[ local_number(*pP) ];

    if ( emark == -1 ) emark = p2_n_bedge++ ; // add edge
    else msg_error("p2_mesh::Reorder()", "corrupted boundary") ;

    if ( vmark == -1 ) vmark = p2_n_bvertex++ ; // add vertex
    if ( pmark == -1 ) pmark = p2_n_bpoly++ ; // add poly

    pV = pE -> p2_v[1] ;

    if ( pV == pA ) continue ; // next boundary

    // find the next contiguous edge
    for ( pbe1 = belist.begin() ; pbe1 != belist.end() ; ++pbe1 ) {
      if ( *pbe1 == NULL ) continue ;
      if ( pV == (*pbe1) -> p2_v[0] ) {
        pE   = *pbe1 ;
        *pbe1 = NULL ;
        goto ok_found ;
      }
    }
    msg_error("p2_mesh::Reorder()", "open boundary" ) ;
  }

  // number the internal vertices, edges, polygons
  typename vector<Integer>::iterator ii ;
  Unsigned j ;
  for ( j = p2_n_bvertex, ii = vlist.begin() ; ii != vlist.end() ; ++ii )
    if ( *ii == -1 ) *ii = j++ ;

  for ( j = p2_n_bedge, ii = elist.begin() ; ii != elist.end() ; ++ii )
    if ( *ii == -1 ) *ii = j++ ;

  for ( j = p2_n_bpoly, ii = plist.begin() ; ii != plist.end() ; ++ii )
    if ( *ii == -1 ) *ii = j++ ;

  P2MESH_MSG("p2_mesh::Reorder() ...ordering pointers internal to edges") ;
  for ( ip = p2_plist.begin() ; ip != p2_plist.end() ; ++ip ) {
    for ( Unsigned nve = 0 ; nve < Size ; ++nve ) {
      ip -> p2_v[ nve ] = &p2_vlist[vlist[local_number( * ip -> p2_v[nve] )]];
      ip -> p2_e[ nve ] = &p2_elist[elist[local_number( * ip -> p2_e[nve] )]];
    }
  }

  P2MESH_MSG("p2_mesh::Reorder() ...ordering pointers internal to polygons") ;
  for ( ie = p2_elist.begin() ; ie != p2_elist.end() ; ++ie ) {

    ie -> p2_v[0] = &p2_vlist[ vlist[ local_number( * ie -> p2_v[0] ) ] ] ;
    ie -> p2_v[1] = &p2_vlist[ vlist[ local_number( * ie -> p2_v[1] ) ] ] ;

    if ( ie -> p2_p[0] != NULL )
      ie -> p2_p[0] = &p2_plist[plist[ local_number( * ie -> p2_p[0] ) ] ] ;

    if ( ie -> p2_p[1] != NULL )
      ie -> p2_p[1] = &p2_plist[plist[ local_number( * ie -> p2_p[1] ) ] ] ;

  }

  P2MESH_MSG("p2_mesh::Reorder() ...ordering vertices") ;
  Unsigned i, jbf ;
  for ( i = 0 ; i < p2_n_vertex ; ++i ) {
    if ( vlist[i] == -1 ) continue ;
    P2V * vi = &p2_vlist[i], buffer ;
    uninitialized_copy( vi, vi + 1, &buffer) ;
    j = i ;
    do {
      jbf = vlist[j] ; vlist[j] = -1 ; j = jbf;
      P2V * vj = &p2_vlist[j], buffer1 ;
      uninitialized_copy( vj, vj + 1, &buffer1) ;
      uninitialized_copy( &buffer, (&buffer) + 1, vj) ;
      uninitialized_copy( &buffer1, (&buffer1) + 1, &buffer) ;
    } while ( i != j ) ;
  }
  
  P2MESH_MSG("p2_mesh::Reorder() ...ordering edges") ;
  for ( i = 0 ; i < p2_n_edge ; ++i ) {
    if ( elist[i] == -1 ) continue ;
    P2E * ei = &p2_elist[i], buffer ;
    uninitialized_copy( ei, ei + 1, &buffer) ;
    j = i ;
    do {
      jbf = elist[j] ; elist[j] = -1 ; j = jbf ;
      P2E * ej = &p2_elist[j], buffer1 ;
      uninitialized_copy( ej, ej + 1, &buffer1) ;
      uninitialized_copy( &buffer, (&buffer) + 1, ej) ;
      uninitialized_copy( &buffer1, (&buffer1) + 1, &buffer) ;
    } while ( i != j ) ;
  }
  
  P2MESH_MSG("p2_mesh::Reorder() ...ordering polygons") ;
  for ( i = 0 ; i < p2_n_poly ; ++i ) {
    if ( plist[i] == -1 ) continue ;
    P2P * pi = &p2_plist[i] ;
    P2P buffer ;
    uninitialized_copy( pi, pi + 1, &buffer) ;
    j = i ;
    do {
      jbf = plist[j] ; plist[j] = -1 ; j = jbf ;
      P2P * pj = &p2_plist[j] ;
      P2P buffer1 ;
      uninitialized_copy( pj, pj + 1, &buffer1) ;
      uninitialized_copy( &buffer, (&buffer) + 1, pj) ;
      uninitialized_copy( &buffer1, (&buffer1) + 1, &buffer) ;
    } while ( i != j ) ;
  }

  // setup internal counting
  p2_n_ivertex = p2_n_vertex - p2_n_bvertex ;
  p2_n_iedge   = p2_n_edge   - p2_n_bedge ;
  p2_n_ipoly   = p2_n_poly   - p2_n_bpoly ;

  // setup for iterators
  p2_VA = p2_vlist.begin() ;
  p2_VB = p2_vlist.begin() + p2_n_bvertex ;
  p2_VC = p2_vlist.end() ;

  p2_EA = p2_elist.begin() ;
  p2_EB = p2_elist.begin() + p2_n_bedge ;
  p2_EC = p2_elist.end() ;

  p2_PA = p2_plist.begin() ;
  p2_PB = p2_plist.begin() + p2_n_bpoly ;
  p2_PC = p2_plist.end() ;

  if ( List ) ReorderList() ;

  P2MESH_MSG("exit p2_mesh::Reorder()") ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::ReorderList() {

  P2MESH_MSG("enter p2_mesh::ReorderList()") ;
  vertex_iterator iv ;
  edge_iterator   ie ;
  poly_iterator   ip ;

  P2MESH_MSG("p2_mesh::ReorderList() ...setup vertices list") ;
  for ( iv = p2_vlist.begin() ; iv != p2_vlist.end() ; ++iv )
    iv -> Reset() ;

  P2MESH_MSG("p2_mesh::ReorderList() ...count vertex-edge list") ;
  unsigned long size_vertex_list = 0 ;
  for ( ie = p2_elist.begin() ; ie != p2_elist.end() ; ++ie ) {
    ie -> p2_v[0] -> IncVertex() ;
    ie -> p2_v[0] -> IncEdge()   ;
    ie -> p2_v[1] -> IncVertex() ;
    ie -> p2_v[1] -> IncEdge()   ;
    size_vertex_list += 4 ;
  }
  P2MESH_MSG("p2_mesh::ReorderList() ...count polygons list") ;
  for ( ip = p2_plist.begin() ; ip != p2_plist.end() ; ++ip ) {
    for ( Unsigned nv = 0 ; nv < Size ; ++nv )
      ip -> p2_v[nv] -> IncPoly() ;
    size_vertex_list += unsigned(Size) ;
  }

  p2_vertex_list.resize(size_vertex_list) ;
  void ** ptr = & p2_vertex_list[0] ;

  P2MESH_MSG("p2_mesh::ReorderList() ...assign memory") ;
  for ( iv = p2_vlist.begin() ; iv != p2_vlist.end() ; ++iv )
    ptr += iv -> Assign(ptr) ;

  P2MESH_MSG("p2_mesh::ReorderList() ...insert edges") ;
  for ( ie = p2_elist.begin() ; ie != p2_elist.end() ; ++ie ) {
    ie -> p2_v[0] -> InsertEdge(&*ie) ;
    ie -> p2_v[1] -> InsertEdge(&*ie) ;
  }

  P2MESH_MSG("p2_mesh::ReorderList() ...insert vertices") ;
  for ( iv = p2_vlist.begin() ; iv != p2_vlist.end() ; ++iv ) {
    for ( Unsigned ne = 0 ; ne < iv -> n_edge() ; ++ne ) {
       P2E & E = iv -> edge(ne) ;
       iv -> InsertVertex( &*iv == E.p2_v[0] ? E.p2_v[1] : E.p2_v[0] ) ;
    }
  }

  P2MESH_MSG("p2_mesh::ReorderList() ...insert polygons") ;
  for ( ip = p2_plist.begin() ; ip != p2_plist.end() ; ++ip ) {
    for ( Unsigned nv = 0 ; nv < Size ; ++nv )
       ip -> p2_v[nv] -> InsertPoly(&*ip) ;
  }

  P2MESH_MSG("exit p2_mesh::ReorderList()") ;
}

//
// typical triangular mesh
// +---+---+---+   +---+---+---+   +-----+--+--+
// | \ | \ | \ |   | / | / | / |   | \ / | \ / |
// +---+---+---+   +---+---+---+   +  +  +  +  +
// | \ | \ | \ |   | / | / | / |   | / \ | / \ |
// +---+---+---+   +---+---+---+   +--+--+--+--+
//   kind = 0        kind = 1        kind = 2
//
// typical quadrilateral mesh
// +---+---+---+
// |   |   |   |
// +---+---+---+
// |   |   |   |
// +---+---+---+

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::tensor_mesh(Real const & xmin, Real const & xmax,
                                Real const & ymin, Real const & ymax,
                                Unsigned const nx, Unsigned const ny,
                                Mark_Vertex mark_vertex,
                                Mark_Edge   mark_edge,
                                Mark_Poly   mark_poly,
                                Unsigned const kind) {
  map_mesh_allocate(nx,ny,kind) ;

  vertex_iterator iv = p2_vlist.begin() ;
  
  Real dx = (xmax - xmin) / nx ;
  Real dy = (ymax - ymin) / ny ;
  for ( Unsigned j = 0 ; j <= ny ; ++j ) {
    Real yy = ymin + j*dy ;
    for ( Unsigned i = 0 ; i <= nx ; ++i ) {
      iv -> p2_x = xmin + i*dx ;
      iv -> p2_y = yy ;
      if ( mark_vertex != NULL ) mark_vertex( *iv, Mark(i,j,nx,ny) ) ;
      ++iv ;
    }
  }

  map_mesh_internal( nx, ny, mark_vertex, mark_edge, mark_poly, kind ) ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::std_tensor_mesh(Unsigned const nx, Unsigned const ny,
                                    Mark_Vertex mark_vertex,
                                    Mark_Edge   mark_edge,
                                    Mark_Poly   mark_poly,
                                    Unsigned const kind) {
  map_mesh_allocate(nx,ny,kind) ;
  vertex_iterator iv = p2_vlist.begin() ;
  for ( Unsigned j = 0 ; j <= ny ; ++j ) {
    Real t = Real(j) / ny ;
    for ( Unsigned i = 0 ; i <= nx ; ++i ) {
      Real s = Real(i) / nx ;
      iv -> p2_x = s ;
      iv -> p2_y = t ;
      if ( mark_vertex != NULL ) mark_vertex( *iv, Mark(i,j,nx,ny) ) ;
      ++iv ;
    }
  }
  map_mesh_internal( nx, ny, mark_vertex, mark_edge, mark_poly, kind ) ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::map_mesh_allocate(Unsigned const nx,
                                      Unsigned const ny,
                                      Unsigned const kind) {
  Unsigned nxp1 = nx + 1 ;
  Unsigned nyp1 = ny + 1 ;
  Unsigned nxy  = nx * ny ;
  Unsigned nxy1 = nxp1 * nyp1 ;

  bool qt = kind == 2 && Size == 3 ;

  p2_n_vertex = nxy1 + ( qt ? nxy : 0 ) ;
  p2_n_edge   = ( nx * nyp1 + nxp1 * ny ) + (4-Size) * nxy * ( qt ? 4 : 1 ) ;
  p2_n_poly   = (5-Size) * nxy * ( qt ? 2 : 1 ) ;

  p2_vlist.resize( p2_n_vertex ) ;
  p2_elist.resize( p2_n_edge ) ;
  p2_plist.resize( p2_n_poly ) ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::map_mesh_internal(Unsigned const nx,
                                      Unsigned const ny,
                                      Mark_Vertex mark_vertex,
                                      Mark_Edge   mark_edge,
                                      Mark_Poly   mark_poly,
                                      Unsigned const kind) {

  P2MESH_MSG("enter p2_mesh::map_mesh_internal()") ;
  Unsigned i, j ;
  Unsigned nxp1 = nx + 1 ;

  // build edge (horizontal)
  edge_iterator ie = p2_elist.begin();
  for ( j = 0 ; j <= ny ; ++j ) {
    Unsigned marker = 0 ;
    if ( j == 0 ) marker = 1 ;
    else if ( j == ny ) marker = 3 ;
    for ( i = 0 ; i < nx ; ++i ) {
      Unsigned va = i + j * nxp1 ;
      Unsigned vb = va + 1 ;
      ie -> p2_v[0] = &p2_vlist[va] ;
      ie -> p2_v[1] = &p2_vlist[vb] ;
      ie -> p2_p[0] = NULL ;
      ie -> p2_p[1] = NULL ;
      if ( mark_edge != NULL ) mark_edge( *ie, marker ) ;
      ++ie ;
    }
  }

  // build edge (vertical)
  for ( i = 0 ; i <= nx ; ++i ) {
    Unsigned marker = 0 ;
    if ( i == 0 ) marker = 4 ;
    else if ( i == nx ) marker = 2 ;
    for ( j = 0 ; j < ny ; ++j ) {
      Unsigned va = i + j * nxp1 ;
      Unsigned vb = va + nxp1 ;
      ie -> p2_v[0] = &p2_vlist[va] ;
      ie -> p2_v[1] = &p2_vlist[vb] ;
      ie -> p2_p[0] = NULL ;
      ie -> p2_p[1] = NULL ;
      if ( mark_edge != NULL ) mark_edge( *ie, marker ) ;
      ++ie ;
    }
  }

  if ( Size == 4 ) {
    map_meshQ(nx, ny, mark_poly) ;
  } else {
    switch ( kind ) {
    case 0 :
      map_meshT0(nx, ny, mark_edge, mark_poly) ;
      break ;
    case 1:
      map_meshT1(nx, ny, mark_edge, mark_poly) ;
      break ;
    case 2:
      map_meshT2(nx, ny, mark_vertex, mark_edge, mark_poly) ;
      break ;
    }
  }

  JointEdges() ;
  Reorder() ;

  P2MESH_MSG("exit p2_mesh::map_mesh_internal()") ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::map_mesh(Shape_Fun      shape_fun,
                             Unsigned const nx,
                             Unsigned const ny,
                             Mark_Vertex    mark_vertex,
                             Mark_Edge      mark_edge,
                             Mark_Poly      mark_poly,
                             Unsigned const kind) {

  P2MESH_MSG("enter p2_mesh::map_mesh(...)") ;

  map_mesh_allocate(nx,ny,kind) ;

  vertex_iterator iv = p2_vlist.begin();
  for ( Unsigned j = 0 ; j <= ny ; ++j ) {
    Real t = Real(j) / ny ;
    for ( Unsigned i = 0 ; i <= nx ; ++i ) {
      Real s = Real(i) / nx ;
      shape_fun(s, t, iv -> p2_x, iv -> p2_y) ;
      if ( mark_vertex != NULL ) mark_vertex( *iv, Mark(i,j,nx,ny) ) ;
      ++iv ;
    }
  }

  map_mesh_internal(nx,ny,mark_vertex,mark_edge,mark_poly,kind) ;

  P2MESH_MSG("exit p2_mesh::map_mesh(...)") ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::read_map_mesh(char const * const file,
                                  Mark_Vertex mark_vertex,
                                  Mark_Edge   mark_edge,
                                  Mark_Poly   mark_poly,
                                  Unsigned const kind) {

  P2MESH_MSG("enter p2_mesh::red_map_mesh(...)") ;

  ifstream file_vertex( file ) ;
  if ( ! file_vertex.good() )
    msg_error( "p2_mesh::read_map_mesh(...)", "cannot open input file" ) ;

  Unsigned nx, ny ;
  file_vertex >> eatcomments >> nx >> ny >> eatline ;

  if ( nx < 1 || ny < 1 )
    msg_error( "p2_mesh::read_map_mesh(...)", "bad grid dimension" ) ;

  --nx ; --ny ;
  map_mesh_allocate(nx,ny,kind) ;

  P2MESH_MSG("enter p2_mesh::red_map_mesh(...) ...read vertices") ;

  vertex_iterator iv = p2_vlist.begin();
  for ( Unsigned j = 0 ; j <= ny ; ++j ) {
    for ( Unsigned i = 0 ; i <= nx ; ++i ) {
      file_vertex >> eatcomments >> iv -> p2_x >> iv -> p2_y >> eatline ;
      if ( mark_vertex != NULL ) mark_vertex( *iv, Mark(i,j,nx,ny) ) ;
      ++iv ;
    }
  }
  file_vertex.close() ;

  map_mesh_internal(nx,ny,mark_vertex,mark_edge,mark_poly,kind) ;

  P2MESH_MSG("enter p2_mesh::red_map_mesh(...)") ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::bbox(Real & xmin, Real & ymin,
                         Real & xmax, Real & ymax) const {
  vertex_const_iterator iv = p2_vlist.begin() ;
  xmin = xmax = iv -> p2_x ;
  ymin = ymax = iv -> p2_y ;
  for ( ++iv ; iv != p2_vlist.end() ; ++iv ) {
    xmin = min(xmin, iv -> p2_x) ; // from STL
    xmax = max(xmax, iv -> p2_x) ;
    ymin = min(ymin, iv -> p2_y) ;
    ymax = max(ymax, iv -> p2_y) ;
  }
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::map_meshQ(Unsigned const nx,
                              Unsigned const ny,
                              Mark_Poly mark_poly) {

  P2MESH_MSG("enter p2_mesh::map_meshQ(...)") ;
  Unsigned const nxp1 = nx + 1 ;
  Unsigned const nxp2 = nx + 2 ;

  P2MESH_MSG("p2_mesh::map_meshQ() ...build poly") ;
  poly_iterator ip = p2_plist.begin() ;
  for ( Unsigned j = 0 ; j < ny ; ++j ) {
    for ( Unsigned i = 0 ; i < nx ; ++i) {
      Unsigned v0 = i + j * nxp1 ;
      ip -> p2_v[0]    = &p2_vlist[v0] ;
      ip -> p2_v[1]    = &p2_vlist[v0 + 1] ;
      ip -> p2_v[2]    = &p2_vlist[v0 + nxp2] ;
      ip -> p2_v[p2_3] = &p2_vlist[v0 + nxp1] ;
      ip -> p2_e[0]    = NULL ;
      ip -> p2_e[1]    = NULL ;
      ip -> p2_e[2]    = NULL ;
      ip -> p2_e[p2_3] = NULL ;
      if ( mark_poly != NULL ) mark_poly( *ip, Mark(i,j,nx-1,ny-1) ) ;
      ++ip ;
    }
  }
  P2MESH_MSG("exit p2_mesh::map_meshQ(...)") ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::map_meshT0(Unsigned const nx,
                               Unsigned const ny,
                               Mark_Edge mark_edge,
                               Mark_Poly mark_poly) {

  P2MESH_MSG("enter p2_mesh::map_meshT0(...)") ;
  Unsigned i, j ;
  Unsigned nxp1 = nx + 1 ;

  // diagonal edge & poly
  edge_iterator ie = p2_elist.begin() + (2 * nx * ny + nx + ny) ;
  poly_iterator ip = p2_plist.begin() ;
  for ( j = 0 ; j < ny ; ++j ) {
    for ( i = 0 ; i < nx ; ++i ) {
      Unsigned v0 = i + j * nxp1 ;
      P2V * pV0 = &p2_vlist[v0] ;
      P2V * pV1 = &p2_vlist[v0 + 1] ;
      P2V * pV2 = &p2_vlist[v0 + nxp1 + 1] ;
      P2V * pV3 = &p2_vlist[v0 + nxp1] ;

      ie -> p2_v[0] = pV1 ;
      ie -> p2_v[1] = pV3 ;
      ie -> p2_p[0] = NULL ;
      ie -> p2_p[1] = NULL ;
      if ( mark_edge != NULL ) mark_edge( *ie, 0 ) ;
      ++ie ;

      ip -> p2_v[0] = pV0 ;
      ip -> p2_v[1] = pV1 ;
      ip -> p2_v[2] = pV3 ;
      ip -> p2_e[0] = NULL ;
      ip -> p2_e[1] = NULL ;
      ip -> p2_e[2] = NULL ;
      if ( mark_poly != NULL ) mark_poly( *ip, Mark(i,j,nx,ny) ) ;
      ++ip ;

      ip -> p2_v[0] = pV1 ;
      ip -> p2_v[1] = pV2 ;
      ip -> p2_v[2] = pV3 ;
      ip -> p2_e[0] = NULL ;
      ip -> p2_e[1] = NULL ;
      ip -> p2_e[2] = NULL ;
      if ( mark_poly != NULL ) mark_poly( *ip, Mark(i+1,j+1,nx,ny) ) ;
      ++ip ;
    }
  }
  P2MESH_MSG("exit p2_mesh::map_meshT0(...)") ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::map_meshT1(Unsigned const nx,
                               Unsigned const ny,
                               Mark_Edge mark_edge,
                               Mark_Poly mark_poly) {

  P2MESH_MSG("enter p2_mesh::map_meshT1(...)") ;
  Unsigned nxp1 = nx + 1 ;
  Unsigned nxp2 = nx + 2 ;

  // diagonal edge & poly
  edge_iterator ie = p2_elist.begin() + (2 * nx * ny + nx + ny) ;
  poly_iterator ip = p2_plist.begin() ;

  for ( Unsigned j = 0 ; j < ny ; ++j ) {
    for ( Unsigned i = 0 ; i < nx ; ++i ) {
      Unsigned v0 = i + j * nxp1 ;
      P2V * pV0 = &p2_vlist[v0] ;
      P2V * pV1 = &p2_vlist[v0 + 1] ;
      P2V * pV2 = &p2_vlist[v0 + nxp2] ;
      P2V * pV3 = &p2_vlist[v0 + nxp1] ;

      ie -> p2_v[0] = pV0 ;
      ie -> p2_v[1] = pV2 ;
      ie -> p2_p[0] = NULL ;
      ie -> p2_p[1] = NULL ;
      if ( mark_edge != NULL ) mark_edge( *ie, 0 ) ;
      ++ie ;

      ip -> p2_v[0] = pV0 ;
      ip -> p2_v[1] = pV1 ;
      ip -> p2_v[2] = pV2 ;
      ip -> p2_e[0] = NULL ;
      ip -> p2_e[1] = NULL ;
      ip -> p2_e[2] = NULL ;
      if ( mark_poly != NULL ) mark_poly( *ip, Mark(i+1,j,nx,ny) ) ;
      ++ip ;

      ip -> p2_v[0] = pV0 ;
      ip -> p2_v[1] = pV2 ;
      ip -> p2_v[2] = pV3 ;
      ip -> p2_e[0] = NULL ;
      ip -> p2_e[1] = NULL ;
      ip -> p2_e[2] = NULL ;
      if ( mark_poly != NULL ) mark_poly( *ip, Mark(i,j+1,nx,ny) ) ;
      ++ip ;
    }
  }
  P2MESH_MSG("exit p2_mesh::map_meshT1(...)") ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::map_meshT2(Unsigned const  nx,
                               Unsigned const  ny,
                               Mark_Vertex mark_vertex,
                               Mark_Edge   mark_edge,
                               Mark_Poly   mark_poly) {

  P2MESH_MSG("enter p2_mesh::map_meshT2(...)") ;
  Unsigned i, j ;
  Unsigned nxp1 = nx + 1 ;
  Unsigned nyp1 = ny + 1 ;
  Unsigned nxy1 = nxp1 * nyp1 ;

  // add extra vertex diagonals & poly
  vertex_iterator iv = p2_vlist.begin() + nxy1 ;
  edge_iterator   ie = p2_elist.begin() + (2 * nx * ny + nx + ny) ;
  poly_iterator   ip = p2_plist.begin() ;

  for ( j = 0 ; j < ny ; ++j ) {
    for ( i = 0 ; i < nx ; ++i ) {
      Unsigned va = i + j * nxp1 ;
      P2V * pVa = & p2_vlist[ va ] ;
      P2V * pVb = & p2_vlist[ va + 1 ] ;
      P2V * pVc = & p2_vlist[ va + nxp1 + 1 ] ;
      P2V * pVd = & p2_vlist[ va + nxp1 ] ;
      P2V * pVe = & p2_vlist[ nxy1 + i + j * nx ] ;

      iv -> p2_x = 0.25 * (pVa -> x() + pVb -> x() + pVc -> x() + pVd -> x());
      iv -> p2_y = 0.25 * (pVa -> y() + pVb -> y() + pVc -> y() + pVd -> y());
      if ( mark_vertex != NULL ) mark_vertex( *iv, 0 ) ;
      ++iv ;

      ie -> p2_v[0] = pVa ;
      ie -> p2_v[1] = pVe ;
      ie -> p2_p[0] = NULL ;
      ie -> p2_p[1] = NULL ;
      if ( mark_edge != NULL ) mark_edge( *ie, 0 ) ;
      ++ie ;

      ie -> p2_v[0] = pVb ;
      ie -> p2_v[1] = pVe ;
      ie -> p2_p[0] = NULL ;
      ie -> p2_p[1] = NULL ;
      if ( mark_edge != NULL ) mark_edge( *ie, 0 ) ;
      ++ie ;

      ie -> p2_v[0] = pVc ;
      ie -> p2_v[1] = pVe ;
      ie -> p2_p[0] = NULL ;
      ie -> p2_p[1] = NULL ;
      if ( mark_edge != NULL ) mark_edge( *ie, 0 ) ;
      ++ie ;

      ie -> p2_v[0] = pVd ;
      ie -> p2_v[1] = pVe ;
      ie -> p2_p[0] = NULL ;
      ie -> p2_p[1] = NULL ;
      if ( mark_edge != NULL ) mark_edge( *ie, 0 ) ;
      ++ie ;

      ip -> p2_v[0] = pVa ;
      ip -> p2_v[1] = pVb ;
      ip -> p2_v[2] = pVe ;
      ip -> p2_e[0] = NULL ;
      ip -> p2_e[1] = NULL ;
      ip -> p2_e[2] = NULL ;
      if ( mark_poly != NULL ) mark_poly( *ip, ( j==0 ? 1 : 0 ) ) ;
      ++ip ;

      ip -> p2_v[0] = pVb ;
      ip -> p2_v[1] = pVc ;
      ip -> p2_v[2] = pVe ;
      ip -> p2_e[0] = NULL ;
      ip -> p2_e[1] = NULL ;
      ip -> p2_e[2] = NULL ;
      if ( mark_poly != NULL ) mark_poly( *ip, ( i==nx-1 ? 2 : 0 ) ) ;
      ++ip ;

      ip -> p2_v[0] = pVc ;
      ip -> p2_v[1] = pVd ;
      ip -> p2_v[2] = pVe ;
      ip -> p2_e[0] = NULL ;
      ip -> p2_e[1] = NULL ;
      ip -> p2_e[2] = NULL ;
      if ( mark_poly != NULL ) mark_poly( *ip, ( j==ny-1 ? 3 : 0 ) ) ;
      ++ip ;

      ip -> p2_v[0] = pVd ;
      ip -> p2_v[1] = pVa ;
      ip -> p2_v[2] = pVe ;
      ip -> p2_e[0] = NULL ;
      ip -> p2_e[1] = NULL ;
      ip -> p2_e[2] = NULL ;
      if ( mark_poly != NULL ) mark_poly( *ip, ( i==0 ? 4 : 0 ) ) ;
      ++ip ;

    }
  }
  P2MESH_MSG("exit p2_mesh::map_meshT2(...)") ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::build_mesh(Unsigned const nv,
                               Real     const *XY,
                               Vmark    const *mv,
                               Mark_Vertex mark_vertex,

                               Unsigned const ne,
                               Unsigned const *E,
                               Emark    const *me,
                               Mark_Edge  mark_edge,

                               Unsigned const np,
                               Unsigned const *P,
                               Pmark    const *mp,
                               Mark_Poly  mark_poly,

                               Unsigned const base) {

  P2MESH_MSG("enter p2_mesh::build_mesh(...)") ;

  p2_vlist.resize( p2_n_vertex = nv ) ;

  for ( vertex_iterator iv = p2_vlist.begin() ;
        iv != p2_vlist.end() ;
        ++iv ) {
    iv -> p2_x = *XY++ ;
    iv -> p2_y = *XY++ ;
    if ( mark_vertex != NULL && mv != NULL )
      { mark_vertex( *iv, *mv ) ; ++mv ; }
  }
  
  p2_elist.resize( p2_n_edge = ne ) ;

  for ( edge_iterator ie = p2_elist.begin() ;
        ie != p2_elist.end() ;
        ++ie ) {
    Unsigned v0 = (*E++ - base) ;
    Unsigned v1 = (*E++ - base) ;

    if ( v0 >= p2_n_vertex || v1 >= p2_n_vertex )
      msg_error("p2_mesh::build_mesh(...)",
                "bad edge definition in edge list") ;

    ie -> p2_v[0] = &p2_vlist[ v0 ] ;
    ie -> p2_v[1] = &p2_vlist[ v1 ] ;
    ie -> p2_p[0] = NULL ;
    ie -> p2_p[1] = NULL ;
    if ( mark_edge != NULL && me != NULL )
      { mark_edge( *ie, *me ) ; ++me ; }
  }

  p2_plist.resize( p2_n_poly = np ) ;

  for ( poly_iterator ip = p2_plist.begin() ;
        ip != p2_plist.end() ;
        ++ip ) {
    for ( Unsigned v = 0 ; v < Size ; ++v ) {
      Unsigned nv_loc = (*P++ - base) ;
      if ( nv_loc >= p2_n_vertex )
        msg_error("p2_mesh::build_mesh(...)",
                  "bad polygon definition in polygon list") ;
      ip -> p2_v[ v ] = &p2_vlist[ nv_loc ] ;
      ip -> p2_e[ v ] = NULL ;
    }
    if ( mark_poly != NULL && mp != NULL )
      { mark_poly( *ip, *mp ) ; ++mp ; }
  }
  
  JointEdges() ;
  Reorder() ;

  P2MESH_MSG("exit p2_mesh::build_mesh(...)") ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::read_mesh(char const * const file_name,
                              Mark_Vertex mark_vertex,
                              Mark_Edge   mark_edge,
                              Mark_Poly   mark_poly,
                              Unsigned const base) {

  P2MESH_MSG("enter p2_mesh::read_mesh(...)") ;

  string file_v = string(file_name) + ".node" ;
  string file_e = string(file_name) + ".edge" ;
  string file_p = string(file_name) + ".ele" ;

  ifstream file_vertex( file_v.c_str() ) ;
  ifstream file_edge  ( file_e.c_str() ) ;
  ifstream file_poly  ( file_p.c_str() ) ;

  if ( ! file_vertex.good() )
    msg_error( "p2_mesh::read_mesh(...)", "cannot open nodes file" ) ;

  if ( ! file_poly.good() )
    msg_error( "p2_mesh::read_mesh(...)", "cannot open polygons file" ) ;

  file_vertex >> eatcomments >> p2_n_vertex >> eatline ;
  p2_vlist.resize(p2_n_vertex) ;

  for ( vertex_iterator iv = p2_vlist.begin() ;
        iv != p2_vlist.end() ;
        ++iv ) {
    Unsigned nv ;
    file_vertex >> eatcomments >> nv >> iv -> p2_x >> iv -> p2_y ;
    nv -= base ;
    if ( nv >= p2_n_vertex || ! file_vertex.good() )
      msg_error( "p2_mesh::read_mesh(...)",
                 "error in reading vertex coordinates" ) ;
    if ( mark_vertex != NULL )
      { Vmark info ; file_vertex >> info ; mark_vertex( *iv, info ) ; }
    file_vertex >> eatline ;
  }
  file_vertex.close() ;

  // read poly
  file_poly >> eatcomments >> p2_n_poly >> eatline ;
  p2_plist.resize(p2_n_poly) ;

  for ( poly_iterator ip = p2_plist.begin() ;
        ip != p2_plist.end() ;
        ++ip ) {
    Unsigned np, nv ;
    file_poly >> eatcomments >> np ;
    np -= base ;
    if ( np >= p2_n_poly || ! file_poly.good() )
    msg_error( "p2_mesh::read_mesh(...)",
               "error in reading polygon definitions" ) ;

    for ( Unsigned v = 0 ; v < Size ; ++v ) {
      file_poly >> nv ; nv -= base ;
      if ( nv >= p2_n_vertex || ! file_poly.good() )
        msg_error( "p2_mesh::read_mesh(...)",
                   "error in reading vertex numbers for the polygon" ) ;
      ip -> p2_v[ v ] = &p2_vlist[nv] ;
      ip -> p2_e[ v ] = NULL ;
    }
    if ( mark_poly != NULL )
      { Pmark info ; file_poly >> info ; mark_poly( *ip, info ) ; }
    file_poly >> eatline ;
  }

  file_poly. close() ;

  if ( file_edge.good() ) {
    file_edge >> eatcomments >> p2_n_edge >> eatline ;
    p2_elist.resize(p2_n_edge) ;
    for ( edge_iterator ie = p2_elist.begin() ;
      ie != p2_elist.end() ;
      ++ie ) {
      Unsigned ne, v0, v1 ;
      file_edge >> eatcomments >> ne >> v0 >> v1 ;
      ne -= base ;
      v0 -= base ;
      v1 -= base ;
      if ( ne >= p2_n_edge || ! file_edge.good() )
         msg_error( "p2_mesh::read_mesh(...)",
                    "error in reading edge definition" ) ;
      if ( v0 >= p2_n_vertex || v1 >= p2_n_vertex )
         msg_error( "p2_mesh::read_mesh(...)",
                    "error in reading vertex numbers for the edge" ) ;
      ie -> p2_v[0] = &p2_vlist[v0] ;
      ie -> p2_v[1] = &p2_vlist[v1] ;
      ie -> p2_p[0] = NULL ;
      ie -> p2_p[1] = NULL ;
      if ( mark_edge != NULL )
        { Emark info ; file_edge >> info ; mark_edge( *ie, info ) ; }
      file_edge >> eatline ;
    }
    file_edge.close() ;
    JointEdges() ;
  } else {
    BuildEdges() ;
  }

  Reorder() ;
  P2MESH_MSG("exit p2_mesh::read_mesh(...)") ;
}

template <typename P2_COMMON>
bool
p2_mesh<P2_COMMON>::test_mesh() const {

  char const * msg_err = "" ;
  Unsigned              i, j, np ;
  vertex_const_iterator iv ;
  edge_const_iterator   ie ;
  poly_const_iterator   ip ;

  for ( ip = poly_begin() ; ip != poly_end() ; ++ip ) {
    for ( Unsigned nve = 0 ; nve < Size ; ++nve ) {
      if ( ip -> p2_e[nve] == NULL || ip -> p2_v[nve] == NULL ) {
        msg_err = "incomplete polygon found" ;
        goto error_found ;
      }
    } 
  }

  for ( ie = edge_begin() ; ie != edge_end() ; ++ie ) {
  
    P2V * pa = ie -> p2_v[0] ;
    P2V * pb = ie -> p2_v[1] ;
    
    if ( pa == NULL && pb == NULL ) {
      msg_err = "incomplete edge found" ;
      goto error_found ;
    }

    for ( np = 0 ; np < Unsigned(2) ; ++np ) {
      if ( ie -> p2_p[np] != NULL ) {
        P2P * pp = ie -> p2_p[np] ;
        bool ok_edge = false ;
        bool ok_va   = false ;
        bool ok_vb   = false ;
        for ( i = 0 ; i < pp -> n_vertex() ; ++i ) {
          if ( pp -> p2_e[i] == &*ie ) ok_edge = true ;
          if ( pp -> p2_v[i] == pa   ) ok_va   = true ;
          if ( pp -> p2_v[i] == pb   ) ok_vb   = true ;
        }
        if ( !(ok_edge && ok_va && ok_vb) ) {
          msg_err = "edge with bad connection found" ;
          goto error_found ;
        } 
      }
    }
  }

  for ( ip = poly_begin() ; ip != poly_end() ; ++ip ) {
    if ( ip -> area() <= 0 ) {
      msg_err = "polygon with negative area found" ;
      goto error_found ;
    }
  }

  if ( List ) {
    for ( iv = vertex_begin() ; iv != vertex_end() ; ++iv ) {
  
      if ( iv -> n_vertex() != iv -> n_edge() ) {
        msg_err = "n_vertex != n_edge on vertex list" ;
        goto error_found ;
      }
      
      for ( i = 0 ; i < iv -> n_vertex() ; ++i ) {
      
        if ( iv -> psV(i) == NULL ) {
          msg_err = "found a NULL pointer in the vertex --> vertex list" ;
          goto error_found ;
        }
      
        if ( iv -> psE(i) == NULL ) {
          msg_err = "found a NULL pointer in the vertex --> edge list" ;
          goto error_found ;
        }

        P2E const * pE = static_cast<P2E*>(iv -> psE(i)) ;
        P2V const * pV = static_cast<P2V*>(iv -> psV(i)) ;

        bool ok1 = pE -> p2_v[0] == &*iv && pE -> p2_v[1] == pV ;
        bool ok2 = pE -> p2_v[0] == pV   && pE -> p2_v[1] == &*iv ;

        if ( ! ( ok1 || ok2 )  ) {
          msg_err = "bad vertex list" ;
          goto error_found ;
        }
      }

      for ( i = 0 ; i < iv -> n_poly() ; ++i ) {
        if ( iv -> psP(i) == NULL ) {
          msg_err = "found a NULL pointer in the vertex --> poly list" ;
          goto error_found ;
        }

        P2P * pP = static_cast<P2P*>(iv -> psP(i)) ;
        bool ok = false ;
        for ( j = 0 ; j < Size ; ++j )
          ok = ok || (pP -> p2_v[j] == &*iv) ;
        if ( !ok ) {
          msg_err = "bad vertex --> polygon list" ;
          goto error_found ;
        }
      }
    }
  }
  
  return true ;

error_found:
  cerr << endl << "p2_mesh::test_mesh() -- " << msg_err << endl ;
  return false ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::report(ostream & s) const {
  s << endl
    << "            p2_mesh statistics" << endl
    << "            Polygon Type = "
    << ( Size == 3 ? "Triangle" : "Quadrilateral" ) << endl
    << "           +----------+----------+----------+" << endl
    << "           | Total    | Internal | Boundary |" << endl
    << "+----------+----------+----------+----------+" << endl
    << "| Vertices | " << setw(8) << p2_n_vertex  << " | "
                       << setw(8) << p2_n_ivertex << " | "
                       << setw(8) << p2_n_bvertex << " | " << endl
    << "+----------+----------+----------+----------+" << endl
    << "| Edges    | " << setw(8) << p2_n_edge    << " | "
                       << setw(8) << p2_n_iedge   << " | "
                       << setw(8) << p2_n_bedge   << " | " << endl
    << "+----------+----------+----------+----------+" << endl
    << "| Polygons | " << setw(8) << p2_n_poly    << " | "
                       << setw(8) << p2_n_ipoly   << " | "
                       << setw(8) << p2_n_bpoly   << " | " << endl
    << "+----------+----------+----------+----------+" << endl ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::print(ostream & s, Unsigned const base) const {
  vertex_const_iterator iv ;
  edge_const_iterator   ie ;
  poly_const_iterator   ip ;

  s << endl
    << "P2MESH STRUCTURE:" << endl << endl
    << "VERTICES" << endl ;

  for ( iv = vertex_begin() ; iv != vertex_end() ; ++iv )
    print(s, *iv, base ) ;

  s << endl << "EDGES" << endl ;
  for ( ie = edge_begin() ; ie != edge_end() ; ++ie )
     print(s, *ie, base ) ;

  s << endl << "POLYS" << endl ;
  for ( ip = poly_begin() ; ip != poly_end() ; ++ip )
    print(s, *ip, base ) ;

  s << endl
    << "polygons = " << p2_n_poly   << endl
    << "edges    = " << p2_n_edge   << endl
    << "vertices = " << p2_n_vertex << endl ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::print(ostream & s, P2V const & V, Unsigned const base)
const {
  Unsigned nv = base + local_number(V) ;
  s << "vertex:" << setw(5) << nv <<
       " ( " << setw(5) << V.x() << " , " <<
                setw(5) << V.y() << " )" << endl ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::print(ostream & s, P2E const & E, Unsigned const base)
const {
  Unsigned ne = base + local_number( E ) ;
  Unsigned na = base + local_number( *E.p2_v[0] ) ;
  Unsigned nb = base + local_number( *E.p2_v[1] ) ;

  s << "edge:" << setw(5) << ne
    << "    V[ " << setw(5) << na << " " << setw(5) << nb << " ]"
    << "    P{ " ;

  if ( E.p2_p[0] == NULL ) s << "  *  " ;
  else                     s << setw(5) << base + local_number( *E.p2_p[0] ) ;

  s << " " ;

  if ( E.p2_p[1] == NULL ) s << "  *  " ;
  else                     s << setw(5) << base + local_number( *E.p2_p[1] ) ;

  s << " }" << endl ;
}

template <typename P2_COMMON>
void
p2_mesh<P2_COMMON>::print(ostream & s, P2P const & P, Unsigned const base)
const {
  s << "poly:" << setw(5) << base + local_number( P ) << "    V[ " ;
  for ( Unsigned nv = 0 ; nv < Size ; ++nv )
    s << setw(5) << base + local_number( * P.p2_v[nv] ) << " " ;
  s << "]    E[ " ;
  for ( Unsigned ne = 0 ; ne < Size ; ++ne )
    s << setw(5) << base + local_number( * P.p2_e[ne] ) << " " ;
  s << "]" << endl ;
}

# endif

// end of file p2mesh.hh
// by Enrico Bertolazzi & Gianmarco Manzini
