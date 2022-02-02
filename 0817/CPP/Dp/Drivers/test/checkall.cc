/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  P2MESH : DRIVER FOR CHECKING THE PACKAGE COMPILATION                    |
 |                                                                          |
 |  date         : 1999, 11 October                                         |
 |  version      : 1.0                                                      |
 |  file         : checkall.cc                                              |
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
 |   instantiate all the methods of P2MESH and check the template           |
 |   capability of the compiler.                                            |
 |                                                                          |
\*--------------------------------------------------------------------------*/

# ifndef SIZE
  # define SIZE 4
# endif

# ifndef LIST
  # define LIST true
# endif

# include "p2mesh.hh"

typedef double   Real ;
typedef unsigned Unsigned ;
typedef int      Integer ;
typedef double   Real ;

typedef long     Vmark ;
typedef long     Emark ;
typedef long     Pmark ;

class Vertex ;
class Edge   ;
class Poly   ;
class Mesh   ;

class Common : public p2_common<Vertex,Edge,Poly,Mesh,
                                SIZE,LIST,Real,Integer,Unsigned,
                                Vmark, Emark, Pmark> { } ; 

class Vertex   : public p2_vertex<Common> {} ;
class Edge     : public p2_edge<Common>   {} ;
class Poly     : public p2_poly<Common>   {} ;
class Mesh     : public p2_mesh<Common>   {} ;

static
void
shape(Real const &, Real const &, Real &, Real &) {}

static
void
mark_vertex(Vertex &, Vmark const &) {}

static
void
mark_edge(Edge &, Emark const &) {}

static
void
mark_poly(Poly &, Pmark const &) {}

int
main() {

  Mesh mesh ;
  
  mesh . read_mesh("pippo", NULL, NULL, NULL, 1) ;

  Real xmin, ymin, xmax, ymax ;
  mesh . bbox(xmin, ymin, xmax, ymax) ;

  cout << mesh . n_vertex() << endl ;
  cout << mesh . n_bvertex() << endl ;
  cout << mesh . n_ivertex() << endl ;

  cout << mesh . n_edge () << endl ;
  cout << mesh . n_bedge() << endl ;
  cout << mesh . n_iedge() << endl ;

  cout << mesh . n_poly () << endl ;
  cout << mesh . n_bpoly() << endl ;
  cout << mesh . n_ipoly() << endl ;
  
  Vertex       & v  = mesh . vertex(0) ;
  Vertex const & cv = mesh . vertex(0) ;

  Edge         & e  = mesh . edge(0) ;
  Edge const   & ce = mesh . edge(0) ;
  
  Poly         & p  = mesh . poly(0) ;
  Poly const   & cp = mesh . poly(0) ;

  Unsigned n1 = mesh . local_number(v) ;
  Unsigned n2 = mesh . local_number(cv) ;
  Unsigned n3 = mesh . local_number(e) ;
  Unsigned n4 = mesh . local_number(ce) ;
  Unsigned n5 = mesh . local_number(p) ;
  Unsigned n6 = mesh . local_number(cp) ;

  // iterators for vertex
  Mesh::vertex_iterator       vb  = mesh . vertex_begin() ;
  Mesh::vertex_const_iterator cvb = mesh . vertex_begin() ;
  Mesh::vertex_iterator       ve  = mesh . vertex_end() ;
  Mesh::vertex_const_iterator cve = mesh . vertex_end() ;

  Mesh::vertex_iterator       vb1  = mesh . bvertex_begin() ;
  Mesh::vertex_const_iterator cvb1 = mesh . bvertex_begin() ;
  Mesh::vertex_iterator       ve1  = mesh . bvertex_end() ;
  Mesh::vertex_const_iterator cve1 = mesh . bvertex_end() ;

  Mesh::vertex_iterator       vb2  = mesh . ivertex_begin() ;
  Mesh::vertex_const_iterator cvb2 = mesh . ivertex_begin() ;
  Mesh::vertex_iterator       ve2  = mesh . ivertex_end() ;
  Mesh::vertex_const_iterator cve2 = mesh . ivertex_end() ;

  // iterators for edge
  Mesh::edge_iterator       eb  = mesh . edge_begin() ;
  Mesh::edge_const_iterator ceb = mesh . edge_begin() ;
  Mesh::edge_iterator       ee  = mesh . edge_end() ;
  Mesh::edge_const_iterator cee = mesh . edge_end() ;

  Mesh::edge_iterator       eb1  = mesh . bedge_begin() ;
  Mesh::edge_const_iterator ceb1 = mesh . bedge_begin() ;
  Mesh::edge_iterator       ee1  = mesh . bedge_end() ;
  Mesh::edge_const_iterator cee1 = mesh . bedge_end() ;

  Mesh::edge_iterator       eb2  = mesh . iedge_begin() ;
  Mesh::edge_const_iterator ceb2 = mesh . iedge_begin() ;
  Mesh::edge_iterator       ee2  = mesh . iedge_end() ;
  Mesh::edge_const_iterator cee2 = mesh . iedge_end() ;

  // iterators for poly
  Mesh::poly_iterator       pb  = mesh . poly_begin() ;
  Mesh::poly_const_iterator cpb = mesh . poly_begin() ;
  Mesh::poly_iterator       pe  = mesh . poly_end() ;
  Mesh::poly_const_iterator cpe = mesh . poly_end() ;

  Mesh::poly_iterator       pb1  = mesh . bpoly_begin() ;
  Mesh::poly_const_iterator cpb1 = mesh . bpoly_begin() ;
  Mesh::poly_iterator       pe1  = mesh . bpoly_end() ;
  Mesh::poly_const_iterator cpe1 = mesh . bpoly_end() ;

  Mesh::poly_iterator       pb2  = mesh . ipoly_begin() ;
  Mesh::poly_const_iterator cpb2 = mesh . ipoly_begin() ;
  Mesh::poly_iterator       pe2  = mesh . ipoly_end() ;
  Mesh::poly_const_iterator cpe2 = mesh . ipoly_end() ;

  mesh . map_mesh(shape, 10, 20, mark_vertex, mark_edge, mark_poly, 0) ;
  mesh . map_mesh(shape, 10, 20, mark_vertex, mark_edge, mark_poly, 1) ;
  mesh . map_mesh(shape, 10, 20, mark_vertex, mark_edge, mark_poly, 2) ;
  mesh . map_mesh(shape, 10, 20, mark_vertex, mark_edge, mark_poly, 3) ;

  mesh . tensor_mesh(0,1,0,1,10,10, mark_vertex, mark_edge, mark_poly, 0) ;
  mesh . tensor_mesh(0,1,0,1,10,10, mark_vertex, mark_edge, mark_poly, 1) ;
  mesh . tensor_mesh(0,1,0,1,10,10, mark_vertex, mark_edge, mark_poly, 2) ;
  mesh . tensor_mesh(0,1,0,1,10,10, mark_vertex, mark_edge, mark_poly, 3) ;

  mesh . std_tensor_mesh(10, 10, mark_vertex, mark_edge, mark_poly, 0) ;
  mesh . std_tensor_mesh(10, 10, mark_vertex, mark_edge, mark_poly, 1) ;
  mesh . std_tensor_mesh(10, 10, mark_vertex, mark_edge, mark_poly, 2) ;
  mesh . std_tensor_mesh(10, 10, mark_vertex, mark_edge, mark_poly, 3) ;

  Real     * XY = NULL ;
  Unsigned * E  = NULL ;
  Unsigned * P  = NULL ;
  Vmark    * mv = NULL ;
  Emark    * me = NULL ;
  Pmark    * mp = NULL ;

  mesh . build_mesh(10, XY, mv, mark_vertex,
                    10, E, me, mark_edge,
                    10, P, mp, mark_poly,
                    0 ) ;

  mesh . read_mesh("pippo", mark_vertex, mark_edge, mark_poly, 0) ;
  mesh . read_map_mesh("pippo", mark_vertex, mark_edge, mark_poly,0) ;

  mesh . report(cout) ;
  mesh . test_mesh() ;
  
  mesh . print(cout, 0) ;
  mesh . print(cout, mesh . edge(0), 0) ;
  mesh . print(cout, mesh . vertex(0), 0) ;
  mesh . print(cout, mesh . poly(0), 0) ;

  // vertex

  cout << v . x()  << " " << v . y()  << endl ;
  cout << cv . x() << " " << cv . y() << endl ;
  cout << v . n_vertex() << endl ;
  cout << v . n_edge() << endl ;
  cout << v . n_poly() << endl ;
  cout << cv . n_vertex() << endl ;
  cout << cv . n_edge() << endl ;
  cout << cv . n_poly() << endl ;

  { Vertex       & vv  = v . vertex(0) ;
    Vertex const & cvv = v . vertex(0) ; }
  { Vertex const & cvv = cv . vertex(0) ; }

  { Edge       & ee  = v . edge(0) ;
    Edge const & cee = v . edge(0) ; }
  { Edge const & cee = cv . edge(0) ; }

  { Poly       & pp  = v . poly(0) ;
    Poly const & cpp = v . poly(0) ; }
  { Poly const & pp  = cv . poly(0) ; }

  { Unsigned lc = v . local_number(v) ; }
  { Unsigned lc = v . local_number(e) ; }
  { Unsigned lc = v . local_number(p) ; }

  { Unsigned lc = v . local_number(cv) ; }
  { Unsigned lc = v . local_number(ce) ; }
  { Unsigned lc = v . local_number(cp) ; }

  { Unsigned lc = cv . local_number(v) ; }
  { Unsigned lc = cv . local_number(e) ; }
  { Unsigned lc = cv . local_number(p) ; }

  { Unsigned lc = cv . local_number(cv) ; }
  { Unsigned lc = cv . local_number(ce) ; }
  { Unsigned lc = cv . local_number(cp) ; }
  
  // edge

  cout << e . n_vertex() << endl ;
  cout << e . n_edge() << endl ;
  cout << e . n_poly() << endl ;
  cout << ce . n_vertex() << endl ;
  cout << ce . n_edge() << endl ;
  cout << ce . n_poly() << endl ;

  { Vertex       & vv  = e . vertex(0) ;
    Vertex const & cvv = e . vertex(0) ; }
  { Vertex const & cvv = ce . vertex(0) ; }

  { Edge       & ee  = e . edge(0) ;
    Edge const & cee = e . edge(0) ; }
  { Edge const & cee = ce . edge(0) ; }

  { Poly       & pp  = e . poly(0) ;
    Poly const & cpp = e . poly(0) ; }
  { Poly const & pp  = ce . poly(0) ; }

  { bool ok = e  . ok_poly(1) ; }
  { bool ok = ce . ok_poly(1) ; }
  
  { Unsigned lc = e . local_number(v) ; }
  { Unsigned lc = e . local_number(e) ; }
  { Unsigned lc = e . local_number(p) ; }

  { Unsigned lc = e . local_number(cv) ; }
  { Unsigned lc = e . local_number(ce) ; }
  { Unsigned lc = e . local_number(cp) ; }

  { Unsigned lc = ce . local_number(v) ; }
  { Unsigned lc = ce . local_number(e) ; }
  { Unsigned lc = ce . local_number(p) ; }

  { Unsigned lc = ce . local_number(cv) ; }
  { Unsigned lc = ce . local_number(ce) ; }
  { Unsigned lc = ce . local_number(cp) ; }

  cout << e . x(0)  << " " << e . y(1) << endl ;
  cout << e . xm() << " " << e . ym() << endl ;
  cout << e . xt(0.1) << " " << e . yt(0.2) << endl ;
  cout << e . nx() << " " << e . ny() << endl ;
  cout << e . tx() << " " << e . ty() << endl ;
  cout << e . length() << endl ;

  cout << ce . x(0)  << " " << ce . y(1) << endl ;
  cout << ce . xm() << " " << ce . ym() << endl ;
  cout << ce . xt(0.1) << " " << ce . yt(0.2) << endl ;
  cout << ce . nx() << " " << ce . ny() << endl ;
  cout << ce . tx() << " " << ce . ty() << endl ;
  cout << ce . length() << endl ;
  
  // poly

  cout << p . n_vertex() << endl ;
  cout << p . n_edge() << endl ;
  cout << p . n_poly() << endl ;
  cout << cp . n_vertex() << endl ;
  cout << cp . n_edge() << endl ;
  cout << cp . n_poly() << endl ;

  { Vertex       & vv  = p . vertex(0) ;
    Vertex const & cvv = p . vertex(0) ; }
  { Vertex const & cvv = cp . vertex(0) ; }

  { Edge       & ee  = p . edge(0) ;
    Edge const & cee = p . edge(0) ; }
  { Edge const & cee = cp . edge(0) ; }

  { Poly       & pp  = p . poly(0) ;
    Poly const & cpp = p . poly(0) ; }
  { Poly const & pp  = cp . poly(0) ; }

  { Unsigned lc = p . local_number(v) ; }
  { Unsigned lc = p . local_number(e) ; }
  { Unsigned lc = p . local_number(p) ; }

  { Unsigned lc = p . local_number(cv) ; }
  { Unsigned lc = p . local_number(ce) ; }
  { Unsigned lc = p . local_number(cp) ; }

  { Unsigned lc = cp . local_number(v) ; }
  { Unsigned lc = cp . local_number(e) ; }
  { Unsigned lc = cp . local_number(p) ; }

  { Unsigned lc = cp . local_number(cv) ; }
  { Unsigned lc = cp . local_number(ce) ; }
  { Unsigned lc = cp . local_number(cp) ; }


  { bool ok = p  . ok_poly(0) ; }
  { bool ok = cp . ok_poly(0) ; }

  { bool ok = p  . ok_oriented(0) ; }
  { bool ok = cp . ok_oriented(0) ; }

  cout << p . x(0) << " " << p . y(0) << endl ;
  cout << p . xc() << " " << p . yc() << endl ;
  cout << p . area() << endl ;
  cout << p . nx(0) << " " << p . ny(1) << endl ;
  cout << p . tx(0) << " " << p . ty(0) << endl ;
  cout << p . xm(1) << " " << p . ym(2) << endl ;
  cout << p . xt(1,0.3) << " " << p . yt(2,0.3) << endl ;

  Real s, t, xx, yy, J[2][2], gxy[2], gst[2] ;
  gxy[0] = gxy[1] = gst[0] = gst[1] = 0 ;  
  s = t = 0 ;
  p . jacobian(s, t, J)         ; cp . jacobian(s, t, J) ;
  p . inverse_jacobian(s, t, J) ; cp . inverse_jacobian(s, t, J) ;
  p . st_to_xy(s, t, xx, yy)    ; cp . st_to_xy(s, t, xx, yy) ;
  p . xy_to_st(xx, yy, s, t)    ; cp . xy_to_st(xx, yy, s, t) ;

  p  . grad_st_to_xy(s, t, gst, gxy) ;
  cp . grad_st_to_xy(s, t, gst, gxy) ;

  p  . grad_xy_to_st(s, t, gxy, gst) ;
  cp . grad_xy_to_st(s, t, gxy, gst) ;
  
  // iterator
  {
    Iterator<Vertex> iv, iv1(mesh), iv2(mesh,2) ;
    Iterator<Edge>   ie, ie1(mesh), ie2(mesh,2) ;
    Iterator<Poly>   ip, ip1(mesh), ip2(mesh,2) ;
  
    iv . set_loop(mesh, 1) ;
    ie . set_loop(mesh, 1) ;
    ip . set_loop(mesh, 1) ;
  
    { iv . begin() ; bool ok = iv . end_of_loop() ; }
    { ie . begin() ; bool ok = ie . end_of_loop() ; }
    { ip . begin() ; bool ok = ip . end_of_loop() ; }
  
    ++iv ; ++ ie ; ++ip ;
  
    iv -> x() ; ie -> x(1) ; ip -> x(1) ;
    (*iv) . x() ; (*ie) . x(1) ; (*ip) . x(1) ;
  }
  
  {
    CIterator<Vertex> iv, iv1(mesh), iv2(mesh,2) ;
    CIterator<Edge>   ie, ie1(mesh), ie2(mesh,2) ;
    CIterator<Poly>   ip, ip1(mesh), ip2(mesh,2) ;
  
    iv . set_loop(mesh, 1) ;
    ie . set_loop(mesh, 1) ;
    ip . set_loop(mesh, 1) ;
  
    { iv . begin() ; bool ok = iv . end_of_loop() ; }
    { ie . begin() ; bool ok = ie . end_of_loop() ; }
    { ip . begin() ; bool ok = ip . end_of_loop() ; }
  
    ++iv ; ++ ie ; ++ip ;
  
    iv -> x() ; ie -> x(1) ; ip -> x(1) ;
    (*iv) . x() ; (*ie) . x(1) ; (*ip) . x(1) ;
  } 
}

// end of file checkall.cc
// by Enrico Bertolazzi & Gianmarco Manzini
