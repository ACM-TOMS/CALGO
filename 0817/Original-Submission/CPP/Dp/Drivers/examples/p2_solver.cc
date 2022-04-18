/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  P2MESH : a polygon-based mesh manager                                   |
 |           example driver program                                         |
 |                                                                          |
 |  date         : 1999, 11 October                                         |
 |  version      : 1.0                                                      |
 |  file         : p2.cc                                                    |
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
 |    numerical solution of the Laplace problem on a square domain using    |
 |    P2 conforming base polynomials.                                       |
 |                                                                          |
 |  references:                                                             |
 |    primer.ps  -- beginner introduction with commented examples and       |
 |                  references therein                                      |
 |                                                                          |
\*--------------------------------------------------------------------------*/

// include p2mesh template library
# include "p2mesh.hh"

// declare the name of the user defined classes
class Vertex ;
class Edge ;
class Triangle ;
class Mesh ;
class Elliptic_Solver ;

typedef double (*pFun)(double const & x, double const & y) ;

// setup common code & data
class Common : public p2_common<Vertex,Edge,Triangle,Mesh> {
protected:
  static unsigned const degree_of_freedom = 6 ;
  
  static void shape(unsigned const,
                    double const &,
                    double const &,
                    double &) ;
                    
  static void shape_grad(unsigned const,
                         double const &,
                         double const &,
                         double [2]) ;
} ;

// define Vertex class
class Vertex : public p2_vertex<Common> {
public:
  unsigned EqNumber    (Mesh const &) const ;
  bool     IsOnBoundary(Mesh const &) const ;
} ;

// define Edge class
class Edge : public p2_edge<Common> {
public:
  unsigned EqNumber    (Mesh const &) const ;
  bool     IsOnBoundary(Mesh const &) const ;
} ;

// define Triangle class
class Triangle : public p2_poly<Common> {
public:
  void   eval_JJT(double[2][2], double &) const ;
  double eval_int_f(unsigned const, pFun, double const &) const ;
  double eval_int_grad(unsigned const, unsigned const,
                       double const [2][2], double const &) const ;
                       
  unsigned EqNumber    (Mesh const &, Unsigned const) const ;
  bool     IsOnBoundary(Mesh const &, Unsigned const) const ;
} ;

// define Mesh class
class Mesh : public p2_mesh<Common> {} ;

// define the solver class
class Elliptic_Solver : public Common {

private: 
  Mesh   mesh ;
  double **mat ;
  double *sol, *rhs ;

public:
  Elliptic_Solver(void) {} ;
  ~Elliptic_Solver(void) {} ;

  void Solve(pFun, pFun, unsigned const, unsigned const) ;
  void Save_Mtv(void) ;
} ;

/* Common CLASS METOHODS
 *
 * BASE POLYNOMIAL P2
 *
 *  4
 *  | \
 *  5   3
 *  |     \
 *  0---1---2
 */

// values of shapes function
void
Common::shape(unsigned const nb,
              double const & s, double const & t,
	      double & res) {
  switch ( nb ) {
  case 0: res = (1-2*(s+t))*(1-(s+t)) ; break ;
  case 1: res = 4*s*(1-(s+t))         ; break ;
  case 2: res = s*(2*s-1)             ; break ;
  case 3: res = 4*s*t                 ; break ;
  case 4: res = (2*t-1)*t             ; break ;
  case 5: res = 4*(1-(s+t))*t         ; break ;
  }
}

// values of gradients of bases function
void
Common::shape_grad(unsigned const nb,
                   double const & s,
                   double const & t,
		   double g[2]) {
  switch ( nb ) {
  case 0: g[0] = 4*(s+t)-3     ; g[1] = 4*(s+t)-3     ; break ;
  case 1: g[0] = 4 - 8*s - 4*t ; g[1] = -4*s          ; break ;
  case 2: g[0] = 4*s-1         ; g[1] = 0             ; break ;
  case 3: g[0] = 4*t           ; g[1] = 4*s           ; break ;
  case 4: g[0] = 0             ; g[1] = 4*t-1         ; break ;
  case 5: g[0] = -4*t          ; g[1] = 4 - 4*s - 8*t ; break ;
  }
}

// Vertex CLASS METHODS
inline
unsigned
Vertex::EqNumber(Mesh const & m) const
{ return m . local_number(*this); }

inline
bool
Vertex::IsOnBoundary(Mesh const & m) const
{ return m . local_number(*this) < m . n_bvertex() ; }

// Edge CLASS METHODS
inline
unsigned
Edge::EqNumber(Mesh const & m) const
{ return m . n_vertex() + m . local_number(*this) ; }

inline
bool
Edge::IsOnBoundary(Mesh const & m) const
{ return m . local_number(*this) < m . n_bedge() ; }

// Triangle CLASS METHODS
void
Triangle::eval_JJT(double JJT[2][2], double & detJ) const {
  double iJ[2][2] ;
  inverse_jacobian(0.0, 0.0, iJ) ;
  JJT[0][0] = iJ[0][0]*iJ[0][0] + iJ[0][1]*iJ[0][1] ;
  JJT[0][1] =
  JJT[1][0] = iJ[0][0]*iJ[1][0] + iJ[0][1]*iJ[1][1] ;
  JJT[1][1] = iJ[1][0]*iJ[1][0] + iJ[1][1]*iJ[1][1] ;
  
  detJ = 1/(iJ[0][0] * iJ[1][1] - iJ[1][0] * iJ[0][1]) ;
}

double
Triangle::eval_int_f(unsigned const i,
                     pFun           func,
                     double const & detJ) const {

  static double s[] = { 0.5, 0.5, 0.0 } ;
  static double t[] = { 0.0, 0.5, 0.5 } ;

  double b ;
  double res = 0 ;
  for ( unsigned k = 0 ; k < 3 ; ++k ) {
    shape(i, s[k], t[k], b) ;
    res += func( xm(k), ym(k) ) * b ;
  }
  return detJ * res / 6 ;
}

double
Triangle::eval_int_grad(unsigned const i,
                        unsigned const j,
                        double const JJT[2][2],
                        double const & detJ) const {

  static double s[] = { 0.5, 0.5, 0.0 } ;
  static double t[] = { 0.0, 0.5, 0.5 } ;

  double gi[2], gj[2] ;

  double res = 0 ;
  for ( unsigned k = 0 ; k < 3 ; ++k ) {
    shape_grad(i, s[k], t[k], gi) ;
    shape_grad(j, s[k], t[k], gj) ;
    res += JJT[0][0] * gi[0] * gj[0] +
           JJT[0][1] * gi[0] * gj[1] +
           JJT[1][0] * gi[1] * gj[0] +
           JJT[1][1] * gi[1] * gj[1] ;
  }
  return detJ * res / 6 ;
}

unsigned
Triangle::EqNumber( Mesh const & m, Unsigned const loc) const {
  if ( loc % 2 == 1 ) return edge(loc/2)   . EqNumber(m) ;
  else                return vertex(loc/2) . EqNumber(m) ;
}

bool
Triangle::IsOnBoundary( Mesh const & m, const Unsigned loc) const {
  if ( loc % 2 == 1 ) return edge(loc/2)   . IsOnBoundary(m) ;
  else                return vertex(loc/2) . IsOnBoundary(m) ;
}

// Elliptic_Solver CLASS METHODS
void
Elliptic_Solver::Solve(pFun f, pFun g,
                       unsigned const nx, unsigned const ny) {
  unsigned i, j, k ;

  // build the mesh
  mesh . std_tensor_mesh( nx, ny, NULL, NULL, NULL ) ;
  
  // allocate memory
  unsigned neq  = mesh . n_vertex() + mesh . n_edge() ;
  unsigned nnum = 2*neq + neq * neq ;
  sol = new double [ nnum ] ;
  mat = new double * [ neq ] ;
  if ( sol == NULL || mat == NULL )
    { cerr << "not enought memory" << endl ; exit(0) ; }
  rhs    = sol + neq ;
  mat[0] = rhs + neq ;
  for ( i = 1 ; i < neq ; ++i ) mat[i] = mat[i-1] + neq ;

  // clean up memory
  for ( i = 0 ; i < nnum ; ++i ) sol[i] = 0 ;

  // build the linear system
  Iterator<Triangle> triangle(mesh) ;
  foreach( triangle ) {
    double detJ, JJT[2][2] ;
    triangle -> eval_JJT(JJT, detJ) ;
    
    for ( i = 0 ; i < degree_of_freedom ; ++i ) {
      if ( triangle -> IsOnBoundary(mesh,i) ) continue ;
      unsigned ig = triangle -> EqNumber(mesh,i) ;
      rhs[ig] += triangle -> eval_int_f(i,f,detJ) ;
      for ( j = 0 ; j < degree_of_freedom ; ++j ) {
        unsigned jg = triangle -> EqNumber(mesh,j) ;
        mat[ig][jg] += triangle -> eval_int_grad(i, j, JJT, detJ) ;
      }
    }
  }

  // setup boundary conditions
  Iterator<Vertex> vertex(mesh,1) ;
  foreach( vertex ) {
    unsigned ig = vertex -> EqNumber(mesh) ;
    mat[ig][ig] = 1 ;
    rhs[ig] = g( vertex -> x(), vertex -> y() ) ;
  }

  Iterator<Edge> edge(mesh,1) ;
  foreach( edge ) {
    unsigned ig = edge -> EqNumber(mesh) ;
    mat[ig][ig] = 1 ;
    rhs[ig] = g( edge -> xm(), edge -> ym() ) ;
  }

  // copy rhs to the solution vector
  for ( i = 0 ; i < neq ; ++i ) sol[i] = rhs[i] ;

  // solve the linear system by modified  Gaussian Elimination
  // without pivoting.
  cout << "Solving a " << neq << "x" << neq << " linear system"
       << endl ;
  for ( i = 0 ; i < neq ; ++i ) {
    for ( k = 0 ; k < neq ; ++k ) {
      if ( k != i ) {
        double bf = mat[k][i]/mat[i][i] ;
        sol[k] -= bf * sol[i] ;  
        for ( j = i+1 ; j < neq ; ++j )
          mat[k][j] -= bf * mat[i][j] ;
      }
    }
  }
  for ( i = 0 ; i < neq ; ++i ) sol[i] /= mat[i][i] ;
}

void
Elliptic_Solver::Save_Mtv(void) {
  cout << "saving data file..." ;
  cout . flush() ;
  ofstream file("p2.mtv") ;
  file << "$ DATA=CONTCURVE" << endl
       << "%contstyle=2 meshplot=true topLabel=P2" << endl ;
  Iterator<Triangle> ip(mesh) ;
  foreach ( ip ) {
    for ( unsigned nv = 0 ; nv < 3 ; ++nv ) {
      Vertex & V = ip -> vertex(nv) ;
      unsigned i = mesh . local_number(V) ;
      file << V . x() << " " << V . y() << " " << sol[i] << endl ;
    }
    file << endl ;
  }
  file << "$ END" << endl ;
  file . close() ;
  cout << "saved" << endl ;
}

// Problem Definition

static
double
f(double const &, double const &)
{ return -4 ; }

static
double
g(double const & x, double const & y)
{ return x*x+y*y ; }

int
main() {
  Elliptic_Solver es ;
  es . Solve( f, g, 8, 8) ;
  es . Save_Mtv() ;
}

// end of file p2.cc
// by Enrico Bertolazzi & Gianmarco Manzini
