/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  P2MESH : a polygon-based mesh manager                                   |
 |           example driver program                                         |
 |                                                                          |
 |  date         : 1999, 11 October                                         |
 |  version      : 1.0                                                      |
 |  file         : q2.cc                                                    |
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
 |    Q2 conforming base polynomials.                                       |
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
class Quad ;
class Mesh ;
class Elliptic_Solver ;

typedef double (*pFun)(double const &, double const &) ;

// setup common code & data
class Common : public p2_common<Vertex,Edge,Quad,Mesh,4> {
  static double p0(double const &) ;
  static double p1(double const &) ;
  static double p2(double const &) ;

  static double dp0(double const &) ;
  static double dp1(double const &) ;
  static double dp2(double const &) ;
protected:
  static unsigned const degree_of_freedom = 9 ;
  
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

// define Quad class
class Quad : public p2_poly<Common> {
  static double detJ[2][2] ;
  static double JJT[2][2][2][2] ;
  static double st[2] ;
public:
  void   eval_JJT(void) const ;
  double eval_int_f(unsigned const, pFun) const ;
  double eval_int_grad(unsigned const, unsigned const) const ;
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

/* Common CLASS METHODS
 *
 * BASE POLYNOMIAL Q2
 *
 *  6---5---4
 *  |       |
 *  7   8   3
 *  |       |
 *  0---1---2
 */
 
inline double Common::p0(double const & x) { return 0.5*(x-1)*x ; }
inline double Common::p1(double const & x) { return (1-x)*(1+x) ; }
inline double Common::p2(double const & x) { return 0.5*(x+1)*x ; }

inline double Common::dp0(double const & x) { return x-0.5 ; }
inline double Common::dp1(double const & x) { return -2*x ; }
inline double Common::dp2(double const & x) { return x+0.5 ; }

// values of shapes function
void
Common::shape(unsigned const nb,
              double const & s, double const & t,
	      double & res) {
  switch ( nb ) {
  case 0: res = p0(s)*p0(t) ; break ;
  case 1: res = p1(s)*p0(t) ; break ;
  case 2: res = p2(s)*p0(t) ; break ;
  case 3: res = p2(s)*p1(t) ; break ;
  case 4: res = p2(s)*p2(t) ; break ;
  case 5: res = p1(s)*p2(t) ; break ;
  case 6: res = p0(s)*p2(t) ; break ;
  case 7: res = p0(s)*p1(t) ; break ;
  case 8: res = p1(s)*p1(t) ; break ;
  }
}

// values of gradients of shapes function
void
Common::shape_grad(unsigned const nb,
                   double const & s,
                   double const & t,
		   double g[2]) {
  switch ( nb ) {
  case 0: g[0] = dp0(s)*p0(t) ; g[1] = p0(s)*dp0(t) ; break ;
  case 1: g[0] = dp1(s)*p0(t) ; g[1] = p1(s)*dp0(t) ; break ;
  case 2: g[0] = dp2(s)*p0(t) ; g[1] = p2(s)*dp0(t) ; break ;
  case 3: g[0] = dp2(s)*p1(t) ; g[1] = p2(s)*dp1(t) ; break ;
  case 4: g[0] = dp2(s)*p2(t) ; g[1] = p2(s)*dp2(t) ; break ;
  case 5: g[0] = dp1(s)*p2(t) ; g[1] = p1(s)*dp2(t) ; break ;
  case 6: g[0] = dp0(s)*p2(t) ; g[1] = p0(s)*dp2(t) ; break ;
  case 7: g[0] = dp0(s)*p1(t) ; g[1] = p0(s)*dp1(t) ; break ;
  case 8: g[0] = dp1(s)*p1(t) ; g[1] = p1(s)*dp1(t) ; break ;
  }
}

// Vertex CLASS METHODS
inline
unsigned
Vertex::EqNumber(Mesh const & m) const
{ return m . local_number(*this) ; }

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

// Quad CLASS METHODS
void
Quad::eval_JJT() const {

  for ( unsigned i = 0 ; i < 2 ; ++i ) {
    for ( unsigned j  = 0 ; j < 2 ; ++j ) {
      double iJ[2][2] ;
      inverse_jacobian(st[i],st[j], iJ) ;
      JJT[i][j][0][0] = iJ[0][0]*iJ[0][0] + iJ[0][1]*iJ[0][1] ;
      JJT[i][j][0][1] =
      JJT[i][j][1][0] = iJ[0][0]*iJ[1][0] + iJ[0][1]*iJ[1][1] ;
      JJT[i][j][1][1] = iJ[1][0]*iJ[1][0] + iJ[1][1]*iJ[1][1] ;

      detJ[i][j] = 1/(iJ[0][0] * iJ[1][1] - iJ[1][0] * iJ[0][1]) ;
    }
  }
}

double
Quad::eval_int_f(unsigned const i, pFun func) const {

  double x, y, b, res = 0 ;
  for ( unsigned ii = 0 ; ii < 2 ; ++ii ) {
    for ( unsigned jj = 0 ; jj < 2 ; ++jj ) {
      shape(i, st[ii], st[jj], b) ;
      st_to_xy(st[ii], st[jj], x, y) ;
      res += detJ[ii][jj] * func(x,y) * b ;
    }
  }
  return res ;
}

double
Quad::eval_int_grad(unsigned const i, unsigned const j) const {
  double gi[2], gj[2], res = 0 ;
  for ( unsigned ii = 0 ; ii < 2 ; ++ii ) {
    for ( unsigned jj = 0 ; jj < 2 ; ++jj ) {
      shape_grad(i, st[ii], st[jj], gi) ;
      shape_grad(j, st[ii], st[jj], gj) ;
      res += detJ[ii][jj] *
	     ( JJT[ii][jj][0][0] * gi[0] * gj[0] +
	       JJT[ii][jj][0][1] * gi[0] * gj[1] +
	       JJT[ii][jj][1][0] * gi[1] * gj[0] +
	       JJT[ii][jj][1][1] * gi[1] * gj[1] ) ;
    }
  }
  return res ;
}

unsigned
Quad::EqNumber( Mesh const & m, Unsigned const loc) const {
  if ( loc == 8 ) {
    return m . n_vertex() + m . n_edge() + m . local_number(*this) ;
  } else {
    if ( loc % 2 == 1 ) return edge(loc/2)   . EqNumber(m) ;
    else                return vertex(loc/2) . EqNumber(m) ;
  }
}

bool
Quad::IsOnBoundary( Mesh const & m, const Unsigned loc) const {
  if ( loc == 8 ) {
    return false ;
  } else {
    if ( loc % 2 == 1 ) return edge(loc/2)   . IsOnBoundary(m) ;
    else                return vertex(loc/2) . IsOnBoundary(m) ;
  }
}

// Elliptic_Solver CLASS METHODS
void
Elliptic_Solver::Solve(pFun f, pFun g,
                       unsigned const nx, unsigned const ny) {
  unsigned i, j, k ;

  // build the mesh
  mesh . std_tensor_mesh( nx, ny, NULL, NULL, NULL ) ;
  
  // allocate memory
  unsigned neq = mesh . n_vertex() + mesh . n_edge() + mesh . n_poly() ;
  unsigned nnum = 2*neq + neq * neq ;
  sol = new double [ nnum ] ;
  mat = new double * [ neq ] ;
  if ( sol == NULL || mat == NULL ) {
    cerr << "not enought memory" << endl ;
    exit(0) ;
  }
  rhs    = sol + neq ;
  mat[0] = rhs + neq ;
  for ( i = 1 ; i < neq ; ++i ) mat[i] = mat[i-1] + neq ;
  
  // clean up memory
  for ( i = 0 ; i < nnum ; ++i ) sol[i] = 0 ;

  // build the linear system
  Iterator<Quad> quad(mesh) ;
  foreach( quad ) {
    quad -> eval_JJT() ;
    
    for ( i = 0 ; i < degree_of_freedom ; ++i ) {
      if ( quad -> IsOnBoundary(mesh,i) ) continue ;
      unsigned ig = quad -> EqNumber(mesh,i) ;
      rhs[ig] += quad -> eval_int_f(i,f) ;
                           
      for ( j = 0 ; j < degree_of_freedom ; ++j ) {
        unsigned jg = quad -> EqNumber(mesh,j) ;
        mat[ig][jg] += quad -> eval_int_grad(i, j) ;
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
  ofstream file("q2.mtv") ;
  file << "$ DATA=CONTCURVE" << endl
       << "%contstyle=2 meshplot=true topLabel=Q2" << endl ;
  Iterator<Quad> ip(mesh) ;
  foreach ( ip ) {
    for ( unsigned nv = 0 ; nv < 4 ; ++nv ) {
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

// allocate static variables
double Quad::detJ[2][2] ;
double Quad::JJT[2][2][2][2] ;
double Quad::st[2] = { -0.577350269189626, +0.577350269189626 } ;

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

// end of file q2.cc
// by Enrico Bertolazzi & Gianmarco Manzini
