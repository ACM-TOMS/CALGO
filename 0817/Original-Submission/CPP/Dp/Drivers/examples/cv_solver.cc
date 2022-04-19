/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  P2MESH : a polygon-based mesh manager                                   |
 |           example driver program                                         |
 |                                                                          |
 |  date         : 1999, 11 October                                         |
 |  version      : 1.0                                                      |
 |  file         : cv.cc                                                    |
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
 |   implementation of a first order cell-vertex finite volume method for   |
 |   the compressible Euler equations on a 2D domain.                       |
 |   Data file for the ramp test case are in the file ramp.*                |
 |                                                                          |
 |  references:                                                             |
 |    primer.ps  -- beginner introduction with commented examples and       |
 |                  references therein                                      |
 |                                                                          |
\*--------------------------------------------------------------------------*/

# include "p2mesh.hh"
# include <math.h>

typedef double Real ;

class Vertex ;
class Edge ;
class Triangle ;
class Mesh ;
class Common ;
class Solver ;

class Common : public p2_common<Vertex,Edge,Triangle,Mesh,
                                3,true,Real> {
protected:

  typedef bool (*PCHECK)   (Real const [4]) ;
  
  typedef void (*PFLUX)    (Real [4],
                            Real [4],
                            Real const [4]) ;
                            
  typedef void (*PNUMFLUX) (Real [4],
                            Real const [4],
                            Real const [4],
                            Real const &,
                            Real const &) ;
                            
  typedef void (*PCFL)     (Real &,
                            Real const &,
                            Real const &,
                            Real const [4]) ;

  typedef enum {
    BC_INTERNAL=0,
    BC_SUPERSONIC_INLET,
    BC_SOLID,
    BC_FREE
  } BC ;

} ;

class Vertex : public p2_vertex<Common> {
private:
  friend class Edge ;
  friend class Solver ;
  Real _area, hxy, sol[4], sol0[4] ;
public:
  void Init(Real const[4]) ;  
  Real const & area(void) const { return _area ; }
  void RK_Setsol(void) ;
  void RK_Update(Real const &, Unsigned const) ;
} ;

class Edge : public p2_edge<Common> {
private:
  friend class Vertex ;
  friend class Solver ;
  BC   ibc ;
  Real num_flux[2][4], nx[2], ny[2], len[2] ;
public:
  void Init(void) ;
  void InternalNumFlux(PNUMFLUX) ;
  void BoundaryNumFlux(PNUMFLUX, Real const [4]) ;
} ;

class Triangle : public p2_poly<Common> {} ;

class Mesh : public p2_mesh<Common> {} ;

class Solver : public Common {
private:
  static void mark_edge(Edge &, Unsigned const &) ;

  Mesh               mesh ;
  Iterator<Vertex>   vertex ;
  Iterator<Edge>     edge, iedge, bedge ;
  Iterator<Triangle> triangle ;

  Unsigned max_iter ;
  Real     CFL_run, Tend, time, dt ;
  
  PCHECK   ok_State ;
  PNUMFLUX NumFlux ;
  PFLUX    Flux ;
  PCFL     CFLxy ;

  Real inlet_state[4], init_state[4] ;

public:
  Solver(PFLUX, PNUMFLUX, PCHECK, PCFL);
  void SetUp(char const *) ;
  void SetTimeStep(bool &, Unsigned const);
  void TimeStep(void) ;
  void Save_Mtv(void) ;
} ;

/* * * * * * * * * * * * *\
 * Vertex Class Methods  *
\* * * * * * * * * * * * */

inline
void
Vertex::Init(Real const state[4]) {
  copy(state, state+4, sol) ;
  unsigned i ;
  _area = 0 ;
  for ( i = 0 ; i < n_poly() ; ++i )
    _area += poly(i).area() ;
  _area /= 3 ;

  hxy = edge(0) . length() ;
  for ( i = 1 ; i < n_edge() ; ++i )
    hxy = min(hxy, edge(i) . length() ) ;

}

inline
void
Vertex::RK_Setsol(void) {
  copy(sol, sol+4, sol0) ;
}

void
Vertex::RK_Update(Real const & dt, Unsigned const irk) {

  // compute residual
  Real res[4] ;
  res[0] = res[1] = res[2] = res[3] = 0 ;
  for ( Unsigned ie = 0 ; ie < n_edge() ; ++ie ) {
    Edge & E = edge(ie) ;
    bool ok_dir = this == &E.vertex(0) ;
    Real len = ok_dir ? E.len[0] : -E.len[0] ;
    res[0] += len * E.num_flux[0][0] ;
    res[1] += len * E.num_flux[0][1] ;
    res[2] += len * E.num_flux[0][2] ;
    res[3] += len * E.num_flux[0][3] ;
    
    len = (ok_dir || E.ibc != BC_INTERNAL) ? E.len[1] : -E.len[1] ;
    res[0] += len * E.num_flux[1][0] ;
    res[1] += len * E.num_flux[1][1] ;
    res[2] += len * E.num_flux[1][2] ;
    res[3] += len * E.num_flux[1][3] ;
  }

  // update
  static Real crk0[2] = {1, 0.5} ;
  static Real crk1[2] = {0, 0.5} ;
  static Real CRKR[2] = {1, 0.5} ;
  Real crkr = CRKR[irk]*dt/area() ;

  sol[0] = crk0[irk] * sol0[0] + crk1[irk] * sol[0] - crkr * res[0] ;
  sol[1] = crk0[irk] * sol0[1] + crk1[irk] * sol[1] - crkr * res[1] ;
  sol[2] = crk0[irk] * sol0[2] + crk1[irk] * sol[2] - crkr * res[2] ;
  sol[3] = crk0[irk] * sol0[3] + crk1[irk] * sol[3] - crkr * res[3] ;

}

/* * * * * * * * * * * *\
 * Edge Class Methods  *
\* * * * * * * * * * * */

void
Edge::Init(void) {
  nx[0] = poly(0).yc() - ym() ;
  ny[0] = xm() - poly(0).xc() ;
  
  if ( ok_poly(1) ) {
    nx[1] = ym() - poly(1).yc() ;
    ny[1] = poly(1).xc() - xm() ;
  } else {
    nx[1] = p2_edge<Common>::nx()/2 ;
    ny[1] = p2_edge<Common>::ny()/2 ;
  }

  for ( Unsigned i = 0 ; i < 2 ; ++i ) {
    len[i] = sqrt( nx[i]*nx[i] + ny[i]*ny[i] ) ;
    nx[i] /= len[i] ;
    ny[i] /= len[i] ;
  }
}

void
Edge::InternalNumFlux(PNUMFLUX NumFlux) {
  Vertex & VA = vertex(0) ;
  Vertex & VB = vertex(1) ;
  NumFlux(num_flux[0], VA.sol, VB.sol, nx[0], ny[0]) ;
  NumFlux(num_flux[1], VA.sol, VB.sol, nx[1], ny[1]) ;
}

void
Edge::BoundaryNumFlux(PNUMFLUX NumFlux, Real const inlet[4]) {
  Vertex & VA = vertex(0) ;
  Vertex & VB = vertex(1) ;

  NumFlux(num_flux[0], VA.sol, VB.sol, nx[0], ny[0]) ;

  Real lsol[4], rsol[4] ;
  lsol[0] = 0.5*(VA.sol[0] + VB.sol[0]) ;
  lsol[1] = 0.5*(VA.sol[1] + VB.sol[1]) ;
  lsol[2] = 0.5*(VA.sol[2] + VB.sol[2]) ;
  lsol[3] = 0.5*(VA.sol[3] + VB.sol[3]) ;

  switch (ibc) {
  case BC_FREE:
    copy(lsol, lsol+4, rsol) ;
  break ;
  case BC_SUPERSONIC_INLET:
    copy(inlet, inlet+4, rsol) ;
  break ;
  case BC_SOLID:
    {
      Real qt = -lsol[1] * ny[1] + lsol[2] * nx[1] ;
      Real qn = 0 ;
      rsol[0] = lsol[0] ;
      rsol[1] = qn * nx[1] - qt * ny[1] ;
      rsol[2] = qn * ny[1] + qt * nx[1]  ;
      rsol[3] = lsol[3] ;
    }
  break ;
  default:
    cerr << "bad boundary " << (int)ibc << endl ;
    exit(0) ;
  }

  NumFlux(num_flux[1], lsol, rsol, nx[1], ny[1]) ;
}

/* * * * * * * * * * * * *\
 * Solver Class Methods  *
\* * * * * * * * * * * * */

Solver::Solver(PFLUX    Flux_,
               PNUMFLUX NumFlux_,
               PCHECK   ok_State_,
               PCFL     Cfl_) {
  Flux     = Flux_ ;
  NumFlux  = NumFlux_ ;
  ok_State = ok_State_ ;
  CFLxy    = Cfl_ ;
}

void Solver::mark_edge(Edge & E, Unsigned const & marker) {
  switch ( marker ) {
  case 0 : E.ibc = BC_INTERNAL         ; break ;
  case 1 : E.ibc = BC_SUPERSONIC_INLET ; break ;
  case 2 : E.ibc = BC_SOLID            ; break ;
  case 3 : E.ibc = BC_FREE             ; break ;
  default:
    cerr << "mark_edge( E, " << marker
         << ") bad boundary condition"  << endl ;
    exit(0) ;
  }
}

void Solver::SetUp(char const * file) {

  char file_par[1024] ;
  strcpy(file_par,file) ;
  strcat(file_par,".inp") ;

  ifstream file_input( file_par ) ;

  if ( ! file_input . good() ) {
    cerr << "error in opening file: " << file_par << endl ;
    exit(0) ;
  }

  time = 0 ;
  file_input
    >> dt 
    >> Tend 
    >> max_iter 
    >> CFL_run
    >> inlet_state[0]
    >> inlet_state[1]
    >> inlet_state[2]
    >> inlet_state[3]
    >> init_state[0]
    >> init_state[1]
    >> init_state[2]
    >> init_state[3] ;

  cout
    << "Parameters" << endl
    << "dt       = " << dt       << endl
    << "Tend     = " << Tend     << endl
    << "max_iter = " << max_iter << endl
    << "CFL_run  = " << CFL_run  << endl
    << endl
    << "Input state:" 
    << " r = " << setw(5) << inlet_state[0]
    << " u = " << setw(5) << inlet_state[1]
    << " v = " << setw(5) << inlet_state[2]
    << " E = " << setw(5) << inlet_state[3]
    << endl
    << "Initial state:" 
    << " r = " << setw(5) << init_state[0]
    << " u = " << setw(5) << init_state[1]
    << " v = " << setw(5) << init_state[2]
    << " E = " << setw(5) << init_state[3]
    << endl << endl ;
    
  file_input . close() ;
  
  // initialize
  mesh     . read_mesh(file, NULL, mark_edge, NULL, 1) ;
  vertex   . set_loop(mesh) ;
  edge     . set_loop(mesh) ;
  bedge    . set_loop(mesh,1) ;
  iedge    . set_loop(mesh,2) ;
  triangle . set_loop(mesh) ;
  
  foreach ( vertex ) vertex -> Init(init_state) ;
  foreach ( edge   ) edge   -> Init() ;

}

void
Solver::SetTimeStep(bool & continue_loop, Unsigned const iter) {
  Real CFL_curr = 0 ;
  foreach(vertex)
    CFLxy( CFL_curr, dt, vertex -> hxy, vertex -> sol ) ;

  Real rapp = min(1.2, CFL_run / CFL_curr) ;
  dt       *= rapp ;
  CFL_curr *= rapp ;

  // chek time step
  Real new_time = time+dt ;
  if ( new_time > Tend ) {
    continue_loop = false ;
    dt   = Tend - time ;
    time = Tend ;
  } else {
    time = new_time ;
    continue_loop = continue_loop && iter < max_iter ;
  }

  cout
    << " iter="       << setw(4) << iter
    << " time (n+1)=" << setw(8) << time
    << " CFL="        << setw(8) << CFL_curr
    << " dt="         << setw(8) << dt
    << endl ;

}

void
Solver::TimeStep(void) {
  foreach(vertex) vertex -> RK_Setsol() ;
  for( Unsigned irk = 0 ; irk < 2 ; ++irk ) {
    foreach(iedge) iedge -> InternalNumFlux(NumFlux) ;
    foreach(bedge) bedge -> BoundaryNumFlux(NumFlux,inlet_state) ;
    foreach(vertex) {
      vertex -> RK_Update(dt,irk) ;
      if ( !ok_State(vertex -> sol) ) {
	cerr << "POSITIVITY_CHECK: negative pressure found" ;
	exit(0) ;
      }
    }
  }
}

void
Solver::Save_Mtv(void) {
  ofstream file("cv.mtv") ;
  if ( ! file . good() ) {
    cerr << "Cannot open for write file: ``cv.mtv''" << endl ;
    exit(0) ;
  }

  file << "$ DATA=CONTCURVE\n%contstyle=2 topLabel=mass"
       << endl ;
  foreach ( triangle ) {
    for ( Unsigned nv = 0 ; nv < triangle -> n_vertex() ; ++nv ) {
      Vertex & V = triangle -> vertex(nv) ;
      file << V.x() << " " << V.y() << " " << V.sol[0] << endl ;
    }
    file << endl ;
  }

  file << "$ END" << endl ;
  file . close() ;
}

# include "eu.hh"

int
main() {

  Solver solver(Euler::Flux,
		Euler::Godunov,
		Euler::ok_State,
		Euler::CFL) ;

  solver . SetUp("ramp") ;

  bool continue_loop = true  ;
  for ( unsigned iter = 1 ; continue_loop ; ++iter ) {
    solver.SetTimeStep(continue_loop, iter) ; // variable time step dt
    solver.TimeStep() ; // update one time step
  } ;

  solver . Save_Mtv() ;
  cout << "End of Program" << endl ;
}

// end of file cv.cc
// by Enrico Bertolazzi & Gianmarco Manzini
