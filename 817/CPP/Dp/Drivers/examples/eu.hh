/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  P2MESH : a polygon-based mesh manager                                   |
 |           example driver program                                         |
 |                                                                          |
 |  date         : 1999, 11 October                                         |
 |  version      : 1.0                                                      |
 |  file         : eu.hh                                                    |
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
 |    header file of the exact Riemann solver class for the compressible    |
 |    Euler equations.                                                      |
 |                                                                          |
 |  references:                                                             |
 |    primer.ps  -- beginner introduction with commented examples and       |
 |                  references therein                                      |
 |                                                                          |
\*--------------------------------------------------------------------------*/

class Euler {
public:
  typedef double Real ;

private:
  typedef struct { Real r, u, c, p ; } State ;
  static void Prefun(Real &, Real &, Real const &, State const &) ;
  static void StateEval(State &, Real const[4], Real const &, Real const &) ;
  static void Riemann(Real &, Real &, State const &, State const &) ;                                         

  static inline Real const & r (Real const val[4]) { return val[0] ; }
  static inline Real const & rU(Real const val[4]) { return val[1] ; }
  static inline Real const & rV(Real const val[4]) { return val[2] ; }
  static inline Real         U (Real const val[4]) { return val[1]/val[0] ; }
  static inline Real         V (Real const val[4]) { return val[2]/val[0] ; }
  static inline Real const & E (Real const val[4]) { return val[3] ; }

  static Real ec(Real const val[4]) ;
  static Real P (Real const val[4]) ;
  static Real C (Real const val[4]) ;

  static inline Real vn(Real const val[4], Real const & nx, Real const & ny)
    { return ( val[1] * nx + val[2] * ny ) / val[0] ; }

  static inline Real vt(Real const val[4], Real const & nx, Real const & ny)
    { return ( val[2] * nx - val[1] * ny ) / val[0] ; }


public:
  static void CFL(Real &, Real const &,  Real const &, Real const [4]) ;
  static bool ok_State(Real const [4]) ;
  static void Sample(Real const &, Real const &, Real const &,
                     State const &, State const &, State &) ;
  
  static void Flux(Real [4], Real [4], Real const [4]) ;
  
  static void LF(Real [4], Real const [4], Real const [4],
                 Real const & , Real const & ) ;

  static void Godunov(Real [4], Real const [4], Real const [4],
                      Real const &, Real const &) ;

} ;

// end of file eu.hh
// by Enrico Bertolazzi & Gianmarco Manzini

