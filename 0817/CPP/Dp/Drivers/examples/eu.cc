/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  P2MESH : a polygon-based mesh manager                                   |
 |           example driver program                                         |
 |                                                                          |
 |  date         : 2000, 8 February                                         |
 |  release date : 1999, 11 October                                         |
 |  version      : 1.1                                                      |
 |  file         : eu.cc                                                    |
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
 |   implementation of the exact Riemann solver for the compressible        |
 |   Euler equations.                                                       |
 |                                                                          |
 |  references:                                                             |
 |    primer.ps  -- beginner introduction with commented examples and       |
 |                  references therein                                      |
 |                                                                          |
\*--------------------------------------------------------------------------*/

# include "eu.hh"

# include <stdlib.h>
# include <math.h>
# include <fstream>
# include <iostream>

using namespace std ;

typedef Euler::Real Real ;

inline Real Abs(Real const & a)
{ return a > 0 ? a : -a ; }

inline Real max(Real const & a, Real const & b)
{ return a > b ? a : b ; }

inline Real min(Real const & a, Real const & b)
{ return a < b ? a : b ; }

static Real const GAMMA = 1.4 ;
static Real const G1 = (GAMMA - 1) / (2 * GAMMA) ;
static Real const G2 = (GAMMA + 1) / (2 * GAMMA) ;
static Real const G3 = 2 * GAMMA / (GAMMA - 1) ;
static Real const G4 = 2 / (GAMMA - 1) ;
static Real const G5 = 2 / (GAMMA + 1) ;
static Real const G6 = (GAMMA - 1) / (GAMMA + 1) ;
static Real const G7 = (GAMMA - 1) / 2 ;
static Real const G8 = 1 / GAMMA ;
static Real const G9 = GAMMA - 1 ;

Real Euler::ec(Real const val[4]) {
  return 0.5 * ( val[1]*val[1] + val[2]*val[2] ) / r(val) ;
}

Real Euler::P(Real const val[4]) {
  Real press = G9 * ( E(val) - ec(val) ) ;
  if ( press <= 0 ) {
    cerr
      << "Euler::P("
      << val[0] << ","
      << val[1] << ","
      << val[2] << ","
      << val[3] << ") found bad pressure p = "
      << press << endl ;
    exit(0) ;
  }
  return press ;
}

Real Euler::C(Real const val[4]) {
  Real C2 = GAMMA * P(val) / r(val) ;
  if ( C2 <= 0 ) {
    cerr
      << "Euler::C("
      << val[0] << ","
      << val[1] << ","
      << val[2] << ","
      << val[3] << ") found bad speed C^2 = "
      << C2 << endl ;
    exit(0) ;
  }
  return sqrt(C2) ;
}

bool Euler::ok_State(Real const val[4]) {
  Real press = G9 * ( E(val) - ec(val) ) ;
  return press > 0 && r(val) > 0 ;
}

void Euler::CFL(Real       & CFL,
		Real const & dt,
		Real const & h,
		Real const   val[4] ) {
  Real c = C(val) ;
  Real u = U(val) ; if ( u < 0 ) u = -u ;
  Real v = V(val) ; if ( v < 0 ) v = -v ;

  CFL = max( dt*(max(u,v)+c)/h, CFL ) ;
}

void Euler::Flux(Real fx[4], Real fy[4], Real const val[4]) {

  Real u = U(val) ;
  Real v = V(val) ;
  
  for ( unsigned i = 0 ; i < 4 ; ++i ) {
    fx[i] = u * val[i] ;
    fy[i] = v * val[i] ;
  }

  Real press = P(val) ;

  fx[1] += press ;
  fx[3] += press * u ;

  fy[2] += press ;
  fy[3] += press * v ;

}

void Euler::LF(Real         nflux[4],
               Real const   lsol[4],
               Real const   rsol[4],
               Real const & nx,
               Real const & ny) {

  Real vl = nx * U(lsol) + ny * V(lsol) ;
  Real vr = nx * U(rsol) + ny * V(rsol) ;
  Real artvisc = max( Abs(vl) + C(lsol), Abs(vr) + C(rsol)) ;
  
  Real lfx[4], lfy[4], rfx[4], rfy[4] ;
  
  Flux(lfx, lfy, lsol) ;
  Flux(rfx, rfy, rsol) ;

  for ( unsigned i = 0 ; i < 4 ; ++i )
    nflux[i] = 0.5 * ( nx * (lfx[i]+rfx[i]) + ny * (lfy[i]+rfy[i])
		       - artvisc * ( rsol[i] - lsol[i] ) ) ;
}

void Euler::StateEval(State      & s,
		      Real const   S[4],
		      Real const & nx,
		      Real const & ny ) {

  s.r = r(S) ;
  s.u = vn(S,nx,ny) ;
  s.p = P(S) ;
  s.c = GAMMA * s.p / r(S) ;
  if ( s.c <= 0 ) {
    cerr << "error in StateEval negative sound speed" << endl ;
    exit(0) ;
  }
  s.c = sqrt( s.c ) ;

}

//----------------------------------------------------------------------

void Euler::Godunov(Real         nflux[4],
                    Real const   lsol[4],
                    Real const   rsol[4],
                    Real const & nx,
                    Real const & ny) {

  if ( !ok_State(lsol) ) {
    cerr << "Euler::Godunov bad left state" << endl ;
    exit(0) ;
  }
  
  if ( !ok_State(rsol) ) {
    cerr << "Euler::Godunov bad right state" << endl ;
    exit(0) ;
  }

  Real const SS = 0 ;
  State sl, sr, sm ;
  StateEval( sl, lsol, nx, ny ) ;
  StateEval( sr, rsol, nx, ny ) ;

  Real PM, UM ;
  // LOCAL RIEMANN PROBLEM RP(I,I+1) IS SOLVED EXACTLY
  Riemann(PM, UM, sl, sr) ;
  // SOLUTION IS SAMPLED AT S=X/T=0 ALONG T-AXIS
  Sample(PM, UM, SS, sl, sr, sm) ;

  Real qn = sm.u ;
  Real qt = qn > 0 ? vt(lsol,nx,ny) : vt(rsol,nx,ny) ;
  Real u  = qn * nx - qt * ny ;
  Real v  = qn * ny + qt * nx ;

  nflux[0] = qn * sm.r ;
  nflux[1] = qn * sm.r * u + sm.p * nx ;
  nflux[2] = qn * sm.r * v + sm.p * ny ;
  nflux[3] = qn * ( sm.p + sm.p / G9 + 0.5 * sm.r *( u*u + v*v ) ) ;
}

//----------------------------------------------------------------------

void Euler::Riemann(Real & P, Real & U, State const & sl, State const & sr) {
  // COMPUTE PRESSURE PM AND PARTICLE VELOCITY UM IN THE MIDDLE
  // PM IS FOUND ITERATIVELY BY A NEWTON-RAPHSON METHOD.

  // COMPUTE GUESS VALUE FROM PVRS RIEMANN SOLVER
  Real PPV = 0.5 * (sl.p + sr.p) - 0.125 * (sr.u - sl.u) * (sl.r + sr.r) * (sl.c + sr.c) ;
  Real PMIN = min(sl.p, sr.p) ;
  Real PMAX = max(sl.p, sr.p) ;
  Real QRAT = PMAX / PMIN  ;

  if ( QRAT <= 2.0 && (PMIN <= PPV && PPV <= PMAX) ) {
    // USE PVRS SOLUTION AS GUESS
    P = PPV ;
  } else {
    if (PPV < PMIN) { // USE TWO-RAREFACTION SOLUTION
      Real PNU = sl.c + sr.c - G7 * (sr.u - sl.u) ;
      Real PDE = sl.c / pow(sl.p, G1 ) + sr.c / pow( sr.p, G1 ) ;
      P = pow(PNU / PDE, G3 ) ;
    } else { // USE TWO-SHOCK APPROXIMATION WITH PPV AS ESTIMATE
      Real GEL = sqrt( (G5 / sl.r) / (G6 * sl.p + PPV) ) ;
      Real GER = sqrt( (G5 / sr.r) / (G6 * sr.p + PPV) ) ;
      P = (GEL * sl.p + GER * sr.p - (sr.u - sl.u) ) / (GEL + GER) ;
    }
   }

  Real const TOL = 1e-6 ;
  Real FL, FR, FLD, FRD ;
  Real P0 = P ;
  Real DU = sr.u - sl.u ;
  for ( unsigned k = 0 ; k < 50 ; ++k ) {
    Prefun(FL, FLD, P, sl) ;
    Prefun(FR, FRD, P, sr) ;
    P -= (FL + FR + DU) / (FLD+FRD) ;
    if ( Abs( (P - P0) / (P + P0) ) <= 0.5 * TOL ) goto fine ;
    P0 = P > 0 ? P : TOL ;
  }
  cout << "Euler::Riemann(...) DIVERGENCE IN NEWTON-RAPHSON ITERATION" << endl ;

fine:
  // COMPUTE U
  U = 0.5 * (sl.u + sr.u + FR - FL) ;
}

void Euler::Prefun(Real & F, Real & FD, Real const & P, State const & s) {
  if (P <= s.p) { // RAREFACTION WAVE
    Real PRAT = P / s.p ;
    F  = G4 * s.c * (pow(PRAT,G1) - 1) ;
    FD = (1.0 / (s.r * s.c) ) * pow(PRAT, -G2 ) ;
  } else { // SHOCK WAVE
    Real AK = G5 / s.r ;
    Real BK = G6 * s.p ;
    Real QRT = sqrt(AK / (BK + P) ) ;
    F  = (P - s.p) * QRT ;
    FD = (1 - 0.5 * (P - s.p) / (BK + P) ) * QRT ;
  }
}

void Euler::Sample(Real const & PM, Real const & UM, Real const & S,
                   State const & sl, State const & sr, State & sm ) {
  if ( S <= UM ) {
    // SAMPLE POINT IS TO THE LEFT OF THE CONTACT
    if ( PM <= sl.p) { // LEFT FAN
      Real SHL = sl.u - sl.c ;
      if ( S <= SHL ) { //LEFT DATA STATE
        sm.r = sl.r ;
        sm.u = sl.u ;
        sm.p = sl.p ;
        sm.c = sl.c ;
      } else {
        Real CML = sl.c * pow(PM / sl.p, G1 ) ;
        Real STL = UM - CML ;
        if ( S > STL ) { // MIDDLE LEFT STATE
          sm.r = sl.r * pow(PM / sl.p, G8 ) ;
          sm.u = UM ;
          sm.p = PM ;
          sm.c = sqrt( GAMMA * sm.p / sm.r ) ;
        } else { // FAN LEFT STATE (INSIDE FAN)
          sm.u = G5 * (sl.c + G7 * sl.u + S)  ;
          sm.c = G5 * (sl.c + G7 * (sl.u - S) )   ;
          sm.r = sl.r * pow(sm.c / sl.c, G4 ) ;
          sm.p = sl.p * pow(sm.c / sl.c, G3 ) ;
        }
      }
    } else { // LEFT SHOCK
       Real PML = PM / sl.p ;
       Real SL = sl.u - sl.c * sqrt(G2 * PML + G1) ;
       if ( S <= SL ) { // LEFT DATA STATE
         sm.r = sl.r ;
         sm.u = sl.u ;
         sm.p = sl.p ;
         sm.c = sl.c ;
       } else { // MIDDLE LEFT STATE (BEHIND SHOCK)
         sm.r = sl.r * (PML + G6 ) / (PML * G6 + 1.0) ;
         sm.u = UM ;
         sm.p = PM ;
         sm.c = sqrt( GAMMA * sm.p / sm.r ) ;
       }
     }
   } else { // RIGHT OF CONTACT
     if ( PM > sr.p ) { // RIGHT SHOCK
       Real PMR = PM / sr.p ;
       Real SR  = sr.u + sr.c * sqrt(G2 * PMR + G1 ) ;
       if (S >= SR) { // RIGHT DATA STATE
         sm.r = sr.r ;
         sm.u = sr.u ;
         sm.p = sr.p ;
         sm.c = sr.c ;
       } else { // MIDDLE RIGHT STATE (BEHIND SHOCK)
         sm.r = sr.r * (PMR + G6 ) / (PMR * G6 + 1.0) ;
         sm.u = UM ;
         sm.p = PM ;
         sm.c = sqrt( GAMMA * sm.p / sm.r ) ;
       }
     } else { // RIGHT FAN
       Real SHR = sr.u + sr.c ;
       if ( S >= SHR ) { // RIGHT DATA STATE
         sm.r = sr.r ;
         sm.u = sr.u ;
         sm.p = sr.p ;
         sm.c = sr.c ;
       } else {
         Real CMR = sr.c * pow(PM / sr.p, G1 ) ;
         Real STR = UM + CMR ;
         if (S <= STR) { // MIDDLE RIGHT STATE
           sm.r = sr.r * pow(PM / sr.p, G8 ) ;
           sm.u = UM ;
           sm.p = PM ;
           sm.c = sqrt( GAMMA * sm.p / sm.r ) ;
         } else { // FAN RIGHT STATE (INSIDE FAN)
           sm.u = G5 * ( - sr.c + G7 * sr.u + S) ;
           sm.c = G5 * (sr.c - G7 * (sr.u - S) ) ;
           sm.r = sr.r * pow(sm.c / sr.c, G4 ) ;
           sm.p = sr.p * pow(sm.c / sr.c, G3 ) ;
         }
       }
     }
   }
}

// end of file eu.cc
// by Enrico Bertolazzi & Gianmarco Manzini
