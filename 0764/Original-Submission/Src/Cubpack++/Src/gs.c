
/////////////////////////////////////////////////////////
//                                                     //
//    Cubpack++                                        //
//                                                     //
//        A Package For Automatic Cubature             //
//                                                     //
//        Authors : Ronald Cools                       //
//                  Dirk Laurie                        //
//                  Luc Pluym                          //
//                                                     //
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
// File : gs.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
#include <gs.h>
////////////////////////////////////////////////////////
GeneralizedSector::GeneralizedSector(real (*f)(real) ,real a,real b,const Point& C)
  :Geometry(2),TheCenter(C),TheAlpha(a),TheBeta(b),TheBoundary(f)
  {
  }
/////////////////////////////////////////////////////////
real
GeneralizedSector::Alpha()
const
  {
  return TheAlpha;
  }
//////////////////////////////////////////////////////////
real
GeneralizedSector::Beta()
const
  {
  return TheBeta;
  }
//////////////////////////////////////////////////////////
const Point&
GeneralizedSector::Center()
const
  {
  return TheCenter;
  }
//////////////////////////////////////////////////////////
real
GeneralizedSector::Boundary(real p) const
  {
  GeneralizedSector* G= (GeneralizedSector*) this;
  return G->TheBoundary(p);
  }
///////////////////////////////////////////////////////////
RealFunction
GeneralizedSector::Boundary()
const
  {
  return TheBoundary;
  }
/////////////////////////////////////////////////////////////
