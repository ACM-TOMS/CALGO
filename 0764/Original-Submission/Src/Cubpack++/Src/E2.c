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
//File E2.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Nov 1994     V0.1d(re-ordering of member initializers)
/////////////////////////////////////////////////////////

#include <E2.h>
#include <E2adapt.h>
#include <invert.h>
#include <S2.h>
/////////////////////////////////////////////////////////

Plane::Plane ()
  : Geometry(2),
    xscale(1.0),
    yscale(1.0),
    TheCenter(0,0)
  {
  }

/////////////////////////////////////////////////////////

Plane::Plane (const Point& center)
  : Geometry(2),
    xscale(1.0),
    yscale(1.0),
    TheCenter(center)
  {
  }
/////////////////////////////////////////////////////////
Plane::Plane(const Point& center, real x, real y)
  : Geometry(2),
    xscale(x),
    yscale(y),
    TheCenter(center)
  {
  }
/////////////////////////////////////////////////////////
real
Plane::ScaleX()
const
   {
    return xscale;
   }
/////////////////////////////////////////////////////////
real
Plane::ScaleY()
const
   {
   return yscale;
   }
/////////////////////////////////////////////////////////
//Processor<Plane>*
//Plane::DefaultProcessor()
//const
  //{
  //return new PlaneAdaptive;
  //}
/////////////////////////////////////////////////////////
const Point&
Plane::Center()
const
  {
  return TheCenter;
  }
/////////////////////////////////////////////////////////
