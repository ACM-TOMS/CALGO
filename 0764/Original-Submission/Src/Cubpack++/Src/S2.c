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
//////////////////////////////////////////
//File S2.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//////////////////////////////////////////

#include <S2.h>
#include <S2rule13.h>
#include <S2adapt.h>
#include <error.h>
#include <C2.h>
#include <polC2.h>
#include <stack.h>
#include <pointer.h>
#include <math.h>

////////////////////////////////////////////
Circle::Circle(const Point& Center, real Radius)
  : Geometry(2),
    TheCenter(Center),
    TheRadius(fabs(Radius))
    {}
///////////////////////////////////////////////
Circle::Circle(const Point& Center, const Point& Boundary)
  : Geometry(2),
    TheCenter(Center),
    TheRadius(sqrt((Center.X()-Boundary.X())*(Center.X()-Boundary.X())
                  +(Center.Y()-Boundary.Y())*(Center.Y()-Boundary.Y())))
    {}
///////////////////////////////////////////////
real
Circle::Volume()
const
  {
  return M_PI*TheRadius*TheRadius;
  }
//////////////////////////////////////////////
const Point&
Circle::Center()
const
  {
  return TheCenter;
  }
/////////////////////////////////////////////
real
Circle::Radius()
const
  {
  return TheRadius;
  }
//////////////////////////////////////////////
//Processor<Circle>*
//Circle::DefaultProcessor()
//const
  //{
  //return new CircleAdaptive(new Circle_Rule13);
  //}
/////////////////////////////////////////////////
