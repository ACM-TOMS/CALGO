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
//File  POLC2.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   14 Feb 1995     V0.1e(bug fix in transformation)
//////////////////////////////////////////

#include <translat.h>
#include <error.h>
#include <C2.h>
#include <math.h>
#include <polC2.h>
//////////////////////////////////////////
#include <polC2prc.h>

////////////////////////////////////////////
PolarRectangle::PolarRectangle(
    const Point& A, const Point& B, const Point& C)
  :Geometry(2)
    {
    if(C==B)
      {
      TheCenter=A;
      }
    else
      {
      real dp = (C-B)*(A-B);
      Error(dp==0,
      "Error: trying to construct a polar rectangle with centre at infinity");
      TheCenter = (A+(A-B)*((C-B)*((C+B)/2-A))/dp);
      };
    TheInnerRadius=(TheCenter-A).Length();
    TheOuterRadius=(TheCenter-B).Length();
    Point BB(B-TheCenter), CC(C-TheCenter);
    TheSmallAngle=atan2(BB.Y(),BB.X());
    TheBigAngle=atan2(CC.Y(),CC.X());
    if (TheBigAngle<=TheSmallAngle) TheBigAngle += 2*M_PI;
    }
///////////////////////////////////////////////////
PolarRectangle::PolarRectangle(const Point& O,
        real r1, real r2, real t1, real t2)
 : Geometry(2)
  {
  TheCenter = O;
  if (r1 > r2)
    {
    TheInnerRadius = r2;
    TheOuterRadius = r1;
    }
  else
    {
    TheInnerRadius = r1;
    TheOuterRadius = r2;
    };
  TheSmallAngle = t1;
  TheBigAngle = t2;
  }
/////////////////////////////////////////////////////


real
PolarRectangle::InnerRadius()
const
  {
  return TheInnerRadius;
  }
////////////////////////////////////////////////////
real
PolarRectangle::OuterRadius()
const
  {
  return TheOuterRadius;
  }
///////////////////////////////////////////////////
real
PolarRectangle::SmallAngle()
const
  {
  return TheSmallAngle;
  }
/////////////////////////////////////////////////////
real
PolarRectangle::BigAngle()
const
  {
  return TheBigAngle;
  }
////////////////////////////////////////////////////
const Point&
PolarRectangle::Center()
const
  {
  return TheCenter;
  }
////////////////////////////////////////////////////
//Processor<PolarRectangle>*
//PolarRectangle::DefaultProcessor()
//const
  //{
  //return new PolarRectangle_Processor;
  //}
////////////////////////////////////////////////////
