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
//File  E2sec.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   14 Feb 1995     V0.1e(bug fix in transformation)
//////////////////////////////////////////

#include <error.h>
#include <math.h>
#include <E2sec.h>
//////////////////////////////////////////
PlaneSector::PlaneSector(
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
      "Error: trying to construct a plane sector with centre at infinity");
      TheCenter = (A+(A-B)*((C-B)*((C+B)/2-A))/dp);
      };
    TheInnerRadius=(TheCenter-A).Length();
    Point BB(B-TheCenter), CC(C-TheCenter);
    TheSmallAngle=atan2(BB.Y(),BB.X());
    TheBigAngle=atan2(CC.Y(),CC.X());
    if (TheBigAngle<=TheSmallAngle) TheBigAngle += 2*M_PI;
    }
///////////////////////////////////////////////////
PlaneSector::PlaneSector(const Point& O,real r, real theta1,
                         real theta2)
  : Geometry(2)
  {
  TheCenter = O;
  TheInnerRadius = r;
  TheSmallAngle = theta1;
  TheBigAngle = theta2;
  if (TheSmallAngle == TheBigAngle) TheBigAngle += 2*M_PI;
  }
/////////////////////////////////////////////////////
real
PlaneSector::InnerRadius()
const
  {
  return TheInnerRadius;
  }
////////////////////////////////////////////////////
real
PlaneSector::SmallAngle()
const
  {
  return TheSmallAngle;
  }
/////////////////////////////////////////////////////
real
PlaneSector::BigAngle()
const
  {
  return TheBigAngle;
  }
////////////////////////////////////////////////////
const Point&
PlaneSector::Center()
const
  {
  return TheCenter;
  }
////////////////////////////////////////////////////
