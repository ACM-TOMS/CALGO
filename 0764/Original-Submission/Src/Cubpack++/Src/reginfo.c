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
//File reginfo.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Nov 1994     V0.1d(re-ordering of member initializers)
/////////////////////////////////////////////////////////
#include <reginfo.h>
////////////////////////////////////////////////////////
RegionInfo::RegionInfo()
  :ReferenceCounting(),
   TheIntegral(0),
   TheAbsoluteError(0),
   IsHopeless(False)
  {
  }
/////////////////////////////////////////////////////////
Boolean
RegionInfo::Hopeless()
const
  {
  return IsHopeless;
  }
////////////////////////////////////////////////////////////
Boolean&
RegionInfo::Hopeless()
  {
  return IsHopeless;
  }
/////////////////////////////////////////////////////////////
real
RegionInfo::Integral()
const
  {
  return TheIntegral;
  }
//////////////////////////////////////////////////////////
real&
RegionInfo::Integral()
  {
  return TheIntegral;
  }
////////////////////////////////////////////////////////
real
RegionInfo::AbsoluteError()
const
  {
  return TheAbsoluteError;
  }
///////////////////////////////////////////////////////////
real&
RegionInfo::AbsoluteError()
  {
  return TheAbsoluteError;
  }
////////////////////////////////////////////////////////////
