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
//////////////////////////////////////////////
//File region.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//////////////////////////////////////////////

#include <iostream.h>
#include <point.h>
#include <region.h>
#include <error.h>

/////////////////////////////////////////////
real
Region::AbsoluteError() const
  {
  return  RI_ptr->AbsoluteError();
  }

///////////////////////////////////////////////
real
Region::Integral() const
  {
  return RI_ptr->Integral();
  }

///////////////////////////////////////////////
Boolean
Region::operator< (const Region& r)
const
  {
  return  (Boolean)(AbsoluteError() < r.AbsoluteError());
  }
///////////////////////////////////////////////
Boolean
Region::operator<=(const Region& r)
const
  {
  return  (Boolean)(AbsoluteError() <= r.AbsoluteError());
  }
///////////////////////////////////////////////
Boolean
Region::operator>(const Region& r)
const
  {
  return  (Boolean)(AbsoluteError() > r.AbsoluteError());
  }
///////////////////////////////////////////////
Boolean
Region::operator>=(const Region& r)
const
  {
  return  (Boolean)(AbsoluteError() >= r.AbsoluteError());
  }
///////////////////////////////////////////////
RegionInfo&
Region::LocalInfo()
  {
  return *RI_ptr;
  }
////////////////////////////////////////////////
Region::Region()
  :RI_ptr( new RegionInfo)
  {
  }
////////////////////////////////////////////////
Boolean
Region::Hopeless()
const
  {
  return RI_ptr->Hopeless();
  }
////////////////////////////////////////////////
