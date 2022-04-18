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
// File : region.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Region
// -------------------------------------------------
//
// BASECLASSES:
//   None
//
// PURPOSE:
//   base class for all regions, providing storage for
//   integral and error approximations
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) Region()
//     -----------
//
//   SELECTORS:
//     1) real Integral() const
//     -----------------------
//
//     2) real AbsoluteError() const
//     -----------------------------
//
//     3) Boolean Hopeless() const
//     ---------------------------
//
//   MODIFIERS:
//     None
//
//   OPERATORS:
//     1) Boolean operator<(const Region&)const
//     ----------------------------------------
//     compares the absolute errors
//
//     2) Boolean operator>(const Region&)const
//     ----------------------------------------
//     compares the absolute errors
//
//     3) Boolean operator<=(const Region&)const
//     ----------------------------------------
//     compares the absolute errors
//
//     4) Boolean operator>=(const Region&)const
//     ----------------------------------------
//     compares the absolute errors
//
//   SPECIAL:
//     None
//
/////////////////////////////////////////////////////////
#ifndef REGION_H
#define REGION_H
//////////////////////////////////////////////////
#include <boolean.h>
#include <real.h>
#include <reginfo.h>
#include <pointer.h>


/////////////////////////////////////////////////

class Region
  {
  public:

  Region();
  real Integral() const;
  real AbsoluteError() const ;
  Boolean Hopeless()const;
  Boolean operator<(const Region&) const;
  Boolean operator<=(const Region&) const;
  Boolean operator>(const Region&) const;
  Boolean operator>=(const Region&) const;

  protected :

  RegionInfo& LocalInfo();

  private:

  Pointer<RegionInfo> RI_ptr;

  };
///////////////////////////////////////////////
#endif
