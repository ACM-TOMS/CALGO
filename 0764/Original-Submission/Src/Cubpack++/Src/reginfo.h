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
// File : reginfo.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS  RegionInfo
// -------------------------------------------------
//
// BASECLASSES:
//   ReferenceCounting
//
// PURPOSE:
//   to store integral and error of a specific region
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) RegionInfo()
//     ---------------
//     after construction Integral() and AbsoluteError()
//     will both return zero
//
//   SELECTORS:
//     1) real Integral() const
//      ------------------------
//
//     2) real AbsoluteError() const
//     -----------------------------

//     3) Boolean Hopeless() const
//     ---------------------------
//
//   MODIFIERS:
//     1) real& Integral()
//     --------------------
//
//     2) real& AbsoluteError()
//     ------------------------
//
//     3) Boolean& Hopeless()
//     ----------------------
//
//   OPERATORS:
//     None
//
//   SPECIAL:
//     None
//
/////////////////////////////////////////////////////////
#ifndef REGINFO_H
#define REGINFO_H
//////////////////////////////////////////
#include <refcount.h>
#include <real.h>
#include <boolean.h>
//////////////////////////////////////////

class RegionInfo : public ReferenceCounting
  {
  public :

  RegionInfo();
  real Integral()const;
  real& Integral();
  real AbsoluteError()const;
  real& AbsoluteError();
  Boolean Hopeless() const;
  Boolean& Hopeless() ;

  private:

  real TheIntegral;
  real TheAbsoluteError;
  Boolean IsHopeless;
  };
/////////////////////////////////////////
#endif
