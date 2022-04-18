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
// File : regcoll.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Jan 1996     V0.1f(typedef introduced)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS REGION_COLLECTION
// -------------------------------------------------
//
// BASECLASSES:
//   COMPOUND_REGION
//
// PURPOSE:
//   groups regions and their eventual offspring. applies
//   the global adaptive algorithm.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//   1) REGION_COLLECTION()
//   ----------------------
//
//   2) REGION_COLLECTION(const REGION_COLLECTION&)
//   ----------------------------------------------
//
//   SELECTORS:
//     None
//
//   MODIFIERS:
//    1) LocalIntegrand(Integrand* I)
//    ------------------------------
//    sets the local integrand of all members to I
//    only to be called before Process().
//
//    2) LocalIntegrand(Function F)
//    ------------------------------
//    sets the local integrand of all members to F
//    only to be called before Process().
//
//   OPERATORS:
//    1) REGION_COLLECTION operator+(const COMPOUND_REGION&)
//    ----------------------------------------------------
//
//    2) REGION_COLLECTION& operator+=(const COMPOUND_REGION&)
//    ----------------------------------------------------
//
//   SPECIAL:
//     None
//
/////////////////////////////////////////////////////////
#ifndef REGCOLL_H
#define REGCOLL_H
/////////////////////////////////////////////
#include <compreg.h>
#include <integran.h>
#include <function.h>
#include <pointer.h>
#include <stack.h>
#include <atomreg.h>

/////////////////////////////////////////////

class REGION_COLLECTION : public COMPOUND_REGION
  {
  public:
  typedef Stack < COMPOUND_REGION> StackCOMPOUND_REGION;

  void LocalIntegrand(Integrand*);
  void LocalIntegrand(Function);
  REGION_COLLECTION& operator+=(const COMPOUND_REGION&) ;
  REGION_COLLECTION operator+(const COMPOUND_REGION&) ;
  REGION_COLLECTION();
  REGION_COLLECTION(const REGION_COLLECTION&);

  protected:

  Pointer < StackCOMPOUND_REGION > SCR_ptr;
  void Preprocess();
  void Improve();
  real MaxAtomicError()const;
  COMPOUND_REGION* NewCopy()const;

  };
///////////////////////////////////////////////
#endif
