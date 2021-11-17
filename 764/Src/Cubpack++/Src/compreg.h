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
// File : compreg.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS COMPOUND_REGION
// -------------------------------------------------
//
// BASECLASSES:
//   Region
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
//     1) COMPOUND_REGION()
//     -------------------
//     an empty compound region is created. the integral
//     and error over it will be 0. calls of Process()
//     will have no effect.
//
//     2) COMPOUNDERGION(const COMPOUND_REGION& C)
//     ------------------------------------------
//     implements sharing semantics !!
//
//   SELECTORS:
//     None
//
//   MODIFIERS:
//     1) void LocalIntegrand (Integrand* I)
//     ------------------------------------
//     establishes *I as the integrand to be integrated
//     on all subregions that are contained in *this.
//     When subsequent calls are made, only the last one
//     has effect. It's an error to call this function
//     after the first call of Process().
//
//
//     2) void Process()
//     -----------------
//     performs one step in the global adaptive algorithm.
//     the first call of Process() will compute an initial
//     approximation; subsequent calls will try to improve
//     this
//
//   OPERATORS:
//
//     1) COMPOUND_REGION& operator=
//             (const COMPOUND_REGION& C)
//     ---------------------------------
//     implements sharing semantics !!
//
//   SPECIAL:
//     None
//
//   NOTE: The other public members are not to be
//         used by normal clients.
/////////////////////////////////////////////////////////
#ifndef COMPREG_H
#define COMPREG_H
/////////////////////////////////////////////
#include <region.h>
#include <integran.h>
#include <function.h>
#include <pointer.h>
#include <atomreg.h>

/////////////////////////////////////////////
/////////////////////////////////////////////

class COMPOUND_REGION : public Region
  {
  public:

  virtual void LocalIntegrand(Integrand*)=0;
  void Process();
  COMPOUND_REGION();
  virtual ~COMPOUND_REGION();
  virtual void Preprocess()=0;
  virtual void Improve()=0;
  virtual real MaxAtomicError()const=0;
  virtual COMPOUND_REGION* NewCopy()const =0;


  private:

  enum Status {Virgin,Active};
  Status TheStatus;
  };
///////////////////////////////////////////////
#endif
