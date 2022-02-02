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
/////////////////////////////////////////////////
//File compreg.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//////////////////////////////////////////////////

#include <compreg.h>
#include <error.h>

/////////////////////////////////////////////////
COMPOUND_REGION::COMPOUND_REGION()
  : Region(),TheStatus(Virgin)
  {
  }
//////////////////////////////////////////////////
void
COMPOUND_REGION::Process()
  {
  Error(Hopeless(),"Processing a hopeless COMPOUND_REGION");
  if (TheStatus==Virgin)
   {
   TheStatus= Active;
   Preprocess();
   }
  else
   {
   Improve();
   }
  }
////////////////////////////////////////////////////
COMPOUND_REGION::~COMPOUND_REGION()
  {
  }
/////////////////////////////////////////////////
