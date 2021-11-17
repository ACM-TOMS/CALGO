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
///////////////////////////////////////////////
//File rule.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   30 Jan 1996     V0.1g(no longer compare signed and unsigned)
/////////////////////////////////////
#include <rule.h>
////////////////////////////////////
template <class GEOMETRY>
Rule<GEOMETRY>::Rule()
  :ReferenceCounting()
   {
   }
////////////////////////////////////////////////
template <class GEOMETRY>
Rule<GEOMETRY>::~Rule()
  {
  }
/////////////////////////////////////////////////
template <class GEOMETRY>
void
Rule<GEOMETRY>::Apply(Integrand& F,GEOMETRY& G, real& i, real&d)
  {
  }
///////////////////////////////////////////////////
template <class GEOMETRY>
void
Rule<GEOMETRY>::ApplyWithDiffs(Integrand& F,GEOMETRY& G, real& r, real&e, Vector<real>& D)
  {
  for (unsigned int i= 0; i<D.Size(); i++)
    {
    D[i] = 0.0;
    };
  Apply(F,G,r,e);
  }
/////////////////////////////////////////////////////
