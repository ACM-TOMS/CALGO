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
/////////////////////////////////////////////////////
//File invert.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
////////////////////////////////////////////////////////
#include <invert.h>

///////////////////////////////////////////////////////
Invert::Invert(Circle* C)
  :Transformation(),
   C_ptr(C),
   RadiusSq(C->Radius()*C->Radius())
   {
   }
//////////////////////////////////////////
void
Invert::Transform (real &w, Point& p)
  {
  const Point diff = p-C_ptr->Center();
  real r=diff.Length();
  real ratio=RadiusSq/(r*r);
  p = C_ptr->Center()+diff*ratio;
  w *= ratio*ratio;
  }
/////////////////////////////////////////////
