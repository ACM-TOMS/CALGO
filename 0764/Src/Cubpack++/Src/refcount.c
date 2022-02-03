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
//File refcount.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Jan 1996     V0.1f(unused parameter removed)
/////////////////////////////////////////////

#include <refcount.h>

//////////////////////////////////////////////
ReferenceCounting::ReferenceCounting()
  {
  numref = 0;
  }
//////////////////////////////////////////////
ReferenceCounting::ReferenceCounting(const ReferenceCounting&)
  {
   //numref isn't copied!
  numref = 0;
  }
/////////////////////////////////////////////
ReferenceCounting&
ReferenceCounting::operator= (const ReferenceCounting&)
   {
   //numref isn't copied!
   return *this;
   }
////////////////////////////////////////////
void
ReferenceCounting::Refer ()
  {
  numref++;
  }
////////////////////////////////////////////
void
ReferenceCounting::UnRefer()
  {
  numref--;
  }
///////////////////////////////////////////////
unsigned int
ReferenceCounting::NumberOfReferences()
const
  {
  return numref;
  }
///////////////////////////////////////////////
