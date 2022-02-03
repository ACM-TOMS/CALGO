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
// File : tools.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// PURPOSE : this file contains some simple tools
//////////////////////////////////////////////////

#ifndef TOOLS_H
#define TOOLS_H
#include <real.h>

///////////////////////////////////////////////////////
inline real
min (real a,real b)
  {
  return (a<b) ? a : b ;
  }


inline real
max (real a,real b)
  {
  return (a>b) ? a : b ;
  }

inline real
sign (real x)
  {
  return (x>0) ? 1.0 : -1.0;
  }

///////////////////////////////////////////////////////
#endif
