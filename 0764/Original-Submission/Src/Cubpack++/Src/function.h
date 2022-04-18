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
// File : function.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// PURPOSE:
//   mapping from a Point to a real
/////////////////////////////////////////////////////////
#ifndef FUNCTION_H
#define FUNCTION_H
/////////////////////////////////////////////////////////
#include <real.h>
#include <point.h>
/////////////////////////////////////////////////////////

typedef real (*Function)(const Point&);
/////////////////////////////////////////////////////////

#endif
