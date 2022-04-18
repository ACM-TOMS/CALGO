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
// File : boolean.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
#ifndef BOOLEAN_H
#define BOOLEAN_H

/////////////////////////////////////////////////
#include <iostream.h>
/////////////////////////////////////////////////
enum Boolean {False, True};

extern ostream& operator<<(ostream&,const Boolean&);

//////////////////////////////////////////////////
#endif
