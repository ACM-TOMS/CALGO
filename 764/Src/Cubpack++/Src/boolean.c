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
////////////////////////////////////////////////
//File boolean.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
////////////////////////////////////////////////

#include <boolean.h>
#include <iostream.h>

ostream& operator << (ostream& os, const Boolean& b)
  {
  if (b==False)
    {
    os<< "False ";
    }
  else
    {
    os<< "True ";
    };
  return os;
  }
