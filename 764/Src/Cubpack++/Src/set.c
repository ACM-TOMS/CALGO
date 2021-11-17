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
///////////////////////////////////////////////////////
//File set.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//////////////////////////////////////////////////////

#include <iostream.h>
#include <set.h>

///////////////////////////////////////////////////////
template <class T>
Set<T>::Set()
  :ReferenceCounting()
  {
  Number = 0;
  }

////////////////////////////////////////////////////////////
template <class T>
Set<T>::~Set()
  {
  }
///////////////////////////////////////////////////////
template <class T>
unsigned int
Set<T>::Size() const
  {
  return(Number);
  }

////////////////////////////////////////////////////////////
template <class T>
Boolean
Set<T>::Empty()
const
  {
  return (Number == 0)? True : False ;
  }
/////////////////////////////////////////////////////////
