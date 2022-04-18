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
/////////////////////////////////////////////
//File vstack.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////
#include <vstack.h>
#include <stdlib.h>
#include <error.h>
#include <iostream.h>
#include <string.h>

template <class T>
VectorStack<T>::VectorStack()
  :Vector<T>()
 {
 }
///////////////////////////////////////////////////
template <class T>
VectorStack<T>::VectorStack(const VectorStack<T>& v)
  :Vector<T>(v)
  {
  }
///////////////////////////////////////////////////
template <class T>
VectorStack<T>&
VectorStack<T>::operator=(const VectorStack<T>& v)
  {
  if (TheSize == 0)
    {
    TheSize = v.TheSize;
    Contents = new T[TheSize];
    Error(!Contents,"Vectorassign:allocation failed");
    };
  Error( (TheSize != v.TheSize),"lengths of arrays incompatible");
  for (unsigned int i =0;i<TheSize;i++)
    {
    Contents[i] = v.Contents[i];
    };
  return *this;
  }

////////////////////////////////////////////////////
template <class T>
void
VectorStack<T>::operator += (const T& t)
  {
  T* New  = new T[TheSize+1];
  for(unsigned int i=0;i<TheSize;i++)
    {
    New[i+1] = Contents[i];
    };
  if (TheSize !=0)
    {
    delete []  Contents;
    };
  Contents = New;
  Contents [0] = t;
  TheSize++;
  }
///////////////////////////////////////////////////
