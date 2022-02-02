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
// File : pointer
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   14 Sep 1994     V0.1c(inlines removed)
//   25 Jan 1996     V0.1f(code moved from .c to .h)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Pointer
// -------------------------------------------------
//
// BASECLASSES:
//   None
//
// PURPOSE:
//   a device for storing pointers, managing reference
//   counts to the pointed objects and controlling the
//   deletion of the pointed objects . Once a T*
//   has been stored as a Pointer<T> you can use it
//   as if C++ had garbage collection.
//   some compilers have problems converting
//   Pointers to normal pointers.
//
// TEMPLATES:
//   a Pointer<T> stores a T*. T should be derived from
//   ReferenceCounting
//
/////////////////////////////////////////////////////////

#ifndef POINTER_H
#define POINTER_H

/////////////////////////////////////////////////////
#include <boolean.h>
/////////////////////////////////////////////////////
template  <class T>
class Pointer
  {
  public:

  Pointer();
  Pointer(const Pointer<T>&);
  Pointer(T* tp);
  Pointer<T>& operator=(const Pointer<T>& tp);
  Pointer<T>& operator=(T* tp);
  Boolean operator==(const Pointer<T>& tp) const;
  Boolean operator!=(const Pointer<T>& tp) const;
  T& operator* () const;
  T* operator->() const;
  T& operator [](int i) const;
  operator T*(){return ptr;}
  ~Pointer();

  private:

  T* ptr;
  };

/////////////////////////////////////////////////////
#include <templist.h>
#ifdef  TEMPLATEINCLUDE
#include <pointer.c>
#endif
/////////////////////////////////////////////////////
#endif
