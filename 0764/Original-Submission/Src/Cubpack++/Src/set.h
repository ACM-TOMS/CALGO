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
// File : set.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Set
// -----------------------
//
// BASECLASSES:
//   ReferenceCounting
//
//
// PURPOSE:
//   provides a pure virtual base class for all
//   collection classes.
//
// TEMPLATES:
//   the type T of the Set elements is templated. it
//   should have a method T* NewCopy () defined.
//
// METHODS:
//   CONSTRUCTORS:
//     1) Set()
//     --------
//      Default constructor
//
//   SELECTORS:
//     1) virtual T* Look() =0
//     ------------------------
//     returns a pointer to an element of the Set.
//
//     2) int NumberOfElements() const
//     --------------------------------
//     returns the number of elements in the set.
//
//   MODIFIERS:
//     1) virtual T* Get()=0
//     -----------------------
//     returns a pointer to an element of the Set and
//     removes them.
//
//     2) virtual void Clear()=0
//     ---------------------------
//     removes all of the elements
//
//   OPERATORS:
//     1) virtual void operator+=(T*) =0
//     -------------------------------------------
//     adds a reference to an element to the set
//
//   SPECIAL:
//     None
///////////////////////////////////////////////////////////



#ifndef SET_H
#define SET_H
#include <templist.h>
#include <refcount.h>
#include <boolean.h>
template <class T>
class Set: public ReferenceCounting
  {
  public:

  Set();
  virtual ~Set();
  virtual void Clear()=0;
  virtual T* Get () =0;
  virtual T* Look() =0;
  virtual void operator+=( T*) =0;
  unsigned int Size() const;
  Boolean Empty() const;

  protected:

  int Number;
  };
#ifdef TEMPLATEINCLUDE
#include <set.c>
#endif

#endif
