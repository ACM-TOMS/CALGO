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
// File : vstack.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS VectorStack
// -------------------------------------------------
//
// BASECLASSES:
//   Vector<T>
//
// PURPOSE:
//   implements a stack of elements T, with individual
//   access possible.
//
// TEMPLATES:
//   The type T of the elements to be stacked is
//   templated. They must have a copy constructor
//
// METHODS:
//   CONSTRUCTORS:
//     1) VectorStack()
//     -----------------
//     creates an empty VectorStack
//
//     2) VectorStack(const VectorStack& V)
//     -------------------------------------
//
//   SELECTORS:
//     None
//
//   MODIFIERS:
//     None
//
//   OPERATORS:
//     1) void operator+=(const T& t) ;
//     ------------------------------
//     pushes a copy of t on top of the stack.
//
//   SPECIAL:
//     None
//
/////////////////////////////////////////////////////////

#ifndef VSTACK_H
#define VSTACK_H
/////////////////////////////////////////////////

#include <boolean.h>
#include <vector.h>

/////////////////////////////////////////////////


template <class T>
class VectorStack: public Vector<T>
  {

  public:

  VectorStack();
  VectorStack(const VectorStack<T>&);
  VectorStack<T>& operator=(const VectorStack<T>&);
  void operator+=(const T&);
  };
///////////////////////////////////////////////////
#include <templist.h>
#ifdef TEMPLATEINCLUDE
#include <vstack.c>
#endif
///////////////////////////////////////////////////////

#endif
