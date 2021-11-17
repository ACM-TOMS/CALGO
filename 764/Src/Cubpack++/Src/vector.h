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
// File : vector.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Vector
// --------------------------
//
// BASECLASSES:
//   None
//
//
// PURPOSE:
//   implements vectors of arbitrary length.
//
// TEMPLATES:
//   The type T of the elements is templated. T can be any
//   class with a copy constructor, a default constructor
//   and an assignment operator.
//
// METHODS:
//   CONSTRUCTORS:
//     1) Vector()
//     -------------
//     default constructor. the length is not defined so
//     the Vector must be initialized afterwards by an
//     assignment operation. use before this assignment
//     will result in an error.
//
//     2) Vector(const Vector&)
//     -----------------------
//     copy constructor.
//
//     3) Vector(unsigned int n)
//     ---------------------------
//     constructs a vector of length n;
//
//     4) Vector(unsigned int n,T* p)
//     -------------------
//     assumes p is pointing to a classical C array and
//     copies its elements. if the array has not the
//     expected dimension n, an error might occur .
//
//
//   SELECTORS:
//     1) unsigned int Size() const
//     ----------------
//     returns the length of a vector.
//
//   MODIFIERS:
//     None
//   OPERATORS:
//     1) Vector<T>& operator=(const Vector<T>&)
//     -----------------------------------------
//     assignment operator. Sizes must be equal,
//      unless the lefthandside hasn't got
//      a Size.
//
//     2) Boolean operator ==(const Vector&)
//     -------------------------------------
//     equality operator. Sizes must be equal
//
//     3) Boolean operator !=(const Vector&)
//     -------------------------------------
//     !(==). Sizes must be equal
//
//     4) T& operator[](unsigned int) const
//     ---------------------------
//     selects the ith component. 0<=i<Size.
//   SPECIAL:
//     None
///////////////////////////////////////////////////////////


#ifndef VECTOR_H
#define VECTOR_H
/////////////////////////////////////////////////

#include <boolean.h>

/////////////////////////////////////////////////


template <class T>
class Vector
  {

  public:

  Vector();
  Vector(unsigned int length);
  Vector(unsigned  int,T*);
  Vector(const Vector<T>&);
  T& operator[](unsigned int) ;
  const T& operator[](unsigned int) const;
  Vector<T>& operator=(const Vector<T>&);
  Boolean operator == (const Vector<T>&) const;
  Boolean operator != (const Vector<T>&) const;
  ~Vector();
  unsigned int Size() const;



  protected:

  unsigned int TheSize;
  T* Contents;
  };
///////////////////////////////////////////////////
#include <templist.h>
#ifdef TEMPLATEINCLUDE
#include <vector.c>
#endif

#endif
