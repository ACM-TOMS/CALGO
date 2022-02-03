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
// File : heap.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Heap
// ------------------------
//
// BASECLASSES:
//   Set<T>
//
//
// PURPOSE:
//   implements a heap. elements can be added and the
//   largest element can be accessed.
//
// TEMPLATES:
//   The type  T of the elements is templated. T must
//   provide a T* NewCopy() as well as relational
//   operators <,>,<= and >=
//
// METHODS:
//   CONSTRUCTORS:
//     1) Heap()
//     ---------
//     default constructor
//   SELECTORS:
//     1) T* Look()
//     ------------
//     returns a pointer to the largest element. this
//     element is not removed from the heap. No copy
//     is made, so modifications of the element are
//     potentially dangerous.
//   MODIFIERS:
//     1) T* Get()
//     -----------
//     returns a pointer to the largest element. this
//     element is removed from the heap.
//
//     2) void Clear()
//     --------------
//     makes the heap empty and destroys all of its
//     elements
//
//   OPERATORS:
//     1) void operator+=(T* R)
//     ------------------------------
//      Puts a pointer to *R into the heap.
//
//     2) void operator+=(Heap<T>& H)
//     ------------------------------
//     puts all elements of H in the heap. no copies are
//     made and the state of H is not defined afterwards.
//
//   SPECIAL:
//     None
//   NOTE:
//     Heap uses a class SubHeap in its implementation.
//     This class should be thought of as being local to
//     Heap. However the current compiler we used did not
//     support nested classes within templates. Direct use
//     of SubHeap is not recommended, since it will
//     eventually become local to Heap.
///////////////////////////////////////////////////////////
#ifndef HEAP_H
#define HEAP_H
#include <templist.h>
//////////////////////////////////////////////
#include <set.h>
#include <iostream.h>
#include <boolean.h>
/////////////////////////////////////////////
template <class T,int CAPACITY>
class SubHeap : public Set<T>
  {
  public:

  SubHeap();
  ~SubHeap();
  T* Get();
  T* Look() ;
  T* Swap( T*);
  T* Bottom() ;
  void Clear();
  void operator +=( T*);
  Boolean Saturated() const;
  void Print() const  ;
  int LeftChild(int) const;
  int RightChild(int) const;
  int FatherOfChild(int) const;

  private:

#ifdef SOLARIS
  T* Contents[256];
  SubHeap<T,255>* Children[256];
#else
  T* Contents[CAPACITY+1];
  SubHeap<T,CAPACITY>* Children[CAPACITY+1];
#endif
  int ActiveChild;
  int LastChild;
  };

////////////////////////////////////////////


template <class T>
class Heap : public Set<T>
  {
  public:

  Heap();
  ~Heap();
  T* Get();
  T* Look() ;
  void Clear();
  void operator+= ( T*);
  void operator+= (Heap<T>&);
  void Print() const;


  private:

  SubHeap<T,255> TheSubHeap;
  };
////////////////////////////////////////////
#ifdef TEMPLATEINCLUDE
#include <heap.c>
#endif

#endif
