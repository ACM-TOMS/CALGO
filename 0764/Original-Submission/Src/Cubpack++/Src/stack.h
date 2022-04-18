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
// File : stack.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Stack
// -------------------------
//
// BASECLASSES:
//   ReferenceCounting
//
// PURPOSE:
//   implements a stack. an iterator over all the elements
//   is provided.
//
// TEMPLATES:
//   The type T of the elements to be stacked is
//   templated. T should have T* NewCopy() defined.
//
//
// METHODS:
//   CONSTRUCTORS:
//     1) Stack()
//     -----------
//      default constructor
//
//   SELECTORS:
//     1) T* Top() const
//     -----------------
//     returns a pointer to the top-element of the stack.
//     this element is not removed so modifications of it
//     are potentially dangerous.
//     not to be called when Size() == 0
//
//     2)  Boolean Empty() const
//     --------------------------
//
//     3)  unsigned int Size() const
//     ------------------------------
//
//     4)  Boolean IteratorAtEnd() const
//     ---------------------------------
//     returns True if the iterator has passed the
//     bottom element of the Stack or if Size() == 0
//
//   MODIFIERS:
//     1) T* Pop()
//     -----------
//     gets the top-element of the stack. it is removed.
//     not to be called when Size() == 0
//
//     2) void MakeEmpty()
//     ---------------
//     Empties the Stack
//
//     3) IteratorReset()
//     ------------------
//     a subsequent call of IteratorNext() will return the
//     top of the stack, unless the stack is empty.
//
//     4) T* IteratorNext()
//     --------------------
//     advances the iterator and returns the next element.
//     not to be called when IteratorAtEnd()
//
//     5) void Push(T* t)
//     ----------------------------
//     pushes a pointer to *t on the stack
//
//     6) void Merge(Stack<T>& S)
//     -------------------------------
//     pushes all elements of S on the Stack. The order in
//     which they are pushed is not defined. The state of
//     S afterwards is not defined.
//
//   OPERATORS:
//     None
//   SPECIAL:
//     None
///////////////////////////////////////////////////////////


#ifndef SELEMENT
#define SELEMENT
template< class T>
struct SElement
    {
    SElement<T>* Next;
    T* Contents;
    };
#endif

///////////////////////////////////////////////////

#ifndef STACK_H
#define STACK_H
///////////////////////////////////////////////////

#include <refcount.h>
#include <boolean.h>
#include <iostream.h>
///////////////////////////////////////////////////
template <class T>
class Stack :public ReferenceCounting
  {
  public:

  Stack();
  ~Stack();
  void MakeEmpty();
  Boolean Empty () const;
  void Push (T*);
  T* Pop();
  T* Top()const;
  void Merge (Stack<T>& Source);
  void IteratorReset();
  Boolean IteratorAtEnd() const;
  T* IteratorNext();
  unsigned int Size() const;


  private:

  SElement<T>* TheTop;
  SElement<T>* Current;
  unsigned int Number;


  };
/////////////////////////////////////////////////
#include <templist.h>
#ifdef TEMPLATEINCLUDE
#include <stack.c>
#endif

#endif
