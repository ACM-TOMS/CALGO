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
//File stack.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   30 Jan 1996     V0.1g(no longer compare signed and unsigned)
//////////////////////////////////////////////////////

#include <iostream.h>
#include <stack.h>
#include <error.h>

////////////////////////////////////////////////////////////
template <class T>
Stack<T>::Stack ()
  :ReferenceCounting()
  {
  TheTop = new SElement<T>;
  Error(!TheTop,"Stack:Allocation Failed");
  TheTop ->Next =TheTop;
  Number = 0;
  }


////////////////////////////////////////////////////////////
template <class T>
void
Stack<T>::Push(   T* r)
  {
  TheTop->Contents = r;
  SElement<T>* w1 = TheTop;
  SElement<T>* w2 = w1->Next;
  w1->Next = new SElement<T>;
  Error(!w1->Next,"Stack:Allocation Failed");
  w1->Next->Contents = TheTop->Contents;
  w1->Next->Next = w2;
  Number++;
  }

////////////////////////////////////////////////////////////
template <class T>
void
Stack<T>::Merge (Stack<T>& r)
  {
  if(&r == this) return;
  SElement<T> *w = r.TheTop->Next;
  for (unsigned int i =0; i<r.Number;i++)
    {
    Push(w->Contents);
    w->Contents = 0;
    w=w->Next;
    };
  r.MakeEmpty();
  }

////////////////////////////////////////////////////////////
template <class T>
T*
Stack<T>::Pop()
  {
  Error(Number==0,"Attempt to pop empty stack");
  Number--;
  T* R=(TheTop->Next->Contents);
  SElement<T> *w =TheTop->Next;
  TheTop->Next = TheTop->Next->Next;
  delete w;
  return R;
  }
/////////////////////////////////////////////////////
template <class T>
T*
Stack<T>::Top()
const
  {
  Error(Number==0,"Attempt to access top of empty Stack");
  return (TheTop->Next->Contents);
  }
////////////////////////////////////////////////////////////
//template <class T>
//ostream&
//operator<<(ostream& os,const Stack<T>& s)
  //{
  //os << "Stack\n";
  //SElement<T> *w = s.TheTop->Next;
  //for (int i =0; i<s.Number;i++)
    //{
    //os<< *(w->Contents);
    //w=w->Next;
    //};
  //return os;
  //}
////////////////////////////////////////////////////////////
template <class T>
void
Stack<T>::MakeEmpty()
  {
  if (Number == 0) return;
  SElement<T>* w1 = TheTop->Next;
  SElement<T>* w2 = w1 ->Next;
  for (unsigned int i=0 ;i<Number;i++)
    {
    delete w1->Contents;
    delete w1;
    w1 = w2;
    w2 = w1->Next;
    };
  TheTop->Next = TheTop;
  Number = 0;
  }

//////////////////////////////////////////////////////////
template <class T>
Stack<T>::~Stack ()
  {
  MakeEmpty();
  delete TheTop;
  }
//////////////////////////////////////////////////////////
template <class T>
unsigned int
Stack<T>::Size()
const
  {
  return Number;
  }
//////////////////////////////////////////////////////
template <class T>
Boolean
Stack<T>::Empty()
const
  {
  return (Number == 0) ? True : False ;
  }
//////////////////////////////////////////////////////
template <class T>
void
Stack<T>::IteratorReset()
  {
  Current  = TheTop->Next;
  }
////////////////////////////////////////////////////
template <class T>
Boolean
Stack<T>::IteratorAtEnd()
const
  {
  return  (Current == TheTop) ? True : False;
  }
////////////////////////////////////////////////////
template <class T>
T*
Stack<T>::IteratorNext()
  {
  T*r = Current->Contents;
  Current = Current->Next;
  return r;
  }
/////////////////////////////////////////////////////
