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
//File pointer.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   11 Sep 1994     V0.1b(bug removed from destructor)
//   14 Sep 1994     V0.1c(inlines removed)
//   25 Jan 1996     V0.1f(code moved from .c to .h)
/////////////////////////////////////////////////////////

#include <pointer.h>

////////////////////////////////////////////////
template <class T>
Pointer<T>::Pointer()
   {
   ptr =  0;
   }
////////////////////////////////////////////////
template <class T>
Pointer<T>::Pointer(const Pointer<T>& p)
  {
  ptr = p.ptr;
  if (ptr !=0) ptr->Refer();
  }
///////////////////////////////////////////////
template <class T>
Pointer<T>::Pointer(T* p)
  {
  ptr = p;
  if (ptr !=0) ptr->Refer();
  }
////////////////////////////////////////////////
template <class T>
Pointer<T>&
Pointer<T>::operator=(const Pointer<T>&  p)
  {
  if (ptr !=0)
    {
    ptr->UnRefer();
    if(ptr->NumberOfReferences() == 0 && (ptr != p.ptr))
      {
      delete ptr;
      };
    };
  ptr = p.ptr;
  if (ptr !=0) ptr->Refer();
  return *this;
  }
////////////////////////////////////////////////
template <class T>
Pointer<T>&
Pointer<T>::operator=(T*  p)
  {
  if (ptr !=0)
    {
    ptr->UnRefer();
    if(ptr->NumberOfReferences() == 0  && (p != ptr))
      {
      delete ptr;
      };
    };
  ptr = p;
  if (ptr !=0) ptr->Refer();
  return *this;
  }
////////////////////////////////////////////////
////////////////////////////////////////////////
template <class T>
T&
Pointer<T>::operator[] (int i)
const
  {
  return ptr[i];
  }
///////////////////////////////////////////////
template <class T>
Pointer<T>::~Pointer()
  {
  if (ptr ==0) return;
  ptr->UnRefer();
  if(ptr->NumberOfReferences() == 0)
    {
    delete ptr;
    };
  }
///////////////////////////////////////////////
template <class T>
Boolean
Pointer<T>::operator == (const Pointer<T>& tp)
const
  {
  return (ptr == tp.ptr) ? True : False;
  }
////////////////////////////////////////////////
template <class T>
Boolean
Pointer<T>::operator != (const Pointer<T>& tp)
const
  {
  return (ptr != tp.ptr)? True : False;
  }
////////////////////////////////////////////////
template <class T>
T&
Pointer<T>::operator *()
const
  {
  return *ptr;
  }
////////////////////////////////////////////////
template <class T>
T*
Pointer<T>::operator ->()
const
  {
  return ptr;
  }
////////////////////////////////////////////////
//template <class T>
//Pointer<T>::operator T*()
//  {
//  return ptr;
//  }
////////////////////////////////////////////////
