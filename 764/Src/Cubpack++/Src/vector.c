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
//File vector.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////
#include <vector.h>
#include <math.h>
#include <error.h>
#include <iostream.h>
////////////////////////////////////////////////////
template<class T>
Vector<T>::Vector(unsigned int s)
  {
  TheSize =s;
  if (s>0)
    {
    Contents = new T[TheSize];
    Error(!Contents,"Vector:allocation failed");
    };
  }

//////////////////////////////////////////////////////
template<class T>
Vector<T>::Vector()
  {
  TheSize = 0;
  }
//////////////////////////////////////////////////////


template<class T>
Vector<T>::~Vector()
  {
  if (TheSize>0) delete [] Contents;
  }

//////////////////////////////////////////////////////

template<class T>
Vector<T>::Vector(const Vector<T>& v)
  {
  TheSize = v.TheSize;
  if (TheSize ==0)
    {
    Contents = 0;
    return;
    };
  Contents = new T[TheSize];
  Error(!Contents,"Vectorcopyconst:allocation failed");
  for (unsigned int i =0;i<TheSize;i++)
    {
    Contents[i] = v.Contents[i];
    };
  }

//////////////////////////////////////////////////////

template<class T>
Vector<T>::Vector( unsigned int s,T* v)
  {
  TheSize = s;
  Contents = new T[TheSize];
  for (unsigned int i =0;i<TheSize;i++)
    {
    Contents[i] = v[i];
    };
  }

//////////////////////////////////////////////////////

template<class T>
const T&
Vector<T>::operator[](unsigned int i)
const
  {
  Error( i>=TheSize,"arrayindex constraint violation");
  return Contents[i];
  }

//////////////////////////////////////////////////////

template<class T>
T&
Vector<T>::operator[](unsigned int i)
  {
  Error( i>=TheSize,"array index constraint violation");
  return Contents[i];
  }


//////////////////////////////////////////////////////

template<class T>
Vector<T>&
Vector<T>::operator=(const Vector<T>& v)
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
Boolean
Vector<T>::operator==(const Vector<T>& v)
const
  {
  if (TheSize !=v.TheSize)
    {
    return False;
    }
  else
    {
    Boolean b = True;
    for (unsigned int i=0;i<TheSize;i++)
      {
      b=(Boolean)(b && (Contents[i]==v.Contents[i]));
      };
      return b;
    };
  }
////////////////////////////////////////////////////
template <class T>
Boolean
Vector<T>::operator!=(const Vector<T>& v)
const
  {
  return (Boolean)!( (*this) == v);
  }
///////////////////////////////////////////////////
template <class T>
unsigned int
Vector<T>::Size()
const
  {
  return TheSize;
  }
//////////////////////////////////////////////
//template <class T>
//ostream&
//operator << (ostream& os,const Vector<T>& v)
  //{
  //os << "Vector";
  //for (unsigned int i=0;i<v.TheSize;i++)
    //{
    //os << v.Contents [i]<<' ';
    //};
  //os<<endl;
  //return os;
  //}
//////////////////////////////////////////////////
