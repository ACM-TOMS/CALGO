//GENIAL - GENeric Image & Array Library
//Copyright (C) 2005  IENT - RWTH Aachen
//
//This program is free software; you can redistribute it and/or
//modify it under the terms of the GNU General Public License
//as published by the Free Software Foundation; either version 2
//of the License, or (at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program; if not, write to the Free Software
//Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#ifndef BITSTREAM_H
#define BITSTREAM_H

#include <iostream>
#include <cassert>

//namespace genial
//{

using namespace std;

typedef unsigned char BYTE;
typedef unsigned short WORD;

template<class E,class Tr> inline basic_istream<E,Tr> &operator>(basic_istream<E,Tr> &is,       BYTE &x) { return is.get((char &)x); }
template<class E,class Tr> inline basic_ostream<E,Tr> &operator<(basic_ostream<E,Tr> &os, const BYTE &x) { return os.put((char &)x); }

inline istream &operator>(istream &is, WORD &x)
{
  BYTE h,l;
  is>h;
  is>l;
  x = (WORD(h)<<8) + l;
  return is;
}

template<class E,class Tr> inline basic_ostream<E,Tr> &operator<(basic_ostream<E,Tr> &os, const WORD &x)
{
  os<(BYTE)(x>>8);
  os<(BYTE)(x);
  return os;
}

template<class T>
class bits
{
  public:
    typedef bits self;

    typedef T value_type;
    typedef int size_type;

  private:
    size_type l;
    value_type v;

  public:
    bits() : l(0), v() {}
    bits(size_type n) : l(n), v() {}
    bits(size_type n, const value_type &x) : l(n), v(x) {}
    bits(const self &r) : l(r.l), v(r.v) {}    
    
    template<class U> self &operator=(const bits<U> &r) { l=r.size(); v=r.v; }

    operator const value_type &() const { return v; }

    size_type size() const { return l; }

    void resize(size_type n) { l=n; }

    self  operator<< (size_type pos) const { return self(size()+pos, v<<pos); }
    self &operator<<=(size_type pos)       { resize(size()+pos); v<<=pos; return *this; }

    self  operator>> (size_type pos) const { if (size()<=pos) return bits(); return self(size()-pos, v>>pos); }
    self &operator>>=(size_type pos)       { if (size()<=pos) return bits(); resize(size()-pos); v>>=pos; return *this; }

    self  operator| (value_type x) const { return self(v|x, size()); }
    self &operator|=(value_type x) const { v|=x; return *this; }

    template<class U> bool operator<(const bits<U> &x) const { return v<x.v;  }

    template<class U> void push_back(const bits<U> &x) { (*this)<<=x.size(); v|=x; }
    void push_back(char x) { (*this)<<=8*sizeof(char); v|=x; }
    
    template<class U> void pop_front(bits<U> &x) { x=(*this)>>(size()-x.size()); resize(size()-x.size()); }
    void pop_front(char &x) { x=(*this)>>(size()-8*sizeof(char)); resize(size()-8*sizeof(char)); }
};

template<class E, class T, class N>
basic_istream<E,T> &operator>(basic_istream<E,T> &is, bits<N> &x)
{
  typedef typename basic_istream<E,T>::char_type char_type;
  
  bits<N> t;
  while (t.size()<x.size()) { char_type c; is>c; t.push_back(c); }
  t.resize(x.size());
  x=t;
  return is;
}

template<class E, class T, class N>
basic_ostream<E,T> &operator<(basic_ostream<E,T> &os, const bits<N> &x)
{
  typedef typename basic_istream<E,T>::char_type char_type; 

  bits<N> t=x;
  t.resize(8*((x.size()+(8*sizeof(char_type))-1)/(8*sizeof(char_type))));

  while (t.size()>0) { char_type c; t.pop_front(c);  os < c; } 
    
  return os;   
}
 
template<class E, class T, class N1, class N2>
basic_istream<E,T> &bitread(basic_istream<E,T> &is, bits<N1> &x, bits<N2> &buf)
{
  typedef typename basic_istream<E,T>::char_type char_type;

  while (buf.size()<x.size()) { char_type c; is>>c; buf.push_back(c); }

  x = buf>>(buf.size()-x.size());
  buf.resize(buf.size()-x.size()); 

  return is;
}

template<class E, class T, class N1, class N2>
basic_ostream<E,T> &bitwrite(basic_ostream<E,T> &os, bits<N1> &x, bits<N2> &buf)
{
  typedef typename basic_istream<E,T>::char_type char_type;
 
  buf.push_back(x);
  
  os < (buf>>(buf.size()%(8*sizeof(char_type))));
  buf.resize(buf.size()%(8*sizeof(char_type)));

  return os;
}

//}

#endif

