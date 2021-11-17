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

#ifndef UTILITY_H
#define UTILITY_H

#include <utility>

class noncopyable
{
  protected:
    noncopyable() {}
    ~noncopyable() {}
    
  private:
    noncopyable(const noncopyable &);
    const noncopyable &operator=(const noncopyable &);
}; 
  
template<class T1,class T2,class T3>
struct triplet
{
  typedef T1 first_type;
  typedef T2 second_type;
  typedef T3 third_type;

  first_type  first;
  second_type second;
  third_type  third;

  triplet(): first(), second(), third() {}
  triplet(const first_type &a, const second_type &b, const third_type &c) : first(a), second(b), third(c) {}
  template<class A,class B,class C> triplet(const triplet<A,B,C> &x) : first(x.first), second(x.second), third(x.third()) {}
};

template<class E,class Tr,class A,class B,class C> inline basic_ostream<E,Tr>& operator<<(basic_ostream<E,Tr> &os, const triplet<A,B,C> &x) { return os << "(" << x.first << "," << x.second << "," << x.third << ")"; }
template<class E,class Tr,class A,class B,class C> inline basic_istream<E,Tr>& operator>>(basic_istream<E,Tr> &is,       triplet<A,B,C> &x) { E c; is >> c&&c=='(' && is>>x.first >> c&&c==',' && is>>x.second >> c&&c==',' >> is>>x.third >> c&&c==')'; }


template<class V1,class V2> V1       &first (pair<V1,V2>       &x) { return x.first;  }
template<class V1,class V2> const V1 &first (const pair<V1,V2> &x) { return x.first;  }
template<class V1,class V2> V1       &second(pair<V1,V2>       &x) { return x.second; }
template<class V1,class V2> const V1 &second(const pair<V1,V2> &x) { return x.second; }
template<class V1,class V2> V1       &third (pair<V1,V2>       &x) { return x.third;  }
template<class V1,class V2> const V1 &third (const pair<V1,V2> &x) { return x.third;  }



template<class T1,class T2>
class attach_pair : public pair<T1,T2>
{
  public:
    typedef attach_pair self;
    typedef pair<T1,T2> base;
    typedef typename base::first_type  first_type;
    typedef typename base::second_type second_type;
  
    typedef first_type        value_type;
    typedef value_type       &reference;
    typedef const value_type &const_reference;
    typedef value_type       *pointer;
    typedef const value_type *const_pointer;
    
    typedef second_type attachment_type;
        
  public:
    attach_pair() : base() {}
    attach_pair(const value_type &v, const attachment_type &a) : base(v,a) {}
    explicit attach_pair(const value_type &v) : base(v,attachment_type()) {}
    explicit attach_pair(const base &x) : base(x) {}
    attach_pair(const self &x) : base(x.value(),x.attachment()) {}
    
    reference       value()       { return first; }
    const_reference value() const { return first; }
    operator reference      ()       { return value; }
    operator const_reference() const { return value(); }
    
    second_type       &attachment()       { return second; }
    const second_type &attachment() const { return second; }
};

template<class T1,class T2> inline attach_pair<T1,T2> attach(const T1 &x, const T2 &y) { return attach_pair<T1,T2>(x,y); }


#endif

