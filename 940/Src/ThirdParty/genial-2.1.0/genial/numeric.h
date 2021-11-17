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

#ifndef NUMERIC2_H
#define NUMERIC2_H

#include <numeric>
#include <utility>
#include <iterator>
#include <cassert>

using namespace std;

// Controls the unrolling size for some calculation loops. Possible values: 1,2,3,4,5,6,7,8
// WARNING: EXPERIMENTAL. Anyway, ICL automatically unrolls these loops...
#ifndef UNROLL_SIZE
#define UNROLL_SIZE 1
#endif


template<class T> inline bool is_even(const T &x) { return !is_odd(x); }
template<class T> inline bool is_odd (const T &x) { return x&1; }

template<class InIt,class T>
pair<T,T> square_accumulate(InIt begin, InIt end, T sum=iterator_traits<InIt>::value_type(0), T sum2=iterator_traits<InIt>::value_type(0))
{
  for ( ; begin!=end; ++begin)
  {
    typename iterator_traits<InIt>::const_reference r = *begin;
    sum+=r;
    sum2+=r*r;
  }
  return pair<T,T>(sum,sum2);
}



template<class T>
T gcd(T a,T b)
{
  T c;
  while (b!=0) { c=a%b; a=b; b=c; }
  return a;
}

template<class T>
T lcm(T a,T b)
{
  return (a/gcd(a,b))*b;
}


template<class InIt, class T>
T accumulate_n(int n, InIt begin, T init)
{
  assert(n>=0);

#if UNROLL_SIZE>1
  for (; n>=UNROLL_SIZE; n-=UNROLL_SIZE, begin+=UNROLL_SIZE)
  {
    init+=*begin;
#if UNROLL_SIZE>1
    init+=*(begin+1);
#endif
#if UNROLL_SIZE>2
    init+=*(begin+2);
#endif
#if UNROLL_SIZE>3
    init+=*(begin+3);
#endif
#if UNROLL_SIZE>4
    init+=*(begin+4);
#endif
#if UNROLL_SIZE>5
    init+=*(begin+5);
#endif
#if UNROLL_SIZE>6
    init+=*(begin+6);
#endif
#if UNROLL_SIZE>7
    init+=*(begin+7);
#endif
  }
#endif

  for (; n>0; --n, ++begin)
    init+=*begin;

  return init;
}

template<int N,class RanIt,class T>
inline T accumulate_n(RanIt begin, T init)
{
  if (N>32) return accumulate_n(N,begin,init);
  if (N> 0) init+=*begin;
  if (N> 1) init+=*(begin+ 1);
  if (N> 2) init+=*(begin+ 2);
  if (N> 3) init+=*(begin+ 3);
  if (N> 4) init+=*(begin+ 4);
  if (N> 5) init+=*(begin+ 5);
  if (N> 6) init+=*(begin+ 6);
  if (N> 7) init+=*(begin+ 7);
  if (N> 8) init+=*(begin+ 8);
  if (N> 9) init+=*(begin+ 9);
  if (N>10) init+=*(begin+10);
  if (N>11) init+=*(begin+11);
  if (N>12) init+=*(begin+12);
  if (N>13) init+=*(begin+13);
  if (N>14) init+=*(begin+14);
  if (N>15) init+=*(begin+15);
  if (N>16) init+=*(begin+16);
  if (N>17) init+=*(begin+17);
  if (N>18) init+=*(begin+18);
  if (N>19) init+=*(begin+19);
  if (N>20) init+=*(begin+20);
  if (N>21) init+=*(begin+21);
  if (N>22) init+=*(begin+22);
  if (N>23) init+=*(begin+23);
  if (N>24) init+=*(begin+24);
  if (N>25) init+=*(begin+25);
  if (N>26) init+=*(begin+26);
  if (N>27) init+=*(begin+27);
  if (N>28) init+=*(begin+28);
  if (N>29) init+=*(begin+29);
  if (N>30) init+=*(begin+30);
  if (N>31) init+=*(begin+31);
  return init;
};


template<class RanIt1, class RanIt2, class T>
T inner_product_n(int n, RanIt1 begin1, RanIt2 begin2, T init)
{
  assert(n>=0);

#if UNROLL_SIZE>1
  for (; n>=UNROLL_SIZE; n-=UNROLL_SIZE, begin1+=UNROLL_SIZE, begin2+=UNROLL_SIZE)
  {
    init+=*begin1 * *begin2;
#if UNROLL_SIZE>1
    init+=*(begin1+1) * *(begin2+1);
#endif
#if UNROLL_SIZE>2
    init+=*(begin1+2) * *(begin2+2);
#endif
#if UNROLL_SIZE>3
    init+=*(begin1+3) * *(begin2+3);
#endif
#if UNROLL_SIZE>4
    init+=*(begin1+4) * *(begin2+4);
#endif
#if UNROLL_SIZE>=5
    init+=*(begin1+5) * *(begin2+5);
#endif
#if UNROLL_SIZE>6
    init+=*(begin1+6) * *(begin2+6);
#endif
#if UNROLL_SIZE>7
    init+=*(begin1+7) * *(begin2+7);
#endif
  }
#endif

  for (; n>0; --n, ++begin1, ++begin2)
    init+=*begin1* *begin2;

  return init;
}


template<class RanIt1, class RanIt2>
RanIt2 copyn(int n, RanIt1 begin1, RanIt2 begin2)
{
  assert(n>=0);

#if UNROLL_SIZE>1
  for (; n>=UNROLL_SIZE; n-=UNROLL_SIZE, begin1+=UNROLL_SIZE, begin2+=UNROLL_SIZE)
  {
    *begin2=*begin1;
#if UNROLL_SIZE>1
    *(begin2+1)=*(begin1+1);
#endif
#if UNROLL_SIZE>2
    *(begin2+2)=*(begin1+2);
#endif
#if UNROLL_SIZE>3
    *(begin2+3)=*(begin1+3);
#endif
#if UNROLL_SIZE>4
    *(begin2+4)=*(begin1+4);
#endif
#if UNROLL_SIZE>5
    *(begin2+5)=*(begin1+5);
#endif
#if UNROLL_SIZE>6
    *(begin2+6)=*(begin1+6);
#endif
#if UNROLL_SIZE>7
    *(begin2+7)=*(begin1+7);
#endif
  }
#endif

  for (; n>0; --n, ++begin1, ++begin2)
    *begin2=*begin1;
  return begin2;
}



template<class RanIt>
typename iterator_traits<RanIt>::value_type min_n(int n, RanIt begin)
{
  assert(n>=1);

  typename iterator_traits<RanIt>::value_type init=*begin;
  --n; ++begin;

#if UNROLL_SIZE>1
  for (; n>=UNROLL_SIZE; n-=UNROLL_SIZE, begin+=UNROLL_SIZE)
  {
    init=min(init,*begin);
#if UNROLL_SIZE>1
    init=min(init,*(begin+1));
#endif
#if UNROLL_SIZE>2
    init=min(init,*(begin+2));
#endif
#if UNROLL_SIZE>3
    init=min(init,*(begin+3));
#endif
#if UNROLL_SIZE>4
    init=min(init,*(begin+4));
#endif
#if UNROLL_SIZE>5
    init=min(init,*(begin+5));
#endif
#if UNROLL_SIZE>6
    init=min(init,*(begin+6));
#endif
#if UNROLL_SIZE>7
    init=min(init,*(begin+7));
#endif
  }
#endif

  for (; n>0; --n,++begin)
    init=min(init,*begin);

  return init;
}


template<class RanIt>
typename iterator_traits<RanIt>::value_type max_n(int n, RanIt begin)
{
  assert(n>=1);

  typename iterator_traits<RanIt>::value_type init=*begin;
  --n; ++begin;

#if UNROLL_SIZE>1
  for (; n>=UNROLL_SIZE; n-=UNROLL_SIZE, begin+=UNROLL_SIZE)
  {
    init=max(init,*begin);
#if UNROLL_SIZE>1
    init=max(init,*(begin+1));
#endif
#if UNROLL_SIZE>2
    init=max(init,*(begin+2));
#endif
#if UNROLL_SIZE>3
    init=max(init,*(begin+3));
#endif
#if UNROLL_SIZE>4
    init=max(init,*(begin+4));
#endif
#if UNROLL_SIZE>5
    init=max(init,*(begin+5));
#endif
#if UNROLL_SIZE>6
    init=max(init,*(begin+6));
#endif
#if UNROLL_SIZE>7
    init=max(init,*(begin+7));
#endif
  }
#endif

  for (; n>0; --n,++begin)
    init=max(init,*begin);

  return init;
}



template<class RanIt1, class RanIt2>
RanIt2 swap_ranges_n(int n, RanIt1 begin1, RanIt2 begin2)
{
  assert(n>=0);

#if UNROLL_SIZE>1
  for (; n>=UNROLL_SIZE; n-=UNROLL_SIZE, begin1+=UNROLL_SIZE, begin2+=UNROLL_SIZE)
  {
    swap(*begin1,*begin2);
#if UNROLL_SIZE>1
    swap(*(begin1+1),*(begin2+1));
#endif
#if UNROLL_SIZE>2
    swap(*(begin1+2),*(begin2+2));
#endif
#if UNROLL_SIZE>3
    swap(*(begin1+3),*(begin2+3));
#endif
#if UNROLL_SIZE>4
    swap(*(begin1+4),*(begin2+4));
#endif
#if UNROLL_SIZE>5
    swap(*(begin1+5),*(begin2+5));
#endif
#if UNROLL_SIZE>6
    swap(*(begin1+6),*(begin2+6));
#endif
#if UNROLL_SIZE>7
    swap(*(begin1+7),*(begin2+7));
#endif
  }
#endif

  for (; n>0; --n, ++begin1, ++begin2)
    swap(*begin1,*begin2);
  return begin2;
}

template<int N,class RanIt1,class RanIt2>
inline RanIt2 swap_ranges_n(RanIt1 begin1, RanIt2 begin2)
{
  if (N>32) return swap_ranges_n(N,begin1,begin2);
  if (N> 0) swap(*begin1     ,*begin2);
  if (N> 1) swap(*(begin1+ 1),*(begin2+ 1));
  if (N> 2) swap(*(begin1+ 2),*(begin2+ 2));
  if (N> 3) swap(*(begin1+ 3),*(begin2+ 3));
  if (N> 4) swap(*(begin1+ 4),*(begin2+ 4));
  if (N> 5) swap(*(begin1+ 5),*(begin2+ 5));
  if (N> 6) swap(*(begin1+ 6),*(begin2+ 6));
  if (N> 7) swap(*(begin1+ 7),*(begin2+ 7));
  if (N> 8) swap(*(begin1+ 8),*(begin2+ 8));
  if (N> 9) swap(*(begin1+ 9),*(begin2+ 9));
  if (N>10) swap(*(begin1+10),*(begin2+10));
  if (N>11) swap(*(begin1+11),*(begin2+11));
  if (N>12) swap(*(begin1+12),*(begin2+12));
  if (N>13) swap(*(begin1+13),*(begin2+13));
  if (N>14) swap(*(begin1+14),*(begin2+14));
  if (N>15) swap(*(begin1+15),*(begin2+15));
  if (N>16) swap(*(begin1+16),*(begin2+16));
  if (N>17) swap(*(begin1+17),*(begin2+17));
  if (N>18) swap(*(begin1+18),*(begin2+18));
  if (N>19) swap(*(begin1+19),*(begin2+19));
  if (N>20) swap(*(begin1+20),*(begin2+20));
  if (N>21) swap(*(begin1+21),*(begin2+21));
  if (N>22) swap(*(begin1+22),*(begin2+22));
  if (N>23) swap(*(begin1+23),*(begin2+23));
  if (N>24) swap(*(begin1+24),*(begin2+24));
  if (N>25) swap(*(begin1+25),*(begin2+25));
  if (N>26) swap(*(begin1+26),*(begin2+26));
  if (N>27) swap(*(begin1+27),*(begin2+27));
  if (N>28) swap(*(begin1+28),*(begin2+28));
  if (N>29) swap(*(begin1+29),*(begin2+29));
  if (N>30) swap(*(begin1+30),*(begin2+30));
  if (N>31) swap(*(begin1+31),*(begin2+31));
  return begin2+N;
};



#endif

