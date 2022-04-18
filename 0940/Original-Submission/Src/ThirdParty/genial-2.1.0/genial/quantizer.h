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

#ifndef QUANTIZER_H
#define QUANTIZER_H

#include <functional>
#include <numeric>


//namespace genial
//{

using namespace std;

template<class Q, class Arg, class Res>
class quantizer_function : public unary_value_function<Arg, Res>
{
  public:
    typedef Q quantizer_type;
    typedef typename quantizer_type::argument_type argument_type;
    typedef typename quantizer_type::value_type value_type;

  protected:
    const quantizer_type &quant;

  protected:
    inline quantizer_function(const quantizer_type &q) : quant(q) {}
    
    inline quantizer_type       &quantizer()       { return quant; }
    inline const quantizer_type &quantizer() const { return quant; }
};

template<class Q>
class quantizer_index_function : public quantizer_function<Q, typename Q::argument_type, typename Q::value_type>
{
  public:
    typedef quantizer_index_function self;
    typedef quantizer_function<Q, typename Q::argument_type, typename Q::value_type> base;
    
    typedef typename base::quantizer_type  quantizer_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;

  public:
    inline quantizer_index_function(const quantizer_type &q) : base(q) {}

    template<class Arg1> inline const_reference operator()(const Arg1 &x) const { return quant(x); } 
};

template<class T, class Q>
struct quantizerIndexArray : public array_value_traits<const T, typename Q::value_type>
{
  typedef const Q quantizer_type;
  typedef array_value_traits<const T, typename Q::value_type> base;
  typedef typename base::array_type array_type;
  typedef quantizer_index_function<quantizer_type> function_type;
  typedef typename UnaryFunctionArray<array_type,function_type>::self self;
};


template<class Q>
class quantizer_value_function : public quantizer_function<Q, typename Q::value_type, typename Q::argument_type>
{
  public:
    typedef quantizer_value_function self;
    typedef quantizer_function<Q, typename Q::value_type, typename Q::argument_type> base;

    typedef typename base::quantizer_type  quantizer_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    
    using base::quantizer;

  public:
    inline quantizer_value_function(const quantizer_type &q) : base(q) {}

    template<class Arg1> inline const_reference operator()(const Arg1 &x) const { return quantizer()[x]; } 
};

template<class T, class Q>
struct quantizerValueArray : public array_value_traits<const T, typename Q::argument_type>
{
  typedef const Q quantizer_type;
  typedef array_value_traits<const T, typename Q::argument_type> base;
  typedef typename base::array_type array_type;
  typedef quantizer_value_function<quantizer_type> function_type;
  typedef typename UnaryFunctionArray<array_type,function_type>::self self;
};

template<class Q>
class quantizer_quantize_function : public quantizer_function<Q, typename Q::argument_type, typename Q::argument_type>
{
  public:
    typedef quantizer_quantize_function self;
    typedef quantizer_function<Q, typename Q::argument_type, typename Q::argument_type> base;
    
    typedef typename base::quantizer_type quantizer_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    
    using base::quantizer;

  public:
    inline quantizer_quantize_function(const quantizer_type &q) : base(q) {}

    template<class Arg1> inline const_reference operator()(const Arg1 &x) const { return quantizer().quantize(x); } 
};

template<class T, class Q>
struct quantizerQuantizeArray : public array_value_traits<const T, typename Q::argument_type>
{
  typedef const Q quantizer_type;
  typedef array_value_traits<const T, typename Q::argument_type> base;
  typedef typename base::array_type array_type;
  typedef quantizer_quantize_function<quantizer_type> function_type;
  typedef typename UnaryFunctionArray<array_type,function_type>::self self;
};

template<class V, class I, I lev, V a, V b>
class fastScalarQuantizer
{
  public:
    typedef V argument_type;
    typedef I value_type;

  public:
    inline fastScalarQuantizer() {}

    inline operator value_type() const { return lev; }    
};

template<class V, class I, I lev, V a, V b, I delta=(b-lev-a+1)/lev+1>
class fastUniformSQ : public fastScalarQuantizer<V,I,lev,a,b>
{
  public:
    typedef fastUniformSQ self;
    typedef fastScalarQuantizer<V,I,lev,a,b> base;
    
    typedef typename base::value_type    value_type;
    typedef typename base::argument_type argument_type;
    
  public:
    inline fastUniformSQ() {}

    inline value_type    operator()(const argument_type &x) const { return (x-a) / delta; }
    inline argument_type operator[](const value_type &i) const { return delta*i+a+delta/2; }

    inline argument_type quantize(const argument_type &x) const { return (*this)[(*this)(x)]; }
};

template<class V, class I=int>
class ScalarQuantizer
{ 
  public:
    typedef V argument_type;
    typedef I value_type;
    
  protected:
    value_type lev;
    argument_type a;
    argument_type b;
    
  public:
    inline ScalarQuantizer(const value_type &l, const argument_type &x, const argument_type &y) : lev(l), a(x), b(y) {}
    
    inline operator value_type() const { return size(); }
    inline value_type size() const { return lev; }

    inline bool has_value(const argument_type   x) const { return (a<=x) && (x<=b); } // no ref because binder do not accept it
    inline bool has_index(const value_type n) const { return (n<size()); }
};


template<class V, class I=int>
class UniformSQ : public ScalarQuantizer<V,I>
{
  public:
    typedef UniformSQ self;
    typedef ScalarQuantizer<V,I> base;
    
    typedef typename base::value_type    value_type;
    typedef typename base::argument_type argument_type;
    
    using base::a;

  private:
    argument_type delta;
    
  public:
    inline UniformSQ(const value_type &l, const argument_type &x, const argument_type &y) : base(l,x,y), delta((y-l-x+1)/l+1) {}
  
    inline value_type  operator()(const argument_type &x) const { return value_type((x-a) / delta); }
    inline argument_type operator[](const value_type  &n) const { return delta*n+a+delta/2; }
    inline argument_type quantize  (const argument_type &x) const { return (*this)[(*this)(x)]; }      

    template<class G> inline typename quantizerIndexArray   <const Matrix<G>, const self>::self operator()(const Matrix<G> &X) const { return quantizerIndexArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Matrix<G>, const self>::self operator[](const Matrix<G> &X) const { return quantizerValueArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Matrix<G>, const self>::self quantize  (const Matrix<G> &X) const { return quantizerQuantizeArray<const Matrix<G>, const self>::self(X, *this); }
};

template<class V, class S=int>
class UniformQuantizer
{
  public:
    typedef UniformQuantizer self;

    typedef V argument_type;
    typedef S value_type;

  private:
    argument_type delta;
    
  public:
    inline UniformQuantizer(const argument_type &d) : delta(d) {}
  
    inline value_type  operator()(const argument_type &x) const { return value_type(x/delta); }
    inline argument_type operator[](const value_type  &n) const { return delta*n+delta/2; }
    inline argument_type quantize  (const argument_type &x) const { return (*this)[(*this)(x)]; }      

    template<class G> inline typename quantizerIndexArray   <const Matrix<G>, const self>::self operator()(const Matrix<G> &X) const { return quantizerIndexArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Matrix<G>, const self>::self operator[](const Matrix<G> &X) const { return quantizerValueArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Matrix<G>, const self>::self quantize  (const Matrix<G> &X) const { return quantizerQuantizeArray<const Matrix<G>, const self>::self(X, *this); }
};


template<class Q>
class MultiQuantizer
{
  public:
    typedef MultiQuantizer self;

    typedef Q quantizer_type;
    typedef typename quantizer_type::argument_type argument_type;
    typedef typename quantizer_type::value_type value_type;

    typedef vector<quantizer_type> vector_type;
    typedef typename vector_type::const_iterator const_iterator;

  private:
    vector_type vector;

  public:
    inline MultiQuantizer() : vector() {}
    inline MultiQuantizer(const quantizer_type &q1                                                    ) : vector() { push_back(q1);                               }
    inline MultiQuantizer(const quantizer_type &q1, const quantizer_type &q2                          ) : vector() { push_back(q1); push_back(q2);                }
    inline MultiQuantizer(const quantizer_type &q1, const quantizer_type &q2, const quantizer_type &q3) : vector() { push_back(q1); push_back(q2); push_back(q3); }
    inline MultiQuantizer(const vector_type &r) : vector(r) {}
    inline MultiQuantizer(const self &r) : vector(r.vector) {}
 
    inline vector_type       &quantizers()       { return vector; }
    inline const vector_type &quantizers() const { return vector; }

    inline const_iterator begin() const { return quantizers().begin(); }
    inline const_iterator end  () const { return quantizers().end  (); }

    inline const quantizer_type &back() const { return quantizers().back(); }

    inline void push_back(const quantizer_type &q) { quantizers().push_back(q); }

    inline value_type size() const { value_type n=0; for (const_iterator p=begin(); p!=end(); ++p) n+=(*p).size(); return n; }     

    inline value_type  operator()(const argument_type &x) const { value_type n=0; for (const_iterator p=begin(); p!=end(); ++p) if ((*p).has_value(x)) return n+(*p)(x); else n+=(*p).size(); return n; assert(0); }     
    inline argument_type operator[](      value_type   n) const {                for (const_iterator p=begin(); p!=end(); ++p) if ((*p).has_index(n)) return (*p)[n];   else n-=(*p).size(); return back()[back().size()]; assert(0); }
    inline argument_type quantize  (const argument_type &x) const { for (const_iterator p=begin(); p!=end(); ++p) if ((*p).has_value(x)) return (*p).quantize(x); return back()[back().size()]; }
//    inline argument_type   quantize  (const argument_type   &x) const { return (*find_if(begin(), end(), bind2nd(mem_fun_ref(&quantizer_type::has_value),x))).quantize(x); }
};

template<class Q, class V=typename Q::argument_type, class I=typename Q::value_type>
class quantizer_quantizer
{
  public:
    typedef quantizer_quantizer self;

    typedef Q quantizer_type;
    typedef V argument_type;
    typedef I value_type;

  private:
    quantizer_type quant;

  public:
    inline quantizer_quantizer() : quant() {}
    inline quantizer_quantizer(const quantizer_type &q) : quant(q) {}
    inline quantizer_quantizer(const self &r) : quant(r.quant) {}
    template<class A> inline quantizer_quantizer(const A &a) : quant(a) {}
    template<class A,class B> inline quantizer_quantizer(const A &a, const B &b) : quant(a,b) {}
    template<class A,class B,class C> inline quantizer_quantizer(const A &a, const B &b, const C &c) : quant(a,b,c) {}
    template<class A,class B,class C,class D> inline quantizer_quantizer(const A &a, const B &b, const C &c, const D &d) : quant(a,b,c,d) {}
    template<class A,class B,class C,class D,class E> inline quantizer_quantizer(const A &a, const B &b, const C &c, const D &d, const E &e) : quant(a,b,c,d,e) {}

    inline quantizer_type       &quantizer()       { return quant; }
    inline const quantizer_type &quantizer() const { return quant; }

    inline value_type size() const { return quantizer().size(); }
};

template<class Q, class V=typename Q::argument_type, class I=typename Q::value_type>
class SymmetryQuantizer : public quantizer_quantizer<Q,V,I>
{
  public:
    typedef SymmetryQuantizer self;
    typedef quantizer_quantizer<Q,V,I> base;
    
    typedef typename base::quantizer_type quantizer_type;
    typedef typename base::value_type     value_type;
    typedef typename base::argument_type  argument_type;
    
    using base::quantizer;
   
  public:
    inline SymmetryQuantizer() : base() {}
    inline SymmetryQuantizer(const quantizer_type &q) : base(q) {}
    inline SymmetryQuantizer(const self &r) : base(static_cast<const base &>(r)) {}
    template<class A> inline SymmetryQuantizer(const A &a) : base(a) {}
    template<class A,class B> inline SymmetryQuantizer(const A &a, const B &b) : base(a,b) {}
    template<class A,class B,class C> inline SymmetryQuantizer(const A &a, const B &b, const C &c) : base(a,b,c) {}
    template<class A,class B,class C,class D> inline SymmetryQuantizer(const A &a, const B &b, const C &c, const D &d) : base(a,b,c,d) {}

    inline value_type size() const { return 2*quantizer().size(); }

    inline value_type    operator()(const argument_type &x) const { if (x<0) return quantizer().size()-quantizer()(-x); else return quantizer().size()+quantizer()(x); }
    inline argument_type operator[](      value_type     n) const { value_type s=quantizer().size(); if (n<s) return -quantizer()[s-n]; else return quantizer()[n-s]; }
    inline argument_type quantize  (const argument_type &x) const { if (x<0) return -quantizer().quantize(-x); else return quantizer().quantize(x); }

    template<class G> inline typename quantizerIndexArray   <const Vector<G>, const self>::self operator()(const Vector<G> &X) const { return quantizerIndexArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Vector<G>, const self>::self operator[](const Vector<G> &X) const { return quantizerValueArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Vector<G>, const self>::self quantize  (const Vector<G> &X) const { return quantizerQuantizeArray<const Vector<G>, const self>::self(X, *this); }

    template<class G> inline typename quantizerIndexArray   <const Matrix<G>, const self>::self operator()(const Matrix<G> &X) const { return quantizerIndexArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Matrix<G>, const self>::self operator[](const Matrix<G> &X) const { return quantizerValueArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Matrix<G>, const self>::self quantize  (const Matrix<G> &X) const { return quantizerQuantizeArray<const Matrix<G>, const self>::self(X, *this); }
};

template<class Q, class V=typename Q::argument_type, class S=typename Q::value_type>
class BoundQuantizer : public quantizer_quantizer<Q,V,S>
{
  public:
    typedef BoundQuantizer self;
    typedef quantizer_quantizer<Q,V,S> base;
    
    typedef typename base::quantizer_type quantizer_type;
    typedef typename base::value_type     value_type;
    typedef typename base::argument_type  argument_type;
    
    using base::delta;
    using base::quantizer;

  private:
    argument_type val1;
    argument_type val2;

  public:
    inline BoundQuantizer() : base(), val1(), val2() {}
    inline BoundQuantizer(const argument_type &x, const argument_type &y) : base(), val1(x), val2(y) {}
    inline BoundQuantizer(const self &r) : base(static_cast<const base &>(r)), val1(r.val1), val2(r.val2) {}
    template<class A> inline BoundQuantizer(const argument_type &x, const argument_type &y, const A &a) : base(a), val1(x), val2(y) {}
    template<class A,class B> inline BoundQuantizer(const argument_type &x, const argument_type &y, const A &a, const B &b) : base(a,b), val1(x), val2(y) {}
    template<class A,class B,class C> inline BoundQuantizer(const argument_type &x, const argument_type &y, const A &a, const B &b, const C &c) : base(a,b,c), val1(x), val2(y) {}
    template<class A,class B,class C,class D> inline BoundQuantizer(const argument_type &x, const argument_type &y, const A &a, const B &b, const C &c, const D &d) : base(a,b,c,d), val1(x), val2(y) {}

    inline operator value_type() const { return size(); }
    inline value_type size() const { return (val2-val1+1-delta)/delta+1; }

    inline bool has_value(const argument_type x) const { return (val1<=x) && (x<=val2); } // no ref because binder do not accept it
    inline bool has_index(const value_type  n) const { return (n<size()); }

    inline value_type    operator()(const argument_type &x) const { if (x<val1) return quantizer()(val1); if (x>val2) return quantizer()(val2); return quantizer()(x); }
    inline argument_type operator[](      value_type     n) const { argument_type x=quantizer()[n]; if (x<val1) return val1; if (x>val2) return val2; return x; }
    inline argument_type quantize  (const argument_type &x) const { if (x<val1) return val1; if (x>val2) return val2; return quantizer().quantize(x); }

    template<class G> inline typename quantizerIndexArray   <const Vector<G>, const self>::self operator()(const Vector<G> &X) const { return quantizerIndexArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Vector<G>, const self>::self operator[](const Vector<G> &X) const { return quantizerValueArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Vector<G>, const self>::self quantize  (const Vector<G> &X) const { return quantizerQuantizeArray<const Vector<G>, const self>::self(X, *this); }

    template<class G> inline typename quantizerIndexArray   <const Matrix<G>, const self>::self operator()(const Matrix<G> &X) const { return quantizerIndexArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Matrix<G>, const self>::self operator[](const Matrix<G> &X) const { return quantizerValueArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Matrix<G>, const self>::self quantize  (const Matrix<G> &X) const { return quantizerQuantizeArray<const Matrix<G>, const self>::self(X, *this); }
};

template<class Q, class V=typename Q::argument_type, class S=typename Q::value_type>
class OpenBoundQuantizer : public quantizer_quantizer<Q,V,S>
{
  public:
    typedef OpenBoundQuantizer self;
    typedef quantizer_quantizer<Q,V,S> base;
    
    typedef typename base::quantizer_type quantizer_type;
    typedef typename base::value_type     value_type;
    typedef typename base::argument_type  argument_type;

    using base::delta;
    using base::quantizer;
    
  private:
    argument_type val1;
    argument_type val2;

  public:
    inline OpenBoundQuantizer() : base(), val1(), val2() {}
    inline OpenBoundQuantizer(const argument_type &x, const argument_type &y) : base(), val1(x), val2(y) {}
    inline OpenBoundQuantizer(const self &r) : base(static_cast<const base &>(r)), val1(r.val1), val2(r.val2) {}
    template<class A> inline OpenBoundQuantizer(const argument_type &x, const argument_type &y, const A &a) : base(a), val1(x), val2(y) {}
    template<class A,class B> inline OpenBoundQuantizer(const argument_type &x, const argument_type &y, const A &a, const B &b) : base(a,b), val1(x), val2(y) {}

    inline operator value_type() const { return size(); }
    inline value_type size() const { return (val2-val1)/delta; }

    inline bool has_value(const argument_type x) const { return (val1<=x) && (x<val2); } // no ref parameter because binder do not accept it
    inline bool has_index(const value_type  n) const { return (n<size()); }

    inline value_type  operator()(const argument_type &x) const { if (x<val1) return quantizer()(val1); if (x>=val2) return quantizer()(val2-numeric_limits<argument_type>::epsilon()); return quantizer()(x); }
    inline argument_type operator[](      value_type   n) const { argument_type x=quantizer()[n]; if (x<val1) return val1; if (x>val2) return val2; return x; }
    inline argument_type quantize  (const argument_type &x) const { if (x<val1) return val1; if (x>val2) return val2; return quantizer().quantize(x); }

    template<class G> inline typename quantizerIndexArray   <const Vector<G>, const self>::self operator()(const Vector<G> &X) const { return quantizerIndexArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Vector<G>, const self>::self operator[](const Vector<G> &X) const { return quantizerValueArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Vector<G>, const self>::self quantize  (const Vector<G> &X) const { return quantizerQuantizeArray<const Vector<G>, const self>::self(X, *this); }

    template<class G> inline typename quantizerIndexArray   <const Matrix<G>, const self>::self operator()(const Matrix<G> &X) const { return quantizerIndexArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Matrix<G>, const self>::self operator[](const Matrix<G> &X) const { return quantizerValueArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Matrix<G>, const self>::self quantize  (const Matrix<G> &X) const { return quantizerQuantizeArray<const Matrix<G>, const self>::self(X, *this); }
};

template<class Q, class V=typename Q::argument_type, class I=typename Q::value_type>
class ShiftQuantizer : public quantizer_quantizer<Q,V,I>
{
  public:
    typedef ShiftQuantizer self;
    typedef quantizer_quantizer<Q,V,I> base;
    
    typedef typename base::value_type value_type;    
    typedef typename base::argument_type argument_type;

    using base::quantizer;

  private:
    argument_type val;

  public:
    inline ShiftQuantizer() : base(), val() {}
    inline ShiftQuantizer(const argument_type &x) : base(), val(x) {}
    inline ShiftQuantizer(const self &r) : base(static_cast<const base &>(r)), val(r.val) {}
    template<class A> inline ShiftQuantizer(const argument_type &x, const A &a) : base(a), val(x) {}
    template<class A,class B> inline ShiftQuantizer(const argument_type &x, const A &a, const B &b) : base(a,b), val(x) {}
    template<class A,class B,class C> inline ShiftQuantizer(const argument_type &x, const A &a, const B &b, const C &c) : base(a,b,c), val(x) {}
    template<class A,class B,class C,class D> inline ShiftQuantizer(const argument_type &x, const A &a, const B &b, const C &c, const D &d) : base(a,b,c,d), val(x) {}

    inline value_type     operator()(const argument_type &x) const { return quantizer()(x-val); }
    inline argument_type  operator[](      value_type     n) const { return quantizer()[n]+val; }
    inline argument_type  quantize  (const argument_type &x) const { (*this)[(*this)(x)]; }

    template<class G> inline typename quantizerIndexArray   <const Vector<G>, const self>::self operator()(const Vector<G> &X) const { return quantizerIndexArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Vector<G>, const self>::self operator[](const Vector<G> &X) const { return quantizerValueArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Vector<G>, const self>::self quantize  (const Vector<G> &X) const { return quantizerQuantizeArray<const Vector<G>, const self>::self(X, *this); }

    template<class G> inline typename quantizerIndexArray   <const Matrix<G>, const self>::self operator()(const Matrix<G> &X) const { return quantizerIndexArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Matrix<G>, const self>::self operator[](const Matrix<G> &X) const { return quantizerValueArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Matrix<G>, const self>::self quantize  (const Matrix<G> &X) const { return quantizerQuantizeArray<const Matrix<G>, const self>::self(X, *this); }
};

template<class Q, class V=typename Q::argument_type, class I=typename Q::value_type>
class ScaleQuantizer : public quantizer_quantizer<Q,V,I>
{
  public:
    typedef ScaleQuantizer self;
    typedef quantizer_quantizer<Q,V,I> base;
    
    typedef typename base::value_type value_type;
    typedef typename base::argument_type argument_type;
    
    using base::quantizer;

  private:
    argument_type val;

  public:
    inline ScaleQuantizer() : base(), val() {}
    inline ScaleQuantizer(const argument_type &x) : base(), val(x) {}
    inline ScaleQuantizer(const self &r) : base(static_cast<const base &>(r)), val(r.val) {}
    template<class A> inline ScaleQuantizer(const argument_type &x, const A &a) : base(a), val(x) {}
    template<class A,class B> inline ScaleQuantizer(const argument_type &x, const A &a, const B &b) : base(a,b), val(x) {}
    template<class A,class B,class C> inline ScaleQuantizer(const argument_type &x, const A &a, const B &b, const C &c) : base(a,b,c), val(x) {}
    template<class A,class B,class C,class D> inline ScaleQuantizer(const argument_type &x, const A &a, const B &b, const C &c, const D &d) : base(a,b,c,d), val(x) {}
    template<class A,class B,class C,class D,class E> inline ScaleQuantizer(const argument_type &x, const A &a, const B &b, const C &c, const D &d, const E &e) : base(a,b,c,d,e), val(x) {}

    inline value_type  operator()(const argument_type &x) const { return quantizer()(x/val); }
    inline argument_type operator[](      value_type   n) const { return quantizer()[n]*val; }
    inline argument_type quantize  (const argument_type &x) const { (*this)[(*this)(x)]; }

    template<class G> inline typename quantizerIndexArray   <const Vector<G>, const self>::self operator()(const Vector<G> &X) const { return quantizerIndexArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Vector<G>, const self>::self operator[](const Vector<G> &X) const { return quantizerValueArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Vector<G>, const self>::self quantize  (const Vector<G> &X) const { return quantizerQuantizeArray<const Vector<G>, const self>::self(X, *this); }

    template<class G> inline typename quantizerIndexArray   <const Matrix<G>, const self>::self operator()(const Matrix<G> &X) const { return quantizerIndexArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Matrix<G>, const self>::self operator[](const Matrix<G> &X) const { return quantizerValueArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Matrix<G>, const self>::self quantize  (const Matrix<G> &X) const { return quantizerQuantizeArray<const Matrix<G>, const self>::self(X, *this); }
};


template<class V, class I, class Q>
class ArrayComponentQuantizer
{
  public:
    typedef ArrayComponentQuantizer self;

    typedef Q quantizer_type;
    typedef V argument_type;
    typedef I value_type;

    typedef typename argument_type::template rebind<quantizer_type>::other quantizers_type;
 
  private:
    quantizers_type qs;

  public:
    inline ArrayComponentQuantizer() {}
    inline ArrayComponentQuantizer(int i) : qs(i) {}
    inline ArrayComponentQuantizer(int i, const quantizer_type &q) : qs(i,q) {}
    inline ArrayComponentQuantizer(const quantizers_type &q) : qs(q) {}
    inline ArrayComponentQuantizer(const self &r) : qs(r.qs) {}

    inline quantizers_type       &quantizers()       { return qs; }
    inline const quantizers_type &quantizers() const { return qs; }

    inline value_type size() const { value_type s(quantizers.size()); typename value_type::const_iterator it1=s.begin(); for(typename quantizers_type::const_iterator it2=quantizers.begin();it2!=quantizers.end();++it1,++it2) *it1=(*it2).size(); return s; }

    inline value_type    operator()(const argument_type &v) const { value_type      s(v.size()); typename  value_type::iterator it1=  s.begin(); typename quantizers_type::const_iterator it2=quantizers().begin(); for(typename argument_type::const_iterator it3=v.begin();(it3!=v.end()&&it2!=quantizers().end());++it1,++it2,++it3) (*it1)=(*it2)(*it3);          return s;  }
    inline argument_type operator[](const value_type    &s) const { argument_type   v(s.size()); typename argument_type::iterator it1=  v.begin(); typename quantizers_type::const_iterator it2=quantizers().begin(); for(typename  value_type::const_iterator it3=s.begin();(it3!=s.end()&&it2!=quantizers().end());++it1,++it2,++it3) (*it1)=(*it2)[*it3];          return v;  }
    inline argument_type quantize  (const argument_type &v) const { argument_type res(v.size()); typename argument_type::iterator it1=res.begin(); typename quantizers_type::const_iterator it2=quantizers().begin(); for(typename  value_type::const_iterator it3=v.begin();(it3!=v.end()&&it2!=quantizers().end());++it1,++it2,++it3) (*it1)=(*it2).quantize(*it3); return res;}

    template<class G> inline typename quantizerIndexArray   <const Vector<G>, const self>::self operator()(const Vector<G> &X) const { return typename quantizerIndexArray   <const Vector<G>, const self>::self(X, *this); } 
    template<class G> inline typename quantizerValueArray   <const Vector<G>, const self>::self operator[](const Vector<G> &X) const { return typename quantizerValueArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Vector<G>, const self>::self quantize  (const Vector<G> &X) const { return typename quantizerQuantizeArray<const Vector<G>, const self>::self(X, *this); }

    template<class G> inline typename quantizerIndexArray   <const Matrix<G>, const self>::self operator()(const Matrix<G> &X) const { return typename quantizerIndexArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Matrix<G>, const self>::self operator[](const Matrix<G> &X) const { return typename quantizerValueArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Matrix<G>, const self>::self quantize  (const Matrix<G> &X) const { return typename quantizerQuantizeArray<const Matrix<G>, const self>::self(X, *this); }
};

template<class V, class I>
class ThresholdQuantizer
{
  public:
     typedef ThresholdQuantizer self;

     typedef V argument_type;
     typedef I value_type;

     typedef typename DenseVector<argument_type>::self table_type;

  private:
     table_type t;

  public:
    inline ThresholdQuantizer() {}
    inline ThresholdQuantizer(int n, const char *s) : t(n,s) {}

    inline table_type       &table()       { return t; }
    inline const table_type &table() const { return t; }
 
    inline value_type  operator()(const argument_type &v) const { typename table_type::const_iterator it=find_if(t.begin(),t.end(),bind2nd(greater<argument_type>(),v)); if(it!=t.end()) return distance(t.begin(),it)-1; else if(v==*(--t.end())) return distance(t.begin(),--t.end())-1; }
    inline argument_type operator[](const value_type  &n) const { return (t[n]+t[n+1])/2; }
    inline argument_type quantize  (const argument_type &v) const { return (*this)[(*this)(v)]; }

    template<class G> inline typename quantizerIndexArray   <const Vector<G>, const self>::self operator()(const Vector<G> &X) const { return typename quantizerIndexArray   <const Vector<G>, const self>::self(X, *this); } 
    template<class G> inline typename quantizerValueArray   <const Vector<G>, const self>::self operator[](const Vector<G> &X) const { return typename quantizerValueArray   <const Vector<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Vector<G>, const self>::self quantize  (const Vector<G> &X) const { return typename quantizerQuantizeArray<const Vector<G>, const self>::self(X, *this); }

    template<class G> inline typename quantizerIndexArray   <const Matrix<G>, const self>::self operator()(const Matrix<G> &X) const { return typename quantizerIndexArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerValueArray   <const Matrix<G>, const self>::self operator[](const Matrix<G> &X) const { return typename quantizerValueArray   <const Matrix<G>, const self>::self(X, *this); }
    template<class G> inline typename quantizerQuantizeArray<const Matrix<G>, const self>::self quantize  (const Matrix<G> &X) const { return typename quantizerQuantizeArray<const Matrix<G>, const self>::self(X, *this); }
};


template<class InIt>
typename DenseVector<int>::self histogram(int n, InIt begin, InIt end)
{
  DenseVector<int>::self hist(n, 0);

  for(InIt i=begin; i!=end; ++i)
    ++hist[*i];

  return hist;
}

template<class G>
inline DenseVector<int>::self histogram(int n, const Vector<G> &X)
{
  return histogram(n, X.begin(), X.end());
}

template<class G>
inline DenseVector<int>::self histogram(int n, const Matrix<G> &X)
{
  return histogram(n, X.begin(), X.end());
}



template<class OutIt>
void binarize(OutIt begin, OutIt end, const typename iterator_traits<OutIt>::value_type &x)
{
  typedef typename iterator_traits<OutIt>::value_type value_type;
  for (OutIt it=begin; it!=end; ++it)
  {
    if (*it<x) 
      *it=0;
    else
      *it=numeric_limits<value_type>::max();
  }
}

template<class G>
void binarize(Matrix<G> &X, const typename Matrix<G>::value_type &x)
{
  binarize(X.begin(), X.end(), x);
}


//}

#endif
