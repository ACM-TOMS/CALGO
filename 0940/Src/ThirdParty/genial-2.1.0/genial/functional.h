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

#ifndef FUNCTIONAL_H
#define FUNCTIONAL_H

#include <functional>
#include <complex>

//namespace genial
//{

using namespace std;


template <class Arg1, class Arg2, class Res>
class pointer_to_binary_ref_function : public binary_function<Arg1,Arg2,Res> 
{
  public:
    typedef binary_function<Arg1,Arg2,Res> base;
    typedef typename base::result_type result_type;
    typedef typename base::first_argument_type first_argument_type;
    typedef typename base::second_argument_type second_argument_type;

  public:
    typedef result_type (*function_type)(first_argument_type &, second_argument_type &);
    
  protected:
    function_type func;
    
  public:
    pointer_to_binary_ref_function() {}
    explicit pointer_to_binary_ref_function(const function_type &f)  : func(f) {}
    result_type operator()(first_argument_type &x, second_argument_type &y) const { return func(x,y); }
};

template<class Arg1,class Arg2,class Res>
pointer_to_binary_ref_function<const Arg1,const Arg2,Res> ptr_ref_fun(Res (*f)(const Arg1 &,const Arg2 &))
{
  typedef Res (*function_type)(const Arg1 &, const Arg2 &);
  return pointer_to_binary_ref_function<const Arg1,const Arg2,Res>((function_type)f);
}

template<class V>
struct traits
{
  typedef typename V::value_type       value_type;
  typedef typename V::reference       &reference;
  typedef const typename V::reference &const_reference;
  typedef typename V::pointer         &pointer;
  typedef const typename V::pointer   &const_pointer;
};

template<class V>
struct traits<const V>
{
  typedef const typename V::value_type value_type; // doit rester const sinon erreur dans value_function ...
  typedef const typename V::reference &reference;
  typedef const typename V::reference &const_reference;
  typedef const typename V::pointer   &pointer;
  typedef const typename V::pointer   &const_pointer;
};


#ifdef HIDE_FROM_DOCJET
#else
#define UNARY_FUNCTION_BASE_TYPES                         \
  typedef typename base::result_type result_type;         \
  typedef typename base::argument_type argument_type;     \
  typedef typename base::value_type value_type;           \
  typedef typename base::reference reference;             \
  typedef typename base::const_reference const_reference; \
  typedef typename base::pointer pointer;                 \
  typedef typename base::const_pointer const_pointer;

#define BINARY_FUNCTION_BASE_TYPES                                  \
  typedef typename base::result_type result_type;                   \
  typedef typename base::first_argument_type  first_argument_type;  \
  typedef typename base::second_argument_type second_argument_type; \
  typedef typename base::value_type value_type;                     \
  typedef typename base::reference reference;                       \
  typedef typename base::const_reference const_reference;           \
  typedef typename base::pointer pointer;                           \
  typedef typename base::const_pointer const_pointer;

#define TERTIARY_FUNCTION_BASE_TYPES                                \
  typedef typename base::result_type result_type;                   \
  typedef typename base::first_argument_type  first_argument_type;  \
  typedef typename base::second_argument_type second_argument_type; \
  typedef typename base::third_argument_type  third_argument_type;  \
  typedef typename base::value_type value_type;                     \
  typedef typename base::reference reference;                       \
  typedef typename base::const_reference const_reference;           \
  typedef typename base::pointer pointer;                           \
  typedef typename base::const_pointer const_pointer;

#endif


template<class Op>
struct unary_reference_function_traits
{
  typedef Op function_type;
  typedef typename function_type::argument_type   argument_type;
  typedef typename function_type::result_type     result_type;
  typedef typename function_type::value_type      value_type;
  typedef typename function_type::reference       reference;
  typedef typename function_type::const_reference const_reference;
  typedef typename function_type::pointer         pointer;
  typedef typename function_type::const_pointer   const_pointer;
};

template<class Op>
struct unary_reference_function_traits<const Op>
{
  typedef Op function_type;
  typedef typename function_type::argument_type   argument_type;
  typedef typename function_type::result_type     result_type;
  typedef typename function_type::value_type      value_type;
  typedef typename function_type::const_reference reference;
  typedef typename function_type::const_reference const_reference;
  typedef typename function_type::const_pointer   pointer;
  typedef typename function_type::const_pointer   const_pointer;
};


template<class Op>
struct binary_reference_function_traits
{
  typedef Op function_type;
  typedef typename function_type::first_argument_type  first_argument_type;
  typedef typename function_type::second_argument_type second_argument_type;
  typedef typename function_type::result_type     result_type;
  typedef typename function_type::value_type      value_type;
  typedef typename function_type::reference       reference;
  typedef typename function_type::const_reference const_reference;
  typedef typename function_type::pointer         pointer;
  typedef typename function_type::const_pointer   const_pointer;
};

template<class Op>
struct binary_reference_function_traits<const Op>
{
  typedef Op function_type;
  typedef typename function_type::first_argument_type  first_argument_type;
  typedef typename function_type::second_argument_type second_argument_type;
  typedef typename function_type::result_type     result_type;
  typedef typename function_type::value_type      value_type;
  typedef typename function_type::const_reference reference;
  typedef typename function_type::const_reference const_reference;
  typedef typename function_type::const_pointer   pointer;
  typedef typename function_type::const_pointer   const_pointer;
};


template<class Op>
struct tertiary_reference_function_traits
{
  typedef Op function_type;
  typedef typename function_type::first_argument_type  first_argument_type;
  typedef typename function_type::second_argument_type second_argument_type;
  typedef typename function_type::third_argument_type  third_argument_type;
  typedef typename function_type::result_type     result_type;
  typedef typename function_type::value_type      value_type;
  typedef typename function_type::reference       reference;
  typedef typename function_type::const_reference const_reference;
  typedef typename function_type::pointer         pointer;
  typedef typename function_type::const_pointer   const_pointer;
};

template<class Op>
struct tertiary_reference_function_traits<const Op>
{
  typedef Op function_type;
  typedef typename function_type::first_argument_type  first_argument_type;
  typedef typename function_type::second_argument_type second_argument_type;
  typedef typename function_type::third_argument_type  third_argument_type;
  typedef typename function_type::value_type      value_type;
  typedef typename function_type::const_reference reference;
  typedef typename function_type::const_reference const_reference;
  typedef typename function_type::const_pointer   pointer;
  typedef typename function_type::const_pointer   const_pointer;
};



template<class Arg, class Res>
struct unary_value_function : public unary_function<Arg,Res>
{
  typedef unary_function<Arg,Res> base;
  typedef typename base::result_type result_type;
  typedef typename base::argument_type argument_type;

  typedef result_type value_type;
  typedef value_type reference;
  typedef value_type const_reference; 
  typedef value_type *pointer;
  typedef const value_type *const_pointer;
};

template<class Arg, class Res, class Ref=Res &, class ConstRef=const Res &>
struct unary_reference_function : public unary_function<Arg,Res>
{
  typedef unary_function<Arg,Res> base;
  typedef typename base::result_type result_type;
  typedef typename base::argument_type argument_type;
  
  typedef result_type value_type;
  typedef Ref reference;
  typedef ConstRef const_reference;
  typedef value_type *pointer;
  typedef const value_type *const_pointer;
};


template<class Arg1, class Arg2, class Res>
struct binary_value_function : public binary_function<Arg1,Arg2,Res>
{
  typedef binary_function<Arg1,Arg2,Res> base;
  typedef typename base::result_type result_type;
  typedef typename base::first_argument_type  first_argument_type;
  typedef typename base::second_argument_type second_argument_type;

  typedef result_type value_type;
  typedef value_type reference;
  typedef value_type const_reference;
  typedef value_type *pointer;
  typedef const value_type *const_pointer;
};

template<class Arg1, class Arg2, class Res, class Ref=Res &, class ConstRef=const Res &>
struct binary_reference_function : public binary_function<Arg1,Arg2,Res>
{
  typedef binary_function<Arg1,Arg2,Res> base;
  typedef typename base::result_type result_type;
  typedef typename base::first_argument_type  first_argument_type;
  typedef typename base::second_argument_type second_argument_type;

  typedef result_type value_type;
  typedef Ref reference;
  typedef ConstRef const_reference;
  typedef value_type *pointer;
  typedef const value_type *const_pointer;
};

template<class Arg1,class Arg2,class Arg3,class Res>
struct tertiary_function
{
  typedef Arg1 first_argument;
  typedef Arg2 second_argumet;
  typedef Arg3 third_argument;
  typedef Res  result_type;
};

template<class Arg1, class Arg2, class Arg3, class Res>
struct tertiary_value_function : public tertiary_function<Arg1,Arg2,Arg3,Res>
{
  typedef tertiary_function<Arg1,Arg2,Arg3,Res> base;
  typedef typename base::result_type result_type;
  typedef typename base::first_argument_type  first_argument_type;
  typedef typename base::second_argument_type second_argument_type;
  typedef typename base::third_argument_type  third_argument_type;

  typedef result_type value_type;
  typedef value_type reference;
  typedef const value_type const_reference;
  typedef value_type *pointer;
  typedef const value_type *const_pointer;
};

template<class Arg1, class Arg2, class Arg3, class Res, class Ref=Res &, class ConstRef=const Res &>
struct tertiary_reference_function : public tertiary_function<Arg1,Arg2,Arg3,Res>
{
  typedef tertiary_function<Arg1,Arg2,Arg3,Res> base;
  typedef typename base::result_type result_type;
  typedef typename base::first_argument_type  first_argument_type;
  typedef typename base::second_argument_type second_argument_type;
  typedef typename base::third_argument_type  third_argument_type;

  typedef result_type value_type;
  typedef Ref reference;
  typedef ConstRef const_reference;
  typedef value_type *pointer;
  typedef const value_type *const_pointer;
};




template<class Arg,class Res=Arg>
class id_function : public unary_reference_function<Arg,Res>
{
  public:
    typedef unary_reference_function<Arg,Res> base;
    typedef typename base::argument_type   argument_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;

  public:
    id_function() {}
    
    template<class T> T &      operator()(      T &x)       { return x; }
    template<class T> const T &operator()(const T &x) const { return x; }
};

template<class Op,class Arg1=typename Op::argument_type::value_type,class Arg2=Arg1>
class unary_to_binary_function : public binary_reference_function<Arg1,Arg2,typename Op::result_type,typename Op::reference,typename Op::const_reference>
{
  public:
    typedef binary_reference_function<Arg1,Arg2,typename Op::result_type,typename Op::reference,typename Op::const_reference> base;
    typedef typename base::first_argument_type   first_argument_type;
    typedef typename base::second_argument_type  second_argument_type;
    typedef typename base::reference             reference;
    typedef typename base::const_reference       const_reference;
  
    typedef Op function_type;
    
  private:
    function_type func;
    typedef typename function_type::argument_type argument_type;
    
  public:
    unary_to_binary_function() : func() {}
    unary_to_binary_function(const function_type &f) : func(f) {}
    
    reference       operator()(const first_argument_type &x, const second_argument_type &y)       { return func(argument_type(x,y)); }
    const_reference operator()(const first_argument_type &x, const second_argument_type &y) const { return func(argument_type(x,y)); }
};

template<class Op> inline unary_to_binary_function<const Op> binary(const Op &op) { return unary_to_binary_function<const Op>(op); }


template <class Arg,class Res>
class pointer_to_unary_value_function : public unary_value_function<Arg,Res>
{
  public:
    typedef unary_value_function<Arg,Res> base;
    typedef typename base::argument_type   argument_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;

  protected:
    reference (*ptr)(argument_type);
    
  public:
    pointer_to_unary_value_function() {}
    explicit pointer_to_unary_value_function(reference (*p)(argument_type)) : ptr(p) {}
    
    const_reference operator()(const argument_type &x) const { return ptr(x); }
};

template <class Arg, class Res>
inline pointer_to_unary_value_function<Arg, Res> ptr_val_fun(Res (*op)(Arg))
{
  return pointer_to_unary_value_function<Arg, Res>(op);
}


template<class Op>
class shift_function : public unary_reference_function<typename Op::argument_type,typename Op::value_type,typename Op::reference,typename Op::const_reference>
{
  public:
    typedef unary_reference_function<typename Op::argument_type,typename Op::value_type,typename Op::reference,typename Op::const_reference> base;
    typedef typename base::argument_type   argument_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;

  private:
    Op op;
    argument_type a;
  public:
    shift_function(const Op &x, const argument_type &i) : op(x), a(i) {}
    
    reference       operator()(const argument_type &x)       { return op(x+a); }
    const_reference operator()(const argument_type &x) const { return op(x+a); }
};
template<class Op> inline shift_function<const Op> shift(const Op &op, const typename Op::argument_type &x) { return shift_function<const Op>(op,x); }

template<class Op>
class shift_binary_function : public binary_reference_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type,typename Op::reference,typename Op::const_reference>
{
  public:
    typedef binary_reference_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type,typename Op::reference,typename Op::const_reference> base;
    typedef typename base::first_argument_type   first_argument_type;
    typedef typename base::second_argument_type  second_argument_type;
    typedef typename base::reference             reference;
    typedef typename base::const_reference       const_reference;

  private:
    Op f;
    first_argument_type a;
    second_argument_type b;
    
  public:
    shift_binary_function(const Op &op, const first_argument_type &x) : f(op), a(x), b(x) {}
    shift_binary_function(const Op &op, const first_argument_type &x, const second_argument_type &y) : f(op), a(x), b(y) {}
    
    reference       operator()(const first_argument_type &x, const second_argument_type &y)       { return f(a+x,b+y); }
    const_reference operator()(const first_argument_type &x, const second_argument_type &y) const 
    {
     return f(a+x,b+y); 
    }
};
template<class Op> inline shift_binary_function<Op> shift(const Op &op, const typename Op::first_argument_type &x) { return shift_binary_function<Op>(op,x); }
template<class Op> inline shift_binary_function<Op> shift(const Op &op, const typename Op::first_argument_type &x, const typename Op::second_argument_type &y) { return shift_binary_function<Op>(op,x,y); }


template<class Op>
class scale_function : public unary_reference_function<typename Op::argument_type,typename Op::value_type,typename Op::reference,typename Op::const_reference>
{
  public:
    typedef unary_reference_function<typename Op::argument_type,typename Op::value_type,typename Op::reference,typename Op::const_reference> base;
    typedef typename base::argument_type   argument_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
  
    typedef scale_function self;
    
  private:
    Op op;
    argument_type a;
    
  public:
    scale_function(const Op &x, const argument_type &i) : op(x), a(i) {}
    
    template<class A> scale_function(const self &x, const A &t) : op(x.op,t), a(x.a) {}
    
    reference       operator()(const argument_type &x)       { return op(a*x); }
    const_reference operator()(const argument_type &x) const { return op(a*x); }
};
template<class Op> inline scale_function<Op> scale(const Op &op, const typename Op::argument_type &x) { return scale_function<Op>(op,x); }

template<class Op>
class scale_binary_function : public binary_reference_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type,typename Op::reference,typename Op::const_reference>
{
  public:
    typedef binary_reference_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type,typename Op::reference,typename Op::const_reference> base;
    typedef typename base::first_argument_type   first_argument_type;
    typedef typename base::second_argument_type  second_argument_type;
    typedef typename base::reference             reference;
    typedef typename base::const_reference       const_reference;
  
    typedef scale_binary_function self;
    
  private:
    Op f;
    first_argument_type a;
    second_argument_type b;
    
  public:
    scale_binary_function(const Op &op, const first_argument_type &x) : f(op), a(x), b(x) {}
    scale_binary_function(const Op &op, const first_argument_type &x, const second_argument_type &y) : f(op), a(x), b(y) {}
    
    template<class A> scale_binary_function(const self &x, const A &t) : f(x.f,t), a(x.a), b(x.b) {}
    template<class A,class B> scale_binary_function(const self &x, const A &t1, const B &t2) : f(x.f,t1,t2), a(x.a), b(x.b) {}

    reference       operator()(const first_argument_type &x, const second_argument_type &y)       { return f(a*x,b*y); }
    const_reference operator()(const first_argument_type &x, const second_argument_type &y) const { return f(a*x,b*y); }
};
template<class Op> inline scale_binary_function<Op> scale(const Op &op, const typename Op::first_argument_type &x) { return scale_binary_function<Op>(op,x); }
template<class Op> inline scale_binary_function<Op> scale(const Op &op, const typename Op::first_argument_type &x, const typename Op::second_argument_type &y) { return scale_binary_function<Op>(op,x,y); }

template<class Op>
class rotate_function : public binary_reference_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type,typename Op::reference,typename Op::const_reference>
{
  public:
    typedef binary_reference_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type,typename Op::reference,typename Op::const_reference> base;
    typedef typename base::first_argument_type   first_argument_type;
    typedef typename base::second_argument_type  second_argument_type;
    typedef typename base::reference             reference;
    typedef typename base::const_reference       const_reference;

    typedef rotate_function self;
    
  private:
    Op op;
    double sin_theta;
    double cos_theta;
    
  public:
    rotate_function(const Op &x, double theta) : op(x), sin_theta(sin(theta)), cos_theta(cos(theta)) {}
    template<class A> rotate_function(const self &x, const A &a) : op(x.op,a), sin_theta(x.sin_theta), cos_theta(x.cos_theta) {}
    
    reference       operator()(const first_argument_type &x, const second_argument_type &y)       { return op(x*cos_theta+y*sin_theta, y*cos_theta-x*sin_theta); }
    const_reference operator()(const first_argument_type &x, const second_argument_type &y) const { return op(x*cos_theta+y*sin_theta, y*cos_theta-x*sin_theta); }
};
template<class Op> inline rotate_function<Op> rotate(const Op &op, double theta) { return rotate_function<Op>(op,theta); }

template<class Op>
class cartesian_function : public binary_reference_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type,typename Op::reference,typename Op::const_reference>
{
  public:
    typedef binary_reference_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type,typename Op::reference,typename Op::const_reference> base;
    typedef typename base::first_argument_type   first_argument_type;
    typedef typename base::second_argument_type  second_argument_type;
    typedef typename base::reference             reference;
    typedef typename base::const_reference       const_reference;
  
    typedef cartesian_function self;
    
  private:
    Op op;
    
  public:
    cartesian_function(const Op &x) : op(x) {}
    template<class A> cartesian_function(const self &x, const A &a) : op(x.op,a) {}
    template<class A,class B> cartesian_function(const self &x, const A &a, const B & b) : op(x.op,a,b) {}
    
    reference       operator()(const first_argument_type &x, const second_argument_type &y)       
    {
      complex<first_argument_type> z(x,y); 
      return op(abs(z),arg(z));
    }
    const_reference operator()(const first_argument_type &x, const second_argument_type &y) const 
    {
      complex<first_argument_type> z(x,y); 
      first_argument_type a = abs(z);
      first_argument_type b = arg(z);
      return op(a,b); 
    }
};
template<class Op> inline cartesian_function<Op> cartesian(const Op &op) { return cartesian_function<Op>(op); }


template<class Op>
class square_function : public unary_value_function<typename Op::argument_type,typename Op::value_type>
{
  public:
    typedef unary_value_function<typename Op::argument_type,typename Op::value_type> base;
    typedef typename base::argument_type   argument_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
  
    typedef square_function self;
    
  private:
    Op op;
    
  public:
    square_function(const Op &x) : op(x) {}
    
    template<class A> square_function(const self &x, const A &a) : op(x.op,a) {}
    
    reference       operator()(const argument_type &x)       { return sqr(op(x)); }
    const_reference operator()(const argument_type &x) const { return sqr(op(x)); }
};
template<class Op> inline square_function<Op> sqr1(const Op &op) { return square_function<Op>(op); }

template<class Op>
class square_binary_function : public binary_value_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type>
{
  public:
    typedef binary_value_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type> base;
    typedef typename base::first_argument_type   first_argument_type;
    typedef typename base::second_argument_type  second_argument_type;
    typedef typename base::reference             reference;
    typedef typename base::const_reference       const_reference;
  
    typedef square_binary_function self;
    
  private:
    Op op;
    
  public:
    square_binary_function(const Op &f) : op(f) {}
    
    template<class A> square_binary_function(const self &x, const A &a) : op(x.op,a) {}
    template<class A,class B> square_binary_function(const self &x, const A &a, const B &b) : op(x.op,a,b) {}

    reference       operator()(const first_argument_type &x, const second_argument_type &y)       { return sqr(op(x,y)); }
    const_reference operator()(const first_argument_type &x, const second_argument_type &y) const { return sqr(op(x,y)); }
};
template<class Op> inline square_binary_function<Op> sqr2(const Op &op) { return square_binary_function<Op>(op); }


template<class Op>
class unary_real_function : public unary_value_function<typename Op::argument_type,typename Op::value_type::value_type>
{
  public:
    typedef unary_value_function<typename Op::argument_type,typename Op::value_type::value_type> base;
    typedef typename base::argument_type   argument_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
  
    typedef unary_real_function self;
    
  private:
    Op op;
    
  public:
    unary_real_function(const Op &x) : op(x) {}
        
    reference       operator()(const argument_type &x)       { return real(op(x)); }
    const_reference operator()(const argument_type &x) const { return real(op(x)); }
};
template<class Op> inline unary_real_function<Op> real1(const Op &op) { return unary_real_function<Op>(op); }

template<class Op>
class binary_real_function : public binary_value_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type::value_type>
{
  public:
    typedef binary_value_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type::value_type> base;
    typedef typename base::first_argument_type   first_argument_type;
    typedef typename base::second_argument_type  second_argument_type;
    typedef typename base::reference             reference;
    typedef typename base::const_reference       const_reference;
    
    typedef binary_real_function self;
    
  private:
    Op op;
    
  public:
    binary_real_function(const Op &f) : op(f) {}

    reference       operator()(const first_argument_type &x, const second_argument_type &y)       { return real(op(x,y)); }
    const_reference operator()(const first_argument_type &x, const second_argument_type &y) const { return real(op(x,y)); }
};
template<class Op> inline binary_real_function<Op> real2(const Op &op) { return binary_real_function<Op>(op); }


template<class Op>
class unary_imag_function : public unary_value_function<typename Op::argument_type,typename Op::value_type::value_type>
{
  public:
    typedef unary_value_function<typename Op::argument_type,typename Op::value_type::value_type> base;
    typedef typename base::argument_type   argument_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;

    typedef unary_imag_function self;
    
  private:
    Op op;
    
  public:
    unary_imag_function(const Op &x) : op(x) {}
        
    reference       operator()(const argument_type &x)       { return imag(op(x)); }
    const_reference operator()(const argument_type &x) const { return imag(op(x)); }
};
template<class Op> inline unary_imag_function<Op> imag1(const Op &op) { return unary_imag_function<Op>(op); }

template<class Op>
class binary_imag_function : public binary_value_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type::value_type>
{
  public:
    typedef binary_value_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type::value_type> base;
    typedef typename base::first_argument_type   first_argument_type;
    typedef typename base::second_argument_type  second_argument_type;
    typedef typename base::reference             reference;
    typedef typename base::const_reference       const_reference;
  
    typedef binary_imag_function self;
    
  private:
    Op op;
    
  public:
    binary_imag_function(const Op &f) : op(f) {}

    reference       operator()(const first_argument_type &x, const second_argument_type &y)       { return imag(op(x,y)); }
    const_reference operator()(const first_argument_type &x, const second_argument_type &y) const { return imag(op(x,y)); }
};
template<class Op> inline binary_imag_function<Op> imag2(const Op &op) { return binary_imag_function<Op>(op); }


template<class Op,class V=typename Op::value_type>
class unary_mul_function : public unary_value_function<typename Op::argument_type,typename Op::value_type>
{
  public:
    typedef unary_value_function<typename Op::argument_type,typename Op::value_type> base;
    typedef typename base::argument_type   argument_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;

  private:
    Op op;
    V val;
    
  public:
    unary_mul_function(const Op &f, const V &x) : op(f), val(x) {} 
    
    const_reference operator()(const argument_type &x) const { return val*op(x); } 
};

template<class Op> inline unary_mul_function<Op> mul1(const Op &op, const typename Op::value_type &v) { return unary_mul_function<Op>(op,v); }
template<class Op> inline unary_mul_function<Op,typename Op::value_type::value_type> mul1(const Op &op, const typename Op::value_type::value_type &v) { return unary_mul_function<Op,typename Op::value_type::value_type>(op,v); }

template<class Op,class V=typename Op::value_type>
class binary_mul_function : public binary_value_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type>
{
  public:
    typedef binary_value_function<typename Op::first_argument_type,typename Op::second_argument_type,typename Op::value_type> base;
    typedef typename base::first_argument_type   first_argument_type;
    typedef typename base::second_argument_type  second_argument_type;
    typedef typename base::reference             reference;
    typedef typename base::const_reference       const_reference;

  private:
    Op op;
    V val;
    
  public:
    binary_mul_function(const Op &f, const V &x) : op(f), val(x) {} 
    
    const_reference operator()(const first_argument_type &x, const second_argument_type &y) const { return val*op(x,y); } 
};

template<class Op> inline binary_mul_function<Op> mul2(const Op &op, const typename Op::value_type &v) { return binary_mul_function<Op>(op,v); }
template<class Op> inline binary_mul_function<Op,typename Op::value_type::value_type> mul2(const Op &op, const typename Op::value_type::value_type &v) { return binary_mul_function<Op,typename Op::value_type::value_type>(op,v); }


template <class Arg>
struct bitwise_and : public binary_function<Arg,Arg,Arg>
{
  typedef binary_function<Arg,Arg,Arg> base;
  typedef typename base::first_argument_type   first_argument_type;
  typedef typename base::second_argument_type  second_argument_type;
  typedef typename base::result_type           result_type;
  
  result_type operator()(const first_argument_type &x, const second_argument_type &y) const { return x&y; }
};

template<class Arg>
struct sin_function : public unary_value_function<Arg,Arg>
{
  typedef unary_value_function<Arg,Arg> base;
  typedef typename base::argument_type argument_type;
  typedef typename base::result_type   result_type;

  result_type operator()(const argument_type &x) const { return sin(x); }
};

template<class Arg>
struct cos_function : public unary_value_function<Arg,Arg>
{
  typedef unary_value_function<Arg,Arg> base;
  typedef typename base::argument_type argument_type;
  typedef typename base::result_type   result_type;

  result_type operator()(const argument_type &x) const { return cos(x); }
};

template<class Arg>
struct polar_function : public unary_value_function<Arg,complex<Arg> >
{
  typedef unary_value_function<Arg,complex<Arg> > base;
  typedef typename base::argument_type argument_type;
  typedef typename base::result_type   result_type;

  result_type operator()(const argument_type &x) const { return result_type(cos(x),sin(x)); }
};

template<class Arg,class Res=Arg>
class geometric_serie_function : public unary_value_function<Arg,Res>
{
  public:
    typedef unary_value_function<Arg,Res> base;
    typedef typename base::argument_type   argument_type;
    typedef typename base::result_type     result_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;

    typedef geometric_serie_function self;
    
  private:
    mutable result_type v;  
    result_type r;
    
  public:
    geometric_serie_function(const result_type &ratio) : v(1), r(ratio) {}
    
    result_type operator()(const argument_type &x) const { result_type t=v*x; v*=r; return t; }
};



template <class Op>
class reference_binder1st : public binary_reference_function_traits<Op> 
{
  public:
    typedef reference_binder1st self;
    typedef binary_reference_function_traits<Op>  base;
    
    typedef Op function_type;
    typedef typename base::second_argument_type argument_type;
    typedef typename base::result_type     result_type;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    
  protected:
    function_type op;
    typename function_type::first_argument_type val;
    
  public:
    template<class V> explicit inline reference_binder1st(const V &v) : op(), val(v) {}
    template<class V> inline reference_binder1st(const function_type &f, const V &v) : op(f), val(v) {}

    inline typename function_type::first_argument_type       &value()       { return val; }
    inline const typename function_type::first_argument_type &value() const { return val; }

    inline reference       operator()(const argument_type &x)       { return op(val,x); }
    inline const_reference operator()(const argument_type &x) const { return op(val,x); }
};

template <class Op,class T>
inline reference_binder1st<const Op> reference_bind1st(const Op &op, const T &x)
{
  return reference_binder1st<const Op>(op,typename Op::first_argument_type(x));
}


template <class Op>
class reference_binder2nd : public binary_reference_function_traits<Op> 
{
  public:
    typedef reference_binder2nd self;
    typedef binary_reference_function_traits<Op>  base;
    
    typedef Op function_type;
    typedef typename base::first_argument_type argument_type;
    typedef typename base::result_type     result_type;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    
  protected:
    function_type op;
    typename function_type::second_argument_type val;
    
  public:
    template<class V> explicit inline reference_binder2nd(const V &v) : op(), val(v) {}
    template<class V> inline reference_binder2nd(const function_type &f, const V &v) : op(f), val(v) {}

    inline typename function_type::first_argument_type       &value()       { return val; }
    inline const typename function_type::first_argument_type &value() const { return val; }

    inline reference       operator()(const argument_type &x)       { return op(x,val); }
    inline const_reference operator()(const argument_type &x) const { return op(x,val); }
};

template <class Op,class T>
inline reference_binder2nd<const Op> reference_bind2nd(const Op &op, const T &x)
{
  return reference_binder2nd<const Op>(op,typename Op::first_argument_type(x));
}



//template<class Arg1, class Arg2, class Arg3, class Res>
//struct tertiary_function
//{
//  typedef Arg1 first_argument_type;
//  typedef Arg2 sevond_argument_type;
//  typedef Arg3 third_argument_type;
//  typedef Res result_type;
//};
//
//template<class Arg1, class Arg2, class Arg3, class Res>
//struct tertiary_value_function : public tertiary_function<Arg1, Arg2, Arg3, Res>
//{
//  typedef result_type value_type;
//  typedef value_type reference;
//  typedef const value_type const_reference;
//  typedef value_type *pointer;
//  typedef const value_type *const_pointer;
//};
//
//template<class Arg1, class Arg2, class Arg3, class Res, class Ref=Res &, class ConstRef=const Res &>
//struct tertiary_reference_function : public tertiary_function<Arg1, Arg2, Arg3, Res>
//{
//  typedef result_type value_type;
//  typedef Ref reference;
//  typedef ConstRef const_reference;
//  typedef value_type *pointer;
//  typedef const value_type *const_pointer;
//};
//
//template<class TerOp>
//class binder2nd3rd : public unary_reference_function<typename TerOp::first_argument_type, typename TerOp::value_type, typename TerOp::reference, typename TerOp::const_reference>
//{
//  protected:
//    TerOp func;
//    typename TerOp::second_argumend_type y;
//    typename TerOp::third_argument_type  z;
//   
//  public:
//    binder2nd3rd(const typename TerOp::second_argumend_type &arg2, const typename TerOp::third_argumend_type &arg3) : func(), y(arg2), z(arg3) {}
//    binder2nd3rd(const TerOp &op, const typename TerOp::second_argumend_type &arg2, const typename TerOp::third_argumend_type &arg3) : func(f), y(arg2), z(arg3) {}
//    
//    reference       operator()(const argument_type &x)       { return func(x,y,z); }
//    const_reference operator()(const argument_type &x) const { return func(x,y,z); }    
//};
//
//template<class TerOp, class Arg2, class Arg3>
//inline binder2nd3rd<TerOp> bind2nd3rd(const TerOp &op, const Arg2 &y, const Arg3 &z)
//{
//  return binder2nd3rd<TerOp>(op,y,z);
//}
//
//template<class Arg1, class Arg2, class Arg3, class Res>
//class pointer_to_tertiary_function : public tertiary_function<Arg1,Arg2,Arg3,Res> 
//{
//  protected:
//    result_type (*ptr)(first_argument_type,second_argument_type,third_argument_type);
//    
//  public:
//    pointer_to_tertiary_function() :ptr() {}
//    explicit pointer_to_tertiary_function(result_type (*op)(first_argument_type,second_argument_type,third_argument_type)) : ptr(op) {}
//    result_type operator()(const first_argument_type &x, const second_argument_type &y, const third_argument_type &z) const { return ptr(x,y,z); }
//};
//
//template <class Arg1, class Arg2, class Arg3, class Res>
//inline pointer_to_tertiary_function<Arg1,Arg2,Arg3,Res> ptr_fun(Res (*op)(Arg1,Arg2,Arg3)) 
//{
//  return pointer_to_tertiary_function<Arg1,Arg2,Arg3,Res>(op);
//}


template<class Arg,class Res=typename traits<Arg>::value_type>
struct value_function : public unary_reference_function<Arg,Res>
{
  typedef unary_reference_function<Arg,Res> base;
  typedef typename base::argument_type   argument_type;
  typedef typename base::result_type     result_type;
  typedef typename base::reference       reference;
  typedef typename base::const_reference const_reference;

  value_function() {}
  reference       operator()(argument_type       &x)       { return x.value(); }
  const_reference operator()(const argument_type &x) const { return x.value(); }
};

template<class T> struct name_traits          { typedef       typename T::name_type name_type; };
template<class T> struct name_traits<const T> { typedef const typename T::name_type name_type; };

template<class Arg>
struct name_function : public unary_reference_function<Arg,typename name_traits<Arg>::name_type>
{
  typedef unary_reference_function<Arg,typename name_traits<Arg>::name_type> base;
  typedef typename base::argument_type   argument_type;
  typedef typename base::result_type     result_type;
  typedef typename base::reference       reference;
  typedef typename base::const_reference const_reference;

  name_function() {}
  reference       operator()(argument_type       &x)       { return x.name(); }
  const_reference operator()(const argument_type &x) const { return x.name(); }
};

template<class Arg>
struct to_string_function : public unary_value_function<Arg,string>
{
  typedef unary_value_function<Arg,string> base;
  typedef typename base::argument_type   argument_type;
  typedef typename base::result_type     result_type;
  typedef typename base::reference       reference;
  typedef typename base::const_reference const_reference;

  to_string_function() {}
  const_reference operator()(const argument_type &x) const { return to_string(x); }
};

class string_begin_equal_to : public unary_value_function<string,bool>
{
  private:
    string str;

  public:
    explicit string_begin_equal_to(const string &s) : str(s) {}
    explicit string_begin_equal_to(const char   *s) : str(s) {}
  
    const_reference operator()(const argument_type &x) const { return x.substr(0,str.size())==str; }
};

//}

#endif
