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

#ifndef MATRIXGENERATOR_H
#define MATRIXGENERATOR_H

//namespace genial
//{

//Group = Matrices functions

using namespace std;

template<class G> class Matrix;
template<class A> struct matrix_traits;

template<class V>
struct value_matrix_generator_traits
{
  typedef matrix_array_tag array_category;

  typedef V value_type;           
  typedef typename promotion_traits<value_type>::value_type const_value_type;
  typedef value_type reference;
  typedef const value_type const_reference;
  typedef value_type *pointer;
  typedef const value_type *const_pointer; 

  //typedef size_t int_type;
  //typedef ptrdiff_t difference_type;
  typedef int int_type;
  typedef int difference_type;
  typedef MatrixIndex<int_type> index_type;
  typedef MatrixSize<int_type> size_type;
  
  template<class V2> struct array_rebind { typedef Matrix<dense_matrix_generator<V2> > other; }; 
  
  template<class A> struct iterator_rebind       { typedef typename raster_iterator_rebind<      A>::other other; };
  template<class A> struct const_iterator_rebind { typedef typename raster_iterator_rebind<const A>::other other; };
};

template<class V, class Ref=V &, class ConstRef=const V &>
struct reference_matrix_generator_traits
{
  typedef matrix_array_tag array_category;

  typedef V value_type;           
  typedef typename promotion_traits<value_type>::value_type const_value_type;
  typedef Ref reference;
  typedef ConstRef const_reference;
  typedef value_type *pointer;
  typedef const value_type *const_pointer; 

  //typedef size_t int_type;
  //typedef ptrdiff_t difference_type;
  typedef int int_type;
  typedef int difference_type;
  typedef MatrixIndex<int_type> index_type;
  typedef MatrixSize<int_type> size_type;

  template<class V2> struct array_rebind { typedef Matrix<dense_matrix_generator<V2> > other; }; 

  template<class A> struct iterator_rebind       { typedef typename raster_iterator_rebind<      A>::other other; };
  template<class A> struct const_iterator_rebind { typedef typename raster_iterator_rebind<const A>::other other; };
};


template<class V>
class value_matrix_generator : public value_matrix_generator_traits<V>
{
  public:
    typedef value_matrix_generator self;
    typedef value_matrix_generator_traits<V> base;
    typedef typename base::value_type      value_type;      
    typedef typename base::reference       reference;       
    typedef typename base::const_reference const_reference; 
    typedef typename base::pointer         pointer;         
    typedef typename base::const_pointer   const_pointer;   
    typedef typename base::index_type      index_type;      
    typedef typename base::size_type       size_type;       
    typedef typename base::difference_type difference_type; 
    typedef typename base::int_type        int_type;
        
  public:
    inline value_matrix_generator() {}

    inline index_type lower_bound() const { return index_type(0,0); }
    inline size_type size() const { return size_type(0,0); }
    inline void resize(const size_type &s) { assert(s==size()); }
    
    inline int_type nrows() const { return size().nrows(); }
    inline int_type ncols() const { return size().ncols(); }
};


template<class V, class Ref, class ConstRef>
class reference_matrix_generator : public reference_matrix_generator_traits<V,Ref,ConstRef>
{
  public:
    typedef reference_matrix_generator self;
    typedef reference_matrix_generator_traits<V,Ref,ConstRef> base;
    typedef typename base::value_type      value_type;      
    typedef typename base::reference       reference;       
    typedef typename base::const_reference const_reference; 
    typedef typename base::pointer         pointer;         
    typedef typename base::const_pointer   const_pointer;   
    typedef typename base::index_type      index_type;      
    typedef typename base::size_type       size_type;       
    typedef typename base::difference_type difference_type; 
    typedef typename base::int_type        int_type;
        
  public:
    inline reference_matrix_generator() {}

    inline index_type lower_bound() const { return index_type(0,0); }
    inline size_type size() const { return size_type(0,0); }
    inline void resize(const size_type &s) { assert(s==size()); }
    
    inline int_type nrows() const { return size().nrows(); }
    inline int_type ncols() const { return size().ncols(); }
};


template<class V>
class sized_value_matrix_generator : public value_matrix_generator_traits<V>
{
  public:
    typedef sized_value_matrix_generator self;
    typedef value_matrix_generator_traits<V> base;
    typedef typename base::value_type      value_type;      
    typedef typename base::reference       reference;       
    typedef typename base::const_reference const_reference; 
    typedef typename base::pointer         pointer;         
    typedef typename base::const_pointer   const_pointer;   
    typedef typename base::index_type      index_type;      
    typedef typename base::size_type       size_type;       
    typedef typename base::difference_type difference_type; 
    typedef typename base::int_type        int_type;
        
  private:
    size_type si; 

  public:
    inline sized_value_matrix_generator() : si(0,0) {}
    inline sized_value_matrix_generator(const size_type &s) : si(s) {}
    inline sized_value_matrix_generator(int_type m, int_type n) : si(m,n) {}

    inline index_type lower_bound() const { return index_type(0,0); }
    inline size_type size() const { return si; }
    inline void resize(const size_type &s) { si=s; }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }
    
    inline int_type nrows() const { return size().nrows(); }
    inline int_type ncols() const { return size().ncols(); }
    
    inline int_type nelms() const { return size().nelms(); }        
};


template<class V,class Ref=V&,class ConstRef=const V &>
class sized_reference_matrix_generator : public reference_matrix_generator_traits<V,Ref,ConstRef>
{
  public:
    typedef sized_reference_matrix_generator self;
    typedef reference_matrix_generator_traits<V,Ref,ConstRef> base;
    typedef typename base::value_type      value_type;      
    typedef typename base::reference       reference;       
    typedef typename base::const_reference const_reference; 
    typedef typename base::pointer         pointer;         
    typedef typename base::const_pointer   const_pointer;   
    typedef typename base::index_type      index_type;      
    typedef typename base::size_type       size_type;       
    typedef typename base::difference_type difference_type; 
    typedef typename base::int_type        int_type;
        
  protected:
    size_type si; 

  public:
    inline sized_reference_matrix_generator() : si(0,0) {}
    inline sized_reference_matrix_generator(const size_type &s) : si(s) {}
    inline sized_reference_matrix_generator(int_type m, int_type n) : si(m,n) {}

    inline index_type lower_bound() const { return index_type(0,0); }
    inline size_type size() const { return si; }
    inline void resize(const size_type &s) { si=s; }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }
    
    inline int_type nrows() const { return size().nrows(); }
    inline int_type ncols() const { return size().ncols(); }
        
    inline int_type nelms() const { return size().nelms(); }    
};


template<class V>
class square_sized_value_matrix_generator : public value_matrix_generator_traits<V>
{
  public:
    typedef square_sized_value_matrix_generator self;
    typedef value_matrix_generator_traits<V> base;
    typedef typename base::value_type      value_type;      
    typedef typename base::reference       reference;       
    typedef typename base::const_reference const_reference; 
    typedef typename base::pointer         pointer;         
    typedef typename base::const_pointer   const_pointer;   
    typedef typename base::index_type      index_type;      
    typedef typename base::size_type       size_type;       
    typedef typename base::difference_type difference_type; 
    typedef typename base::int_type        int_type;
        
  protected:
    int_type N; 

  public:
    inline square_sized_value_matrix_generator(int_type n) : N(n) {}

    inline index_type lower_bound() const { return index_type(0,0); }
    inline size_type size() const { return size_type(N,N); }
    inline void resize(const size_type &s) { assert(s.nrows()==s.ncols()); N=s.nrows(); }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }
    
    inline int_type nrows() const { return size().nrows(); }
    inline int_type ncols() const { return size().ncols(); } 
            
    inline int_type nelms() const { return size().nelms(); }       
};


template<class V>
class data_matrix_generator : public sized_reference_matrix_generator<V>
{
  public:
    typedef data_matrix_generator self;
    typedef sized_reference_matrix_generator<V> base;
    typedef typename base::value_type      value_type;      
    typedef typename base::reference       reference;       
    typedef typename base::const_reference const_reference; 
    typedef typename base::pointer         pointer;         
    typedef typename base::const_pointer   const_pointer;   
    typedef typename base::index_type      index_type;      
    typedef typename base::size_type       size_type;       
    typedef typename base::difference_type difference_type; 
    typedef typename base::int_type        int_type;
        
  public:
    value_type *data;
      
  public:
    inline data_matrix_generator() : base(0,0), data(NULL) { }
    inline data_matrix_generator(const size_type &d, value_type *p) : base(d), data(p) { }
    inline data_matrix_generator(int_type m,int_type n, value_type *p) : base(m,n), data(p) { }
    template<class G> inline data_matrix_generator(const Matrix<G> &X) : base(X.nrows(),(X.ncols()*sizeof(typename Matrix<G>::value_type))/sizeof(value_type)), data((pointer)&*X.begin()) {}
    
    using base::size;
    using base::nrows;
    using base::ncols;
    using base::lower_bound;
    using base::nelms;
    
    inline int_type row_stride() const { return ncols(); }

    inline reference       operator[](const index_type &p)       { assert(size().withinbounds(p)); return data[p.i*size().ncols()+p.j]; }
    inline const_reference operator[](const index_type &p) const { assert(size().withinbounds(p)); return data[p.i*size().ncols()+p.j]; }

    inline reference       operator()(int_type i, int_type j)       { return (*this)[index_type(i,j)]; }
    inline const_reference operator()(int_type i, int_type j) const { return (*this)[index_type(i,j)]; }
    
    inline void set_data(value_type *p) { data=p; }   
   
    inline pointer       begin()       { return data; }
    inline const_pointer begin() const { return data; }
    inline pointer       end  ()       { return data+nelms(); }
    inline const_pointer end  () const { return data+nelms(); } 
       
    inline void swap(self &x) { base::swap(x); std::swap(data,x.data); }    
};

template<class V> inline void swap(data_matrix_generator<V> &x, data_matrix_generator<V> &y) { x.swap(y); }

template<class V>
struct DataMatrix : public reference_matrix_generator_traits<V>    
{
  typedef V value_type;
  typedef data_matrix_generator<value_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self;
}; 

template<class V> typename DataMatrix<      V>::self data_matrix(int m, int n,       V *p) { return typename DataMatrix<      V>::self(m,n,p); }
template<class V> typename DataMatrix<const V>::self data_matrix(int m, int n, const V *p) { return typename DataMatrix<const V>::self(m,n,p); }



template<class V,class A>
class dense_matrix_generator : public data_matrix_generator<V>
{
  public:
    typedef dense_matrix_generator self;
    typedef data_matrix_generator<V> base;
    typedef typename base::value_type      value_type;      
    typedef typename base::reference       reference;       
    typedef typename base::const_reference const_reference; 
    typedef typename base::pointer         pointer;         
    typedef typename base::const_pointer   const_pointer;   
    typedef typename base::index_type      index_type;      
    typedef typename base::size_type       size_type;       
    typedef typename base::difference_type difference_type; 
    typedef typename base::int_type        int_type;

  protected:
    A alloc;

  public:
    using base::data;
    using base::si;
    
    using base::size;
    using base::nrows;
    using base::ncols;
    using base::lower_bound;
    using base::nelms;
    using base::begin;
    using base::end;
    
    inline dense_matrix_generator() : base() {}
    inline explicit dense_matrix_generator(const size_type &d)                 : base(d  ,alloc.allocate(d.nelms())) { uninitialized_fill(begin(),end()); }
    inline dense_matrix_generator(const size_type &d,     const value_type &v) : base(d  ,alloc.allocate(d.nelms())) { uninitialized_fill(begin(),end(),v           ); }
    dense_matrix_generator(const size_type &d,     const char       *s) : base(d  ,alloc.allocate(d.nelms())) { istringstream iss(s        ); for (pointer i=begin(); i!=end(); ++i) iss>>*i; }
    dense_matrix_generator(const size_type &d,     const string     &s) : base(d  ,alloc.allocate(d.nelms())) { istringstream iss(s.c_str()); for (pointer i=begin(); i!=end(); ++i) iss>>*i; }

    inline dense_matrix_generator(int_type m, int_type n                     ) : base(m,n,alloc.allocate(m*n      )) { uninitialized_fill(begin(),end()); }
    inline dense_matrix_generator(int_type m, int_type n, const value_type &v) : base(m,n,alloc.allocate(m*n      )) { uninitialized_fill(begin(),end(),v           ); }
    dense_matrix_generator(int_type m, int_type n, const char       *s) : base(m,n,alloc.allocate(m*n      )) { istringstream iss(s        ); for (pointer i=begin(); i!=end(); ++i) iss>>*i; }
    dense_matrix_generator(int_type m, int_type n, const string     &s) : base(m,n,alloc.allocate(m*n      )) { istringstream iss(s.c_str()); for (pointer i=begin(); i!=end(); ++i) iss>>*i; }

    inline dense_matrix_generator                           (const self       &x) : base(x.size(),alloc.allocate(x.nelms())) { copy(x.begin(),x.end(),raw_storage_iterator<pointer,value_type>(begin())); }
    template<class G> inline explicit dense_matrix_generator(const Matrix<G>  &x) : base(x.size(),alloc.allocate(x.nelms())) { copy(x.begin(),x.end(),raw_storage_iterator<pointer,value_type>(begin())); }
    template<class G> inline explicit dense_matrix_generator(const Array<2,G> &x) : base(x.size(),alloc.allocate(x.nelms())) { copy(x.begin(),x.end(),raw_storage_iterator<pointer,value_type>(begin())); }

    inline ~dense_matrix_generator() { alloc.deallocate(data,nelms()); }

    inline int_type stride() const { return ncols(); }

    inline void swap(self &x) { base::swap(x); std::swap(alloc,x.alloc); }
       
    inline reference       operator[](const index_type &p)       { assert(size().withinbounds(p)); return data[p.i*size().ncols()+p.j]; }
    inline const_reference operator[](const index_type &p) const { assert(size().withinbounds(p)); return data[p.i*size().ncols()+p.j]; }

    inline reference       operator()(int_type i, int_type j)       { return (*this)[index_type(i,j)]; }
    inline const_reference operator()(int_type i, int_type j) const { return (*this)[index_type(i,j)]; }
    
    inline void inv(const value_type &v, const index_type &p) { assert(size().withinbounds(p)); (*this)[p]=v; }
    
    inline void resize(const size_type &s) { if (s.nelms()==nelms()) { si=s; return; } alloc.deallocate(data,nelms()); si=s; data=alloc.allocate(nelms()); uninitialized_fill(begin(),end()); }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }

    inline bool operator==(const self &x) const { return size()==x.size() &&  equal(begin(), end(), x.begin()); }
    inline bool operator!=(const self &x) const { return size()!=x.size() || !equal(begin(), end(), x.begin()); }
};

template<class V,class A> inline void swap(dense_matrix_generator<V,A> &x, dense_matrix_generator<V,A> &y) { x.swap(y); }

//{unsecret}
//{group:Matrices Interfaces}
//Summary: Dense matrix type
//Arguments:
//  V - Type of the elements
//  A - Optional allocator
//Example:
//  DenseMatrix<int>::self X(3,3,"1 2 3 4 5 6 7 8 9");
//  typedef DenseMatrix<complex<float> >::self MyMatrix;
template<class V,class A> 
struct DenseMatrix : public reference_matrix_generator_traits<V>  
{
  typedef V value_type;
  typedef A allocator_type;
  typedef dense_matrix_generator<value_type,allocator_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self;
}; 


template<int M,int N,class V>
class tiny_matrix_generator : public reference_matrix_generator_traits<V>
{
  public:
    typedef tiny_matrix_generator self;
    typedef reference_matrix_generator_traits<V> base;
    typedef typename base::value_type      value_type;      
    typedef typename base::reference       reference;       
    typedef typename base::const_reference const_reference; 
    typedef typename base::pointer         pointer;         
    typedef typename base::const_pointer   const_pointer;   
    typedef typename base::index_type      index_type;      
    typedef typename base::size_type       size_type;       
    typedef typename base::difference_type difference_type; 
    typedef typename base::int_type        int_type;
        
    template<class V2> struct array_rebind { typedef Matrix<tiny_matrix_generator<M,N,V2> > other; };

  protected:
    __aligned value_type data[M*N];

  public:
    inline tiny_matrix_generator()                             { }
    inline explicit tiny_matrix_generator(const value_type &v) { uninitialized_fill(begin(),end(),v           ); }
    inline explicit tiny_matrix_generator(const char   *s) { istringstream iss(s        ); for (pointer i=begin();i!=end();++i) iss>>*i; }
    inline explicit tiny_matrix_generator(const string &s) { istringstream iss(s.c_str()); for (pointer i=begin();i!=end();++i) iss>>*i; }

    inline tiny_matrix_generator                           (const self       &x) { assert(size()==x.size()); copy(x.begin(),x.end(),raw_storage_iterator<pointer,value_type>(begin())); }
    template<class G> inline explicit tiny_matrix_generator(const Matrix<G>  &x) { assert(size()==x.size()); copy(x.begin(),x.end(),raw_storage_iterator<pointer,value_type>(begin())); }
    template<class G> inline explicit tiny_matrix_generator(const Array<2,G> &x) { assert(size()==x.size()); copy(x.begin(),x.end(),raw_storage_iterator<pointer,value_type>(begin())); }

    inline reference       operator[](const index_type &p)       { assert(size().withinbounds(p)); return data[p.i*size().ncols()+p.j]; }
    inline const_reference operator[](const index_type &p) const { assert(size().withinbounds(p)); return data[p.i*size().ncols()+p.j]; }
    
    inline void inv(const value_type &y, const int_type &p) { assert(size().withinbounds(p)); (*this)[p]=y; }

    inline size_type  size       () const { return size_type (M,N); }
    inline index_type lower_bound() const { return index_type(0,0); }
    inline int_type   nelms      () const { return M*N; }
    inline int_type   stride     () const { return N; }

    inline pointer       begin()       { return data; }
    inline const_pointer begin() const { return data; }
    inline pointer       end  ()       { return data+nelms(); }
    inline const_pointer end  () const { return data+nelms(); }

    inline void resize(const size_type &s) { assert(size()==s); }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }
  
    inline bool operator==(const self &r) const { return size()==r.size() &&  equal(begin(), end(), r.begin()); }
    inline bool operator!=(const self &r) const { return !(*this==r); }
};

//{unsecret}
//{group:Matrices Interfaces}
//Summary: Tiny dense matrix type
//Arguments:
//  M - Number of rows
//  N - Number of columns
//  V - Type of the elements
//Example:
//  TinyMatrix<3,3,int>::self X("1 2 3 4 5 6 7 8 9");
//  typedef TinyDenseMatrix<3,3,complex<float> >::self MyMatrix;
template<int M,int N,class V>
struct TinyMatrix : public reference_matrix_generator_traits<V>   
{ 
  typedef V value_type;
  typedef tiny_matrix_generator<M,N,value_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class V,class A=alignment_allocator<V> >
class lower_triangle_dense_matrix_generator : public reference_matrix_generator_traits<V,V&,const V>
{
  public:
    typedef lower_triangle_dense_matrix_generator self;
    typedef reference_matrix_generator_traits<V,V&,const V> base;
    typedef typename base::value_type      value_type;      
    typedef typename base::reference       reference;       
    typedef typename base::const_reference const_reference; 
    typedef typename base::pointer         pointer;         
    typedef typename base::const_pointer   const_pointer;   
    typedef typename base::index_type      index_type;      
    typedef typename base::size_type       size_type;       
    typedef typename base::difference_type difference_type; 
    typedef typename base::int_type        int_type;
        
    template<class T> struct iterator_rebind       { typedef typename lower_triangle_iterator_rebind<      T>::other other; };
    template<class T> struct const_iterator_rebind { typedef typename lower_triangle_iterator_rebind<const T>::other other; };        

  protected:
    int_type N;
    value_type *data;
    A alloc;

  public:
    inline lower_triangle_dense_matrix_generator() : N(0), data(0) {}

    inline explicit lower_triangle_dense_matrix_generator(int_type n            ) : N(n), data() { data=alloc.allocate(nelms()); uninitialized_fill(begin(),end()); }
    inline lower_triangle_dense_matrix_generator(int_type n, const value_type &v) : N(n), data() { data=alloc.allocate(nelms()); uninitialized_fill(begin(),end(),v           ); }
    inline lower_triangle_dense_matrix_generator(int_type n, const char       *s) : N(n), data() { data=alloc.allocate(nelms()); istringstream iss(s        ); for (pointer i=begin(); i!=end(); ++i) iss>>*i; }
    inline lower_triangle_dense_matrix_generator(int_type n, const string     &s) : N(n), data() { data=alloc.allocate(nelms()); istringstream iss(s.c_str()); for (pointer i=begin(); i!=end(); ++i) iss>>*i; }

    inline lower_triangle_dense_matrix_generator                           (const self       &x) : N(x.nrows()), data() { assert(x.nrows()==x.ncols()); data=alloc.allocate(nelms()); copy(x.begin(),x.end(),raw_storage_iterator<pointer,value_type>(begin())); }
    template<class G> inline explicit lower_triangle_dense_matrix_generator(const Matrix<G>  &x) : N(x.nrows()), data() { assert(x.nrows()==x.ncols()); data=alloc.allocate(nelms()); copy(x.begin(),x.end(),raw_storage_iterator<pointer,value_type>(begin())); }
    template<class G> inline explicit lower_triangle_dense_matrix_generator(const Array<2,G> &x) : N(x.nrows()), data() { assert(x.nrows()==x.ncols()); data=alloc.allocate(nelms()); copy(x.begin(),x.end(),raw_storage_iterator<pointer,value_type>(begin())); }

    inline ~lower_triangle_dense_matrix_generator() { alloc.deallocate(data,nelms()); }
       
    inline reference       operator[](const index_type &p)       { assert(size().withinbounds(p)); assert(p.i>=p.j); return data[(p.i*(p.i+1))/2+p.j]; }
    inline const_reference operator[](const index_type &p) const { assert(size().withinbounds(p)); if    (p.i>=p.j)  return data[(p.i*(p.i+1))/2+p.j]; else return 0; }
    
    inline void inv(const value_type &v, const index_type &p) { (*this)[p]=v; }

    inline index_type lower_bound() const { return index_type(0,0); }
    inline size_type size() const { return size_type(N,N); }
    inline int_type nelms() const { return (N*(N+1))/2; }

    inline int_type nrows() const { return size().nrows(); }
    inline int_type ncols() const { return size().ncols(); }

    inline pointer       begin()       { return data; }
    inline const_pointer begin() const { return data; }
    inline pointer       end  ()       { return data+nelms(); }
    inline const_pointer end  () const { return data+nelms(); }

    inline void resize(int_type n) { if (N==n) return; alloc.deallocate(data,nelms()); N=n; data=alloc.allocate(nelms()); uninitialized_fill(begin(),end()); }
    inline void resize(const size_type &s) { assert(s.nrows()==s.ncols()); resize(s.nrows()); }

    inline bool operator==(const self &x) const { return nelms()==x.nelms() &&  equal(begin(),end(),x.begin()); }
    inline bool operator!=(const self &x) const { return nelms()!=x.nelms() || !equal(begin(),end(),x.begin()); }
};

template<class V,class A=alignment_allocator<V> > 
struct LowerTriangleDenseMatrix : reference_matrix_generator_traits<V> 
{
  typedef V value_type;
  typedef A alocator_type;
  typedef lower_triangle_dense_matrix_generator<value_type,alocator_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self;
};


template<class A1,class A2=A1,class V=PROMOTE2(typename A1::value_type,typename A2::value_type),int C=0>
class outer_product_matrix_generator : public two_arrays_value_generator<A1,A2,V,C>
{
  public:
    typedef outer_product_matrix_generator self;
    typedef two_arrays_value_generator<A1,A2,V,C> base;
    typedef matrix_array_tag array_category;
    typedef typename base::first_array_type  first_array_type;
    typedef typename base::second_array_type second_array_type; 
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;   
    typedef typename base::int_type int_type;
    typedef MatrixSize<int_type> size_type;
    typedef MatrixIndex<int_type> index_type;
 
  public:
    inline outer_product_matrix_generator() : base() {}
    inline outer_product_matrix_generator(first_array_type &x, second_array_type &y) : base(x,y) { }
    inline outer_product_matrix_generator(const self &x) : base(x) {}
    template<class T1,class T2,class V2,int C2> inline outer_product_matrix_generator(const outer_product_matrix_generator<T1,T2,V2,C2> &x) : base(x) {}

    inline first_array_type        &first_array ()       { return base::first_array (); }
    inline const first_array_type  &first_array () const { return base::first_array (); }
    inline second_array_type       &second_array()       { return base::second_array(); }
    inline const second_array_type &second_array() const { return base::second_array(); }

    inline const_reference operator[](const index_type &z) const { return first_array()[z.i]*second_array()[z.j]; }
    
    inline index_type lower_bound() const { return index_type(first_array().lower_bound(), second_array().lower_bound()); }
    inline size_type size() const { return size_type(first_array().size(),second_array().size()); }  
    inline void resize(const size_type &s) { first_array().resize(s.nrows()); second_array().resize(s.ncols()); }
};

template<class A1,class A2=A1,class V=PROMOTE2(typename A1::value_type,typename A2::value_type),int C=0>
struct OuterProductMatrix : public two_arrays_value_traits<A1,A2,V>
{
  typedef A1 first_array_type;
  typedef A2 second_array_type;
  typedef V value_type;
  typedef outer_product_matrix_generator<first_array_type,second_array_type,value_type,C> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
}; 

template<class G1,class G2>
typename OuterProductMatrix<const Vector<G1>,const Vector<G2> >::self outer_product(const Vector<G1> &X, const Vector<G2> &Y)
{
  return typename OuterProductMatrix<const Vector<G1>,const Vector<G2> >::self(X,Y);
}


template<class A>
class row_split_generator : public array_value_generator<A,typename MatrixRow<A>::self>
{
  public:
    typedef row_split_generator self;
    typedef array_value_generator<A,typename MatrixRow<A>::self> base;
    typedef vector_array_tag array_category;
    typedef typename base::array_type      array_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;       
    typedef typename base::int_type int_type;
    typedef typename base::int_type size_type;
    typedef typename base::int_type index_type;
    
    template<class V2> struct array_rebind { typedef typename matrixrow_rebind<typename array_type::template rebind<V2>::other>::other other; };

  public:
    inline row_split_generator(array_type &x) : base(x) {}

    using base::array;
    
    inline reference       operator[](const index_type &n)       { return row(array(),n); }
    inline const_reference operator[](const index_type &n) const { return row(array(),n); }

    inline index_type lower_bound() const { return array().row_lower_bound(); }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }
    inline size_type size () const { return array().nrows(); }
    inline int_type  nelms() const { return array().nrows(); }  
    inline void resize(size_type n) { assert(size()==n); }
};

template<class A>
struct RowSplit
{
  typedef A array_type;
  typedef row_split_generator<array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

//{unsecret}
//{noAutoLink}
//Summary : Isolates each row
//Parameters:
//  X : The matrix
//Return: A vector representing each row: a vector of vectors.
//Example:
//  DenseMatrix<int>::self X(2,2,"1 2 3 4");
//  cout << rows(X)[0] << endl; // [1 2]
//See: ^row^, ^col^, ^cols^
template<class G> inline typename RowSplit<      Matrix<G> >::self rows(      Matrix<G> &X) { return typename RowSplit<      Matrix<G> >::self(X); }
//{unsecret}
template<class G> inline typename RowSplit<const Matrix<G> >::self rows(const Matrix<G> &X) { return typename RowSplit<const Matrix<G> >::self(X); }


template<class A>
class col_split_generator : public array_value_generator<A,typename MatrixCol<A>::self>
{
  public:
    typedef col_split_generator self;
    typedef array_value_generator<A,typename MatrixCol<A>::self> base;
    typedef vector_array_tag array_category;
    typedef typename base::array_type array_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;           
    typedef typename base::int_type int_type;
    typedef typename base::int_type size_type;
    typedef typename base::int_type index_type;
    
    template<class V2> struct array_rebind { typedef typename matrixrow_rebind<typename array_type::template rebind<V2>::other>::other other; };

  public:
    inline col_split_generator(array_type &x) : base(x) {}

    inline array_type       &array()       { return base::array(); }
    inline const array_type &array() const { return base::array(); }    

    inline reference       operator[](const index_type &n)       { return col(array(),n); }
    inline const_reference operator[](const index_type &n) const { return col(array(),n); }

    inline index_type lower_bound() const { return array().col_lower_bound(); }
    inline void set_lower_bound(const index_type &i) { assert(i==lower_bound()); }
    inline size_type size () const { return array().ncols(); }
    inline int_type  nelms() const { return array().ncols(); }  
    inline void resize(size_type n) { assert(size()==n); }
};

template<class A>
struct ColSplit
{
  typedef A array_type;
  typedef col_split_generator<array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

//{unsecret}
//Summary : Isolates each column
//Parameters:
//  X : The matrix
//Return: A vector representing each column: a vector of vectors
//Example:
//  DenseMatrix<int>::self X(2,2,"1 2 3 4");
//  cout << cols(X)[0] << endl; // [1 3]
//See: ^cols^, ^row^, ^rows^
template<class G> inline typename ColSplit<      Matrix<G> >::self cols(      Matrix<G> &X) { return typename ColSplit<      Matrix<G> >::self(X); }
//{unsecret}
template<class G> inline typename ColSplit<const Matrix<G> >::self cols(const Matrix<G> &X) { return typename ColSplit<const Matrix<G> >::self(X); }


template<class T>
class col_join_matrix_generator : public array_reference_generator<T,typename T::value_type::value_type,typename T::value_type::reference,typename T::value_type::const_reference>
{
  public:
    typedef col_join_matrix_generator self;
    typedef array_reference_generator<T,typename T::value_type::value_type,typename T::value_type::reference,typename T::value_type::const_reference> base;
    typedef matrix_array_tag array_category;
    typedef typename base::array_type array_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;           
    typedef typename base::int_type int_type;
    typedef MatrixSize <int_type> size_type;
    typedef MatrixIndex<int_type> index_type;
    
    //todo: array_rebind

  public:
    inline col_join_matrix_generator(array_type &x) : base(x) {}

    using base::array;

    inline reference       operator[](const index_type &n)       { return array()[n.j][n.i]; }
    inline const_reference operator[](const index_type &n) const { return array()[n.j][n.i]; }

    inline index_type lower_bound() const { return index_type(row_lower_bound(), col_lower_bound()); }
    inline int_type row_lower_bound() const { return array().front().lower_bound(); }
    inline int_type col_lower_bound() const { return array().lower_bound(); }

    inline size_type size () const { return size_type(nrows(),ncols()); }
    inline int_type nrows() const { return array().front().size(); }
    inline int_type ncols() const { return array().size(); }

    inline void resize(size_type n) { assert(size()==n); }
};

template<class A>
struct ColJoinMatrix
{
  typedef A array_type;
  typedef col_join_matrix_generator<array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G> typename ColJoinMatrix<      Vector<G> >::self col_join(      Vector<G> &X) { return typename ColJoinMatrix<      Vector<G> >::self(X); }
template<class G> typename ColJoinMatrix<const Vector<G> >::self col_join(const Vector<G> &X) { return typename ColJoinMatrix<const Vector<G> >::self(X); }



template<class T>
class row_join_matrix_generator : public array_reference_generator<T,typename T::value_type::value_type,typename T::value_type::reference,typename T::value_type::const_reference>
{
  public:
    typedef row_join_matrix_generator self;
    typedef array_reference_generator<T,typename T::value_type::value_type,typename T::value_type::reference,typename T::value_type::const_reference> base;
    typedef matrix_array_tag array_category;
    typedef typename base::array_type array_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;           
    typedef typename base::int_type int_type;
    typedef MatrixSize <int_type> size_type;
    typedef MatrixIndex<int_type> index_type;
    
    //todo: array_rebind

  public:
    inline row_join_matrix_generator(array_type &x) : base(x) {}
    
    using base::array;
    
    inline reference       operator[](const index_type &n)       { return array()[n.i][n.j]; }
    inline const_reference operator[](const index_type &n) const { return array()[n.i][n.j]; }

    inline index_type lower_bound() const { return index_type(row_lower_bound(), col_lower_bound()); }
    inline int_type row_lower_bound() const { return array().lower_bound(); }
    inline int_type col_lower_bound() const { return array().front().lower_bound(); }

    inline size_type size () const { return size_type(nrows(),ncols()); }
    inline int_type nrows() const { return array().size(); }
    inline int_type ncols() const { return array().front().size(); }

    inline void resize(size_type n) { assert(size()==n); }
};

template<class A>
struct RowJoinMatrix
{
  typedef A array_type;
  typedef row_join_matrix_generator<array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G> typename RowJoinMatrix<      Vector<G> >::self row_join(      Vector<G> &X) { return typename RowJoinMatrix<      Vector<G> >::self(X); }
template<class G> typename RowJoinMatrix<const Vector<G> >::self row_join(const Vector<G> &X) { return typename RowJoinMatrix<const Vector<G> >::self(X); }


template<class T>
class row_flip_matrix_generator : public array_generator<T>
{
  public:
    typedef row_flip_matrix_generator self;
    typedef array_generator<T> base;
    ARRAY_BASE_TYPES

  public:
    inline row_flip_matrix_generator(array_type &x) : base(x) {}
    
    inline array_type       &array()       { return base::array(); }
    inline const array_type &array() const { return base::array(); }        
     
    inline reference       operator[](const index_type &p)       { return array()(array().row_upper_bound()+array().row_lower_bound()-p.i,p.j); }
    inline const_reference operator[](const index_type &p) const { return array()(array().row_upper_bound()+array().row_lower_bound()-p.i,p.j); }
};

template<class T>
struct rowFlipArray 
{
  typedef T array_type;
  typedef row_flip_matrix_generator<array_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self; 
}; 

template<class G> inline typename rowFlipArray<      Matrix<G> >::self row_flip(      Matrix<G> &X) { return typename rowFlipArray<      Matrix<G> >::self(X); }
template<class G> inline typename rowFlipArray<const Matrix<G> >::self row_flip(const Matrix<G> &X) { return typename rowFlipArray<const Matrix<G> >::self(X); }


template<class T>
class col_flip_matrix_generator : public array_generator<T>
{
  public:
    typedef col_flip_matrix_generator self;
    typedef array_generator<T> base;
    ARRAY_BASE_TYPES

  public:
    inline col_flip_matrix_generator(array_type &x) : base(x) {}

    inline array_type       &array()       { return base::array(); }
    inline const array_type &array() const { return base::array(); }        
     
    inline reference       operator[](const index_type &p)       { return array()(p.i,array().col_upper_bound()+array().col_lower_bound()-p.j); }
    inline const_reference operator[](const index_type &p) const { return array()(p.i,array().col_upper_bound()+array().col_lower_bound()-p.j); }
};

template<class T>
struct colFlipArray 
{
  typedef T array_type;
  typedef col_flip_matrix_generator<array_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self; 
}; 

template<class G> inline typename colFlipArray<      Matrix<G> >::self col_flip(      Matrix<G> &X) { return typename colFlipArray<      Matrix<G> >::self(X); }
template<class G> inline typename colFlipArray<const Matrix<G> >::self col_flip(const Matrix<G> &X) { return typename colFlipArray<const Matrix<G> >::self(X); }


template<class T>
class row_mirror_matrix_generator : public array_generator<T>
{
  public:
    typedef row_mirror_matrix_generator self;
    typedef array_generator<T> base;
    ARRAY_BASE_TYPES

  public:
    row_mirror_matrix_generator(array_type &x) : base(x) {}

    using base::array;
         
    inline reference       operator[](const index_type &p)       { return array()(-p.i,p.j); }
    inline const_reference operator[](const index_type &p) const { return array()(-p.i,p.j); }

    inline index_type lower_bound() const { return index_type(-array().row_upper_bound(),array().col_lower_bound()); }
};

template<class T>
struct rowMirrorArray 
{
  typedef T array_type;
  typedef row_mirror_matrix_generator<array_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self; 
}; 

template<class G> inline typename rowMirrorArray<      Matrix<G> >::self row_mirror(      Matrix<G> &X) { return typename rowMirrorArray<      Matrix<G> >::self(X); }
template<class G> inline typename rowMirrorArray<const Matrix<G> >::self row_mirror(const Matrix<G> &X) { return typename rowMirrorArray<const Matrix<G> >::self(X); }


template<class T>
class col_mirror_matrix_generator : public array_generator<T>
{
  public:
    typedef col_mirror_matrix_generator self;
    typedef array_generator<T> base;
    ARRAY_BASE_TYPES

  public:
    inline col_mirror_matrix_generator(array_type &x) : base(x) {}

    using base::array;          
     
    inline reference       operator[](const index_type &p)       { return array()(p.i,-p.j); }
    inline const_reference operator[](const index_type &p) const { return array()(p.i,-p.j); }

    inline index_type lower_bound() const { return index_type(array().row_lower_bound(),-array().col_upper_bound()); }
};

template<class T>
struct colMirrorArray 
{
  typedef T array_type;
  typedef col_mirror_matrix_generator<array_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self; 
}; 

template<class G> inline typename colMirrorArray<      Matrix<G> >::self col_mirror(      Matrix<G> &X) { return typename colMirrorArray<      Matrix<G> >::self(X); }
template<class G> inline typename colMirrorArray<const Matrix<G> >::self col_mirror(const Matrix<G> &X) { return typename colMirrorArray<const Matrix<G> >::self(X); }


template<class V>
class uniform_matrix_generator:public sized_value_matrix_generator<V>
{
  public:
	  typedef uniform_matrix_generator self;
    typedef sized_value_matrix_generator<V> base;
    typedef typename base::value_type      value_type;      
    typedef typename base::reference       reference;       
    typedef typename base::const_reference const_reference; 
    typedef typename base::pointer         pointer;         
    typedef typename base::const_pointer   const_pointer;   
    typedef typename base::index_type      index_type;      
    typedef typename base::size_type       size_type;       
    typedef typename base::difference_type difference_type; 
    typedef typename base::int_type        int_type;
      
  protected:
	  value_type val;
  
  public:
    inline uniform_matrix_generator() : base(), val() {}
    inline uniform_matrix_generator(const size_type &s, const value_type &v) : base(s), val(v) {}
    inline uniform_matrix_generator(int_type m, int_type n, const value_type &v) : base(m,n), val(v) {}
    
    using base::size; 
     
    inline const_reference operator[](const index_type &p)  const { assert(size().withinbounds(p)); return val;}  
};

template<class V>
struct UniformMatrix
{
	  typedef V value_type;
  	typedef uniform_matrix_generator<value_type> generator_type;
	  typedef typename GeneratorArray<generator_type>::self self;
};


template<class V>
class identity_matrix_generator : public square_sized_value_matrix_generator<V>
{
  public:
    typedef identity_matrix_generator self;
    typedef square_sized_value_matrix_generator<V> base;
    typedef typename base::value_type      value_type;      
    typedef typename base::reference       reference;       
    typedef typename base::const_reference const_reference; 
    typedef typename base::pointer         pointer;         
    typedef typename base::const_pointer   const_pointer;   
    typedef typename base::index_type      index_type;      
    typedef typename base::size_type       size_type;       
    typedef typename base::difference_type difference_type; 
    typedef typename base::int_type        int_type;
        
  public:
    inline identity_matrix_generator() : base() {}
    inline identity_matrix_generator(int_type n) : base(n) {}
    
    using base::size;

    inline const_reference operator[](const index_type &p) const { assert(size().withinbounds(p)); if (p.i==p.j)  return 1; else return 0; }
};

template<class V> 
struct IdentityMatrix   
{
  typedef V value_type;
  typedef identity_matrix_generator<value_type> generator_type; 
  typedef typename GeneratorArray<generator_type>::self self;
}; 


template<class T,int C=0>
class trans_matrix_generator : public array_generator<T,C>
{
  public:
    typedef trans_matrix_generator self;
    typedef array_generator<T> base;
    ARRAY_BASE_TYPES

  public:
    inline trans_matrix_generator(array_type &x) : base(x) {}

    inline array_type       &array()       { return base::array(); }
    inline const array_type &array() const { return base::array(); }            

    inline reference         operator[](const index_type &p)       { return array()[trn(p)]; }
    inline const_reference   operator[](const index_type &p) const { return array()[trn(p)]; }

    inline size_type  size() const { return trn(array().size()); }
    inline index_type lower_bound() const { return trn(array().lower_bound()); }
};

template<class T,int C=0>
struct TransMatrix
{
  typedef T array_type;
  typedef trans_matrix_generator<array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G> inline typename TransMatrix<      Matrix<G> >::self trn(      Matrix<G> &X) { return typename TransMatrix<      Matrix<G> >::self(X); }
template<class G> inline typename TransMatrix<const Matrix<G> >::self trn(const Matrix<G> &X) { return typename TransMatrix<const Matrix<G> >::self(X); }



template<class A>
class rot90_matrix_generator : public array_generator<A>
{
  public:
    typedef rot90_matrix_generator self;
    typedef array_generator<A> base;
    ARRAY_BASE_TYPES

    template<class A2> struct iterator_rebind       { typedef typename shift_iterator_rebind<      A2>::other other; };
    template<class A2> struct const_iterator_rebind { typedef typename shift_iterator_rebind<const A2>::other other; };

    using base::array;

  public:
    inline rot90_matrix_generator(array_type &x) : base(x) {}

    inline reference         operator[](const index_type &p)       { return array()(p.j,-p.i); }
    inline const_reference   operator[](const index_type &p) const { return array()(p.j,-p.i); }

    inline size_type  size() const { return trn(array().size()); }
    inline index_type lower_bound() const { return index_type(-array().col_upper_bound(),array().row_lower_bound()); }
};

template<class A>
struct Rot90Matrix
{
  typedef A array_type;
  typedef rot90_matrix_generator<array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G> inline typename Rot90Matrix<      Matrix<G> >::self rot90(      Matrix<G> &X) { return typename Rot90Matrix<      Matrix<G> >::self(X); }
template<class G> inline typename Rot90Matrix<const Matrix<G> >::self rot90(const Matrix<G> &X) { return typename Rot90Matrix<const Matrix<G> >::self(X); }



template<class A>
class rot270_matrix_generator : public array_generator<A>
{
  public:
    typedef rot270_matrix_generator self;
    typedef array_generator<A> base;
    ARRAY_BASE_TYPES

    template<class A2> struct iterator_rebind       { typedef typename shift_iterator_rebind<      A2>::other other; };
    template<class A2> struct const_iterator_rebind { typedef typename shift_iterator_rebind<const A2>::other other; };

    using base::array;

  public:
    inline rot270_matrix_generator(array_type &x) : base(x) {}

    inline reference         operator[](const index_type &p)       { return array()(-p.j,p.i); }
    inline const_reference   operator[](const index_type &p) const { return array()(-p.j,p.i); }

    inline size_type  size() const { return trn(array().size()); }
    inline index_type lower_bound() const { return index_type(array().col_lower_bound(),-array().row_upper_bound()); }
};

template<class A>
struct Rot270Matrix
{
  typedef A array_type;
  typedef rot270_matrix_generator<array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

template<class G> inline typename Rot270Matrix<      Matrix<G> >::self rot270(      Matrix<G> &X) { return typename Rot270Matrix<      Matrix<G> >::self(X); }
template<class G> inline typename Rot270Matrix<const Matrix<G> >::self rot270(const Matrix<G> &X) { return typename Rot270Matrix<const Matrix<G> >::self(X); }



template<class A1, class A2>
class multiply_matrix_generator : public two_arrays_value_generator<const A1,const A2>
{
  public:
    typedef multiply_matrix_generator self;
    typedef two_arrays_value_generator<const A1,const A2> base;
    TWO_ARRAYS_BASE_TYPES

  public:
    inline multiply_matrix_generator(first_array_type &x, second_array_type &y) : base(x,y) {}
    
    using base::first_array;
    using base::second_array;

    inline const_reference   operator[] (const index_type &p) const { return inner_product(row(first_array(),p.i),col(second_array(),p.j)); }

    inline size_type  size() const { return size_type(first_array().nrows(),second_array().ncols() ); }
    inline index_type lower_bound() const { return index_type(first_array().row_lower_bound(),second_array().col_lower_bound()); }
};

template<class A1, class A2>
struct MultiplyMatrix
{
  typedef A1 first_array_type;
  typedef A2 second_array_type;
  typedef multiply_matrix_generator<first_array_type, second_array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

//{unsecret}
//Summary: Matrix multiplication
//Parameters:
//  X : The first matrix
//  Y : The second array
//Return: An array representing the multiplied matrices
//Example:
//  DenseMatrix<int>::self X(2,2,"1 2 3 4");
//  DenseMatrix<int>::self Y(2,2,"5 6 7 8");
//  DenseVector<int>::self Z(2,"5 6");
//  cout << mul(X,Y) << endl;
//  cout << mul(X,Z) << endl;
//See: ^operator*^
template<class G1,class G2> inline typename MultiplyMatrix<const Matrix<G1>,const Matrix<G2> >::self mul(const Matrix<G1> &X,const Matrix<G2> &Y) { return typename MultiplyMatrix<const Matrix<G1>, const Matrix<G2> >::self(X,Y); }



template<class A2, class A1>
class multiply_matrix_vector_generator : public two_arrays_value_generator<const A1,const A2>
{
  public:
    typedef multiply_matrix_vector_generator self;
    typedef two_arrays_value_generator<const A1,const A2> base;
    TWO_ARRAYS_BASE_TYPES

  public:
    inline multiply_matrix_vector_generator(second_array_type &y, first_array_type &x) : base(x,y) {}
    
    using base::first_array;
    using base::second_array;
    
    inline const_reference   operator[] (const index_type &i) const { return inner_product(row(second_array(),i),first_array()); }

    inline size_type  size() const { return second_array().nrows(); }
    inline index_type lower_bound() const { return second_array().row_lower_bound(); }
};

template<class A1, class A2>
struct MultiplyMatrixVector
{
  typedef A1 first_array_type;
  typedef A2 second_array_type;
  typedef multiply_matrix_vector_generator<first_array_type, second_array_type> generator_type;
  typedef typename GeneratorArray<generator_type>::self self;
};

//{unsecret}
template<class G1,class G2> inline typename MultiplyMatrixVector<const Matrix<G1>,const Vector<G2> >::self mul(const Matrix<G1> &X,const Vector<G2> &Y) { return typename MultiplyMatrixVector<const Matrix<G1>, const Vector<G2> >::self(X,Y); }



//}

#endif

