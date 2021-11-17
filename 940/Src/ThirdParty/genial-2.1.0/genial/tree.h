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

#ifndef TREE_H
#define TREE_H

#include "util.h"

#include "memory.h"

//namespace genial
//{

using namespace std;

template<class Node> class bidirectional_tree_node_reference;
template<class Node> class tree;

template<class Node>
struct node_traits
{
  typedef Node node_type;
  typedef typename node_type::value_type value_type;
  typedef typename node_type::reference reference;
  typedef typename node_type::const_reference const_reference;
  typedef typename node_type::pointer pointer;
  typedef typename node_type::const_pointer const_pointer;

  typedef typename node_type::child_type child_type;
  typedef typename node_type::children_type children_type;
  typedef typename node_type::child_iterator       child_iterator;
  typedef typename node_type::const_child_iterator const_child_iterator;

  typedef typename node_type::size_type size_type;
  typedef typename node_type::difference_type difference_type;

  typedef typename node_type::iterator iterator;
  typedef typename node_type::const_iterator const_iterator;
  typedef typename node_type::reverse_iterator reverse_iterator;
  typedef typename node_type::const_reverse_iterator const_reverse_iterator;

  typedef typename node_type::value_iterator value_iterator;
  typedef typename node_type::const_value_iterator const_value_iterator;
  typedef typename node_type::reverse_value_iterator reverse_value_iterator;
  typedef typename node_type::const_reverse_value_iterator const_reverse_value_iterator;
};

template<class Node>
struct node_traits<const Node>
{
  typedef const Node node_type;
  typedef const typename node_type::value_type value_type;
  typedef typename node_type::const_reference reference;
  typedef typename node_type::const_reference const_reference;
  typedef typename node_type::const_pointer pointer;
  typedef typename node_type::const_pointer const_pointer;

  typedef const typename node_type::child_type child_type;
  typedef const typename node_type::children_type children_type;
  typedef typename node_type::const_child_iterator child_iterator;
  typedef typename node_type::const_child_iterator const_child_iterator;

  typedef typename node_type::size_type size_type;
  typedef typename node_type::difference_type difference_type;

  typedef typename node_type::const_iterator iterator;
  typedef typename node_type::const_iterator const_iterator;
  typedef typename node_type::const_reverse_iterator reverse_iterator;
  typedef typename node_type::const_reverse_iterator const_reverse_iterator;

  typedef typename node_type::const_value_iterator value_iterator;
  typedef typename node_type::const_value_iterator const_value_iterator;
  typedef typename node_type::const_reverse_value_iterator reverse_value_iterator;
  typedef typename node_type::const_reverse_value_iterator const_reverse_value_iterator;
};

template<class Tree>
struct tree_traits : public node_traits<typename Tree::node_type> 
{
  typedef Tree tree_type;
};

template<class Tree>
struct tree_traits<const Tree> : public node_traits<const typename Tree::node_type> 
{
  typedef const Tree tree_type;
};


template<class Node>
class node_iterator
{
  public:
    typedef node_iterator self;

    typedef Node node_type;
    typedef typename node_type::value_type value_type;
    typedef typename node_type::reference reference;
    typedef typename node_type::const_reference const_reference;
    typedef typename node_type::pointer pointer;
    typedef typename node_type::const_pointer const_pointer;

    typedef typename node_type::child_type child_type;
    typedef typename node_type::children_type children_type;
    typedef typename node_type::child_iterator       child_iterator;
    typedef typename node_type::const_child_iterator const_child_iterator;

    typedef typename node_type::size_type size_type;
    typedef typename node_type::difference_type difference_type;

    typedef typename node_type::iterator iterator;
    typedef typename node_type::const_iterator const_iterator;
    typedef typename node_type::reverse_iterator reverse_iterator;
    typedef typename node_type::const_reverse_iterator const_reverse_iterator;

    typedef typename node_type::value_iterator value_iterator;
    typedef typename node_type::const_value_iterator const_value_iterator;
    typedef typename node_type::reverse_value_iterator reverse_value_iterator;
    typedef typename node_type::const_reverse_value_iterator const_reverse_value_iterator;

    typedef typename iterator_traits<child_iterator>::iterator_category iterator_category;

  protected:
    child_iterator it;

  public:
    node_iterator() : it() {}
    explicit node_iterator(const child_iterator &i) : it(i) {}

    node_iterator(const self &x) : it(x.it) {}
//    template<class Node2> node_iterator(const node_iterator<Node2> &x) : it(x.base()) {}

    self &operator=(const self &r) { it=r.it; return *this; }
//    template<class A> self &operator=(const node_iterator<A> &r) { it=r.base(); return *this; }

    child_iterator       &base()       { return it; }
    const child_iterator &base() const { return it; }

    reference       operator* ()       { return reference(iterator(it)); }
    const_reference operator* () const { return  *it; }
    pointer         operator->()       { node_type       &n=*it; return &n; } 
    const_pointer   operator->() const { const node_type &n=*it; return &n; }

    self &operator++() { ++it; return *this; }
    self &operator--() { --it; return *this; }
    self  operator++(int) { self t=*this; ++*this; return t; }
    self  operator--(int) { self t=*this; --*this; return t; }
    self  operator+ (const difference_type &n) const { return self(it+n); }
    self  operator- (const difference_type &n) const { return self(it-n); }
    self &operator+=(const difference_type &n) { it+=n; return *this; }
    self &operator-=(const difference_type &n) { it-=n; return *this; }

    bool operator==(const self &r) const { return it==r.it; }
    bool operator!=(const self &r) const { return it!=r.it; }  

    bool operator!() const { return !it; }
};

template<class Node>
class node_iterator<const Node>
{
  public:
    typedef node_iterator self;

    typedef const Node node_type; // laisser le const car sinon node_value() ne marche pas pour const_iterator
    typedef const typename node_type::value_type value_type;
    typedef typename node_type::const_reference reference;
    typedef typename node_type::const_reference const_reference;
    typedef typename node_type::const_pointer pointer;
    typedef typename node_type::const_pointer const_pointer;

    typedef const typename node_type::child_type child_type;
    typedef const typename node_type::children_type children_type;
    typedef typename node_type::const_child_iterator child_iterator;
    typedef typename node_type::const_child_iterator const_child_iterator;

    typedef typename node_type::size_type size_type;
    typedef typename node_type::difference_type difference_type;

    typedef typename node_type::const_iterator iterator;
    typedef typename node_type::const_iterator const_iterator;
    typedef typename node_type::const_reverse_iterator reverse_iterator;
    typedef typename node_type::const_reverse_iterator const_reverse_iterator;

    typedef typename node_type::const_value_iterator value_iterator;
    typedef typename node_type::const_value_iterator const_value_iterator;
    typedef typename node_type::const_reverse_value_iterator reverse_value_iterator;
    typedef typename node_type::const_reverse_value_iterator const_reverse_value_iterator;

    typedef typename iterator_traits<child_iterator>::iterator_category iterator_category;

  protected:
    child_iterator it;

  public:
    node_iterator() : it() {}
    explicit node_iterator(const child_iterator &i) : it(i) {}

    node_iterator(const self &x) : it(x.it) {}
    template<class Node2> node_iterator(const node_iterator<Node2> &x) : it(x.base()) {}

    self &operator=(const self &r) { it=r.it; return *this; }
    template<class A> self &operator=(const node_iterator<A> &r) { it=r.base(); return *this; }

    child_iterator       &base()       { return it; }
    const child_iterator &base() const { return it; }

    reference       operator* ()       { return  *it; }
    const_reference operator* () const { return  *it; }
    pointer         operator->()       { node_type       &n=*it; return &n; } 
    const_pointer   operator->() const { const node_type &n=*it; return &n; }

    self &operator++() { ++it; return *this; }
    self &operator--() { --it; return *this; }
    self  operator++(int) { self t=*this; ++*this; return t; }
    self  operator--(int) { self t=*this; --*this; return t; }
    self  operator+ (const difference_type &n) const { return self(it+n); }
    self  operator- (const difference_type &n) const { return self(it-n); }
    self &operator+=(const difference_type &n) { it+=n; return *this; }
    self &operator-=(const difference_type &n) { it-=n; return *this; }

    bool operator==(const self &r) const { return it==r.it; }
    bool operator!=(const self &r) const { return it!=r.it; }  

    bool operator!() const { return !it; }
};

template<class Node> 
inline function_iterator<node_iterator<Node>, value_function<const typename node_iterator<Node>::node_type> > //const ou pas const ?
node_value(const node_iterator<Node> &i) 
{
  return function_iterator<node_iterator<Node>, value_function<const typename node_iterator<Node>::node_type> >(i); //const ou pas const ?
}


template<class Node>
class bidirectional_tree_node_reference : public node_traits<Node>
{
  public:
    typedef bidirectional_tree_node_reference self;
    typedef node_traits<Node> base;
    typedef typename base::node_type       node_type;
    typedef typename base::children_type   children_type;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::iterator        iterator;
    typedef typename base::const_iterator  const_iterator;
    typedef typename base::size_type       size_type;

    //typedef typename node_traits<node_type>::iterator iterator;

  private:
    iterator it;

  public:
    bidirectional_tree_node_reference() : it() {}
    explicit bidirectional_tree_node_reference(const iterator &i) : it(i) {}
    bidirectional_tree_node_reference(const self &r) : it(r.it) {}

    self &operator=(const node_type       &x) { node()=x;        return *this; }
    self &operator=(const tree<node_type> &x) { node()=x.root(); return *this; }
    self &operator=(const self            &x) { node()=x.node(); return *this; }

    operator node_type       &()       { return node(); }
    operator const node_type &() const { return node(); }

    node_type       &node()       { return *it.base(); }
    const node_type &node() const { return *it.base(); }

    children_type       &nodes()       { return node().nodes(); }
    const children_type &nodes() const { return node().nodes(); }

    value_type       &value()       { return node().value(); }
    const value_type &value() const { return node().value(); }

    iterator       &parent()       { return node().parent(); }
    const_iterator &parent() const { return node().parent(); }

    iterator       begin()       { return node().begin(); }
    const_iterator begin() const { return node().begin(); }
    iterator       end  ()       { return node().end  (); }
    const_iterator end  () const { return node().end  (); }

    reference       operator[](const size_type &n)       { return node()[n]; }
    const_reference operator[](const size_type &n) const { return node()[n]; }

    reference       front()       { return node().front(); }
    const_reference front() const { return node().front(); }
    reference       back ()       { return node().back (); }
    const_reference back () const { return node().back (); }

    void pop_front() { node().pop_front(); }

    void push_back(const value_type      &x) { node().push_back(node_type(x,        it)); update(); }
    void push_back(const node_type       &x) { node().push_back(node_type(x,        it)); update(); }
//    void push_back(const reference        x) { push_back(x.node()); }
//    void push_back(const_reference        x) { node().push_back(node_type(x,        it)); update(); }

    void push_back(const tree<node_type> &x) { push_back(x.root()); }
    
    template<class V> void push_back(const V &x) { push_back(value_type(x)); } 

    template<class Pred> iterator       find_value_if(Pred pr)       { return node().find_value_if(pr); }
    template<class Pred> const_iterator find_value_if(Pred pr) const { return node().find_value_if(pr); }

    void update(const iterator &p) { node().update(p); } 

  private:
    void update() { iterator e=end(); for (iterator i=begin(); i!=e; ++i) i->update(i); }
};

template<class E,class Tr,class Node> inline basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os, const bidirectional_tree_node_reference<Node> &X) { return os << X.node(); }

template<class Node> inline bidirectional_tree_node_reference<Node> operator<<(      bidirectional_tree_node_reference<Node> X, const typename bidirectional_tree_node_reference<Node>::value_type &x) { X.push_back(x); return X; }
template<class Node> inline bidirectional_tree_node_reference<Node> operator<<(      bidirectional_tree_node_reference<Node> X, const typename bidirectional_tree_node_reference<Node>::node_type  &x) { X.push_back(x); return X; }
template<class Node> inline bidirectional_tree_node_reference<Node> operator<<(      bidirectional_tree_node_reference<Node> X, const tree<Node>                                                   &x) { X.push_back(x); return X; }

template<class Node> inline bidirectional_tree_node_reference<Node> operator>>(const bidirectional_tree_node_reference<Node> X,       typename bidirectional_tree_node_reference<Node>::value_type &x) { x=X.front().value(); X.pop_front(); return X; }
template<class Node> inline bidirectional_tree_node_reference<Node> operator>>(const bidirectional_tree_node_reference<Node> X,       typename bidirectional_tree_node_reference<Node>::node_type  &x) { x=X.front();         X.pop_front(); return X; }
template<class Node> inline bidirectional_tree_node_reference<Node> operator>>(const bidirectional_tree_node_reference<Node> X,       tree<Node>                                                   &x) { x=X.front();         X.pop_front(); return X; }


template<class V>
class bidirectional_tree_node
{
  public:
    typedef bidirectional_tree_node self;

    typedef V value_type;

    typedef self node_type;
    typedef bidirectional_tree_node_reference<node_type> reference;
    typedef const node_type &const_reference;
    typedef node_type *pointer;
    typedef const node_type *const_pointer;

//    typedef node_type child_type;
//    typedef auto_value<node_type> child_type;
    typedef shared_value<node_type> child_type;
    typedef list<child_type> children_type;
    typedef typename children_type::iterator       child_iterator;
    typedef typename children_type::const_iterator const_child_iterator;

    typedef typename children_type::size_type size_type;
    typedef typename children_type::difference_type difference_type;

    typedef node_iterator<self>       iterator;
    typedef node_iterator<const self> const_iterator;

    typedef std::reverse_iterator<      iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    typedef function_iterator<      iterator, value_function<      node_type> > value_iterator;
    typedef function_iterator<const_iterator, value_function<const node_type> > const_value_iterator;

    typedef std::reverse_iterator<      value_iterator> reverse_value_iterator;
    typedef std::reverse_iterator<const_value_iterator> const_reverse_value_iterator;

  protected:
    value_type val;
    iterator par; 
    children_type children;

  public:
    bidirectional_tree_node() : val(), par(), children() {}
    explicit bidirectional_tree_node(const iterator &p) : val(), par(p), children() {}
    bidirectional_tree_node(const value_type &v, const iterator &p) : val(v), par(p), children() {}

    bidirectional_tree_node(const self &r) : val(r.val), par(r.par), children(r.children) {}
    bidirectional_tree_node(const self &r, const iterator &p) : val(r.val), par(p), children(r.children) {}
//    bidirectional_tree_node (const value_type &v) : val(v), par(), children() { }

    self &operator=(const value_type &v) { val=v; children.resize(0); return *this; }
//    self &operator=(const reference   r) { val=r.value(); children=r.nodes(); return *this; }
//    self &operator=(const self &r) { val=r.val; children=r.children; return *this; }

    value_type &      value()       { return val; }
    const value_type &value() const { return val; }

    iterator       &parent()       { return par; }
    const_iterator &parent() const { return par; }

    children_type       &nodes()       { return children; }
    const children_type &nodes() const { return children; }

    size_type size() const { return nodes().size(); };
    bool empty() const { return nodes().empty(); }

    iterator       begin()       { return iterator      (child_iterator      (nodes().begin())); }
    const_iterator begin() const { return const_iterator(const_child_iterator(nodes().begin())); }
    iterator       end  ()       { return iterator      (child_iterator      (nodes().end  ())); }
    const_iterator end  () const { return const_iterator(const_child_iterator(nodes().end  ())); }

    reference       operator[](const size_type &n)       { iterator       it(begin()); advance(it,n); return *it; }
    const_reference operator[](const size_type &n) const { const_iterator it(begin()); advance(it,n); return *it; }

    reference       front()       { return *begin(); }
    const_reference front() const { return *begin(); }
    reference       back ()       { return *--end(); }
    const_reference back () const { return *--end(); }

    void pop_front() { nodes().pop_front(); }

    void push_back(const node_type  &x) { nodes().push_back(child_type(x)); }

    template<class Pred> iterator       find_value_if(Pred pr)       { return std::find_if(node_value(begin()), node_value(end()), pr).base(); }
    template<class Pred> const_iterator find_value_if(Pred pr) const { return std::find_if(node_value(begin()), node_value(end()), pr).base(); }

    void update(const iterator &p) { iterator e=end(); for (iterator i=begin(); i!=e; ++i) i->parent()=p; } // inutil si les iterateurs sont permanents
};

template<class E,class T,class V> inline basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os, const bidirectional_tree_node<V> &x)
{
  os << x.value(); 
  if (x.nodes().size()) os << "->" << x.nodes(); 
  return os; 
}


template<class Tree>
class tree_iterator : public tree_traits<Tree>
{
  public:
    typedef tree_iterator self;
    typedef tree_traits<Tree> base;    
    
    typedef typename base::tree_type       tree_type;
    typedef typename base::node_type       node_type;
    typedef typename base::children_type   children_type;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::iterator        iterator;
    typedef typename base::const_iterator  const_iterator;
    typedef typename base::size_type       size_type;

    typedef bidirectional_iterator_tag iterator_category;

  protected:
    tree_type &tr;
    iterator it;

  public:
    tree_iterator(tree_type &t, const iterator &i) : tr(t), it(i) {}
    tree_iterator(const self &r) : tr(r.tr), it(r.it) {}

    self &operator=(const self &r) { tr=r.tr; it=r.it; }

    tree_type       &tree()       { return tr; }
    const tree_type &tree() const { return tr; }

    node_type       &node()       { return *(it.base()); }
    const node_type &node() const { return *(it.base()); }

    children_type       &nodes()       { return node().nodes(); }
    const children_type &nodes() const { return node().nodes(); }

    size_type      size() const { return node().size(); }
    bool empty() const { return node().empty(); }

    iterator       begin()       { return node().begin(); }
    const_iterator begin() const { return node().begin(); }
    iterator       end  ()       { return node().end  (); }
    const_iterator end  () const { return node().end  (); }

    reference       operator* ()       { return  reference(it); }
    const_reference operator* () const { return  *it; }

    bool operator==(const self &r) const { return it==r.it; }
    bool operator!=(const self &r) const { return it!=r.it; }    
};


template <class Tree>
class basic_bidirectional_tree_iterator : public tree_iterator<Tree>
{
  public:
    typedef basic_bidirectional_tree_iterator self;
    typedef tree_iterator<Tree> base;

    typedef typename base::tree_type       tree_type;
    typedef typename base::node_type       node_type;
    typedef typename base::children_type   children_type;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::iterator        iterator;
    typedef typename base::const_iterator  const_iterator;
    typedef typename base::size_type       size_type;
    
    //typedef typename node_traits<node_type>::iterator iterator;

  public:
    basic_bidirectional_tree_iterator(tree_type &t, const iterator &i) : base(t, i) {}
    basic_bidirectional_tree_iterator(const self &r) : base(r) {}

/*    self &operator++()
    {
      iterator next=begin();
      while (next==end())
      {
        if (it==tree().parent()) return *this;
        next = it;
        it = (*it).parent();
        ++next;
      }
      it = next;
      return *this;
    }*/
};

//{unsecret}
//Summary: Tree Container
template<class Node>
class tree : public node_traits<Node>
{
  public:
    typedef tree self;
    typedef node_traits<Node> base;

    typedef typename base::node_type       node_type;
    typedef typename base::children_type   children_type;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::iterator        iterator;
    typedef typename base::const_iterator  const_iterator;
    typedef typename base::size_type       size_type;
    typedef typename base::child_iterator       child_iterator;
    typedef typename base::const_child_iterator const_child_iterator;

    //typedef typename node_traits<node_type>::iterator iterator;

    typedef basic_bidirectional_tree_iterator<self      > basic_iterator;
    typedef basic_bidirectional_tree_iterator<const self> const_basic_iterator;

    typedef std::reverse_iterator<      basic_iterator> reverse_basic_iterator;
    typedef std::reverse_iterator<const_basic_iterator> const_reverse_basic_iterator;

  private:
    children_type header;

  public:
    tree() : header(1) { (*parent()).push_back(value_type()); }
    explicit tree(const value_type &x) : header(1) { (*parent()).push_back(x); }
    explicit tree(const node_type  &x) : header(1) { (*parent()).push_back(x); }
    tree(const self &r) : header(r.header) {}

    self &operator=(const self &r) { header=r.header; return *this; }
    self &operator=(const value_type &v) { resize(0); root().value()=v; return *this; }

    value_type       &value()       { return root().value(); }
    const value_type &value() const { return root().value(); }

    children_type       &nodes()       { return root().nodes(); }
    const children_type &nodes() const { return root().nodes(); }

    iterator       parent()       { return iterator      (child_iterator      (header.begin())); }
    const_iterator parent() const { return const_iterator(const_child_iterator(header.begin())); }

    reference       root()       { return (*parent()).front(); }
    const_reference root() const { return (*parent()).front(); }

    size_type size() const { return nodes().size(); }
    void resize(const size_type &n) { nodes().resize(n); }

    iterator       begin()       { return root().begin(); }
    const_iterator begin() const { return root().begin(); }
    iterator       end  ()       { return root().end  (); }
    const_iterator end  () const { return root().end  (); }

    reference       operator[](const size_type &n)       { return root()[n]; }
    const_reference operator[](const size_type &n) const { return root()[n]; }

    reference       front()       { return root().front(); }
    const_reference front() const { return root().front(); }
    reference       back ()       { return root().back (); }
    const_reference back () const { return root().back (); }

    void pop_front() { root().pop_front(); }

    void push_back(const value_type      &x) { root().push_back(x); }
    void push_back(const reference        x) { root().push_back(x); }
    void push_back(const_reference        x) { root().push_back(x); }
    void push_back(const tree<node_type> &x) { root().push_back(x); }
    
    template<class Pred> iterator       find_value_if(Pred pr)       { return root().find_value_if(pr); }
    template<class Pred> const_iterator find_value_if(Pred pr) const { return root().find_value_if(pr); }
};

template<class E,class Tr,class Node> inline basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os, const tree<Node> &tr) { os << tr.root(); return os; }

template<class Node> inline tree<Node> &operator<<(      tree<Node> &tr, const typename tree<Node>::value_type &x) { tr.root()<<x; return tr; }
template<class Node> inline tree<Node> &operator<<(      tree<Node> &tr, const typename tree<Node>::node_type  &x) { tr.root()<<x; return tr; }
template<class Node> inline tree<Node> &operator<<(      tree<Node> &tr, const tree<Node>                      &x) { tr.root()<<x; return tr; }

template<class Node> inline tree<Node> &operator>>(const tree<Node> &tr,       typename tree<Node>::value_type &x) { tr.root()>>x; return tr; }
template<class Node> inline tree<Node> &operator>>(const tree<Node> &tr,       typename tree<Node>::node_type  &x) { tr.root()>>x; return tr; }
template<class Node> inline tree<Node> &operator>>(const tree<Node> &tr,       tree<Node>                      &x) { tr.root()>>x; return tr; }

template<class V>
struct BidirectionalTree
{
  typedef bidirectional_tree_node<V> node_type;
  typedef tree<node_type> self;
};


//}

#endif
