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

#ifndef XML_H
#define XML_H

#undef MSXML

#if (defined(__ICL) || defined(_MSC_VER))
#define MSXML
//#define LIBXML
#else
#define LIBXML
#endif


#ifdef MSXML
#import <msxml3.dll>// named_guids
//using namespace MSXML2;
#undef max
#undef min
#elif defined(LIBXML)
#include <ios>
namespace libxml2
{
  #include <libxml/xmlmemory.h>
  #include <libxml/parser.h>
}
#endif

#include "tree.h"
#include <algorithm>
#include <fstream>
#include <map>
#include <cassert>


//namespace genial
//{

using namespace std;

template<class Ch,class Tr,class A> void erase_first_last(basic_string<Ch,Tr,A> &s,const Ch *t)
{
  s.resize(s.find_last_not_of(t)+1);
  s.erase(0,s.find_first_not_of(t));
  return;
};

#ifdef MSXML
inline void check_hr(HRESULT hr) { if (FAILED(hr)) throw error(_com_error(hr).ErrorMessage()); }
static class CoInitializeClass { public: CoInitializeClass() {::CoInitialize(NULL);} ~CoInitializeClass() {::CoUninitialize();} } CoInitialise;
#endif

//{unsecret}
//{group:XML Elements}
//Summary: Base class for XML Elements
class XMLElemBase
{
	public:
		typedef XMLElemBase self;

		typedef string name_type;
		typedef string data_type;
		typedef pair<name_type, data_type> value_type;
		typedef string attribute_first_type;
		typedef string attribute_second_type;
		typedef pair<attribute_first_type,attribute_second_type> attribute_type;
		typedef map <attribute_first_type,attribute_second_type> attributes_type;

	public:
		XMLElemBase() {}
	
	  //{secret}
    virtual self *copy() const { return new self(*this); }

	  //Summary: Name access
		virtual name_type       &name()       { throw error("Problem"); }       
		virtual const name_type &name() const { throw error("Problem"); } 

	  //Summary: Data access
		virtual data_type       &data()       { throw error("Problem"); }       
		virtual const data_type &data() const { throw error("Problem"); }

	  //Summary: Attributes access
		virtual attributes_type       &attributes()       { throw error("Problem"); }   
		virtual const attributes_type &attributes() const { throw error("Problem"); }

    //{secret}
		virtual string to_string()     const { return "<Problem>";  }
		//{secret}
		virtual string to_end_string() const { return "</Problem>"; }
};

inline XMLElemBase *copy(const XMLElemBase *p) { return p->copy(); }
string to_string    (const XMLElemBase &x) { return x.to_string    (); }
string to_end_string(const XMLElemBase &x) { return x.to_end_string(); }


template<class E,class T> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os, const XMLElemBase &x) { return os<<to_string(x); }


//{unsecret}
//{group:XML Elements}
//Summary: XML Comment
class XMLComment : public XMLElemBase
{
  public:
		typedef XMLComment self;

		typedef string name_type;
		typedef string data_type;
		typedef pair<name_type, data_type> value_type;

  public:
    static const name_type _name;
    
  protected:
    data_type _data;

  public:
		XMLComment() {}

		explicit XMLComment(const char      *n) : _data(n) {}
		explicit XMLComment(const data_type &n) : _data(n) {}
	
		template<class A> XMLComment(const A &a) : _data(data_type(a)){}
		
		XMLComment(const self &x) : _data(x.data()) {}
		
	  virtual self *copy() const { return new self(*this); }

		const name_type &name() const { return _name; }

		data_type       &data()       { return _data; }
		const data_type &data() const { return _data; }

		virtual string to_string() const {  return "<!--" + data() + "-->"; }
    virtual string to_end_string() const { return ""; }

#ifdef MSXML
		explicit XMLComment(MSXML2::IXMLDOMNode *x) { data() = x->Gettext(); }
#elif defined(LIBXML)
	  explicit XMLComment(libxml2::xmlNode *x){ data() = (char *)x->content; }   
#endif

};

const XMLComment::name_type XMLComment::_name = "#comment";


//{unsecret}
//{group:XML Elements}
//Summary: XML start declaration
class XMLStartDeclaration : public XMLElemBase
{
  public:
		typedef XMLStartDeclaration self;

		typedef string name_type;
		typedef string data_type;
		typedef pair<name_type, data_type> value_type;
		typedef string attribute_second_type;
		typedef pair<attribute_first_type,attribute_second_type> attribute_type;
		typedef map <attribute_first_type,attribute_second_type> attributes_type;

  protected:
	  attributes_type _attr;
    static const data_type _data;
    static const name_type _name;
  
  public:
		XMLStartDeclaration() {}
		
		template<class A1,class B1> XMLStartDeclaration(const pair<A1,B1> &p1                                              ) : _attr() { attributes()[attribute_first_type(p1.first)]=attribute_second_type(p1.second); }
		template<class A1,class B1> XMLStartDeclaration(const pair<A1,B1> &p1, const pair<A1,B1> &p2                       ) : _attr() { attributes()[attribute_first_type(p1.first)]=p1.second; attributes()[attribute_first_type(p2.first)]=p2.second; }
		template<class A1,class B1> XMLStartDeclaration(const pair<A1,B1> &p1, const pair<A1,B1> &p2, const pair<A1,B1> &p3) : _attr() { attributes()[attribute_first_type(p1.first)]=p1.second; attributes()[attribute_first_type(p2.first)]=p2.second; attributes()[attribute_first_type(p3.first)]=p3.second; }

	  virtual self *copy() const { return new self(*this); }
		
		const name_type &name() const { return _name; }

		const data_type &data() const { return _data; }

		attributes_type       &attributes()       { return _attr; }
		const attributes_type &attributes() const { return _attr; }

		virtual string to_string() const
		{
		  attributes_type::const_iterator it;
      it = attributes().find("version"); if (it==attributes().end()) return "";
			string s;     
      s+="<?"+data()+" version="+"\""+(*it).second+"\"";
      it = attributes().find("encoding"  ); if(it!=attributes().end()) s+=" "+(*it).first + "=" + "\"" +(*it).second + "\"";
      it = attributes().find("standalone"); if(it!=attributes().end()) s+=" "+(*it).first + "=" + "\"" +(*it).second + "\"";
			s+="?>";
			return s;
		}
	  
    virtual string to_end_string() const
    {
      if(attributes().find("version")==attributes().end()) return ""; 
      return "\n";
    }

  public:
#ifdef MSXML
		explicit XMLStartDeclaration(MSXML2::IXMLDOMNode *x)
		{
			MSXML2::IXMLDOMNamedNodeMapPtr pAttr = x->Getattributes();
			if (pAttr) 
				for (int i=0; i<pAttr->Getlength(); ++i) 
					_attr.insert(attribute_type((const char *)pAttr->Getitem(i)->GetnodeName(),(const char *)pAttr->Getitem(i)->Gettext()));
     }
#elif defined(LIBXML)
    explicit XMLStartDeclaration(libxml2::xmlNode *x)
    {
      _attr.insert(attributes_type::value_type("version",(char *)x->doc->version));
      if(x->doc->encoding)   _attr.insert(attributes_type::value_type("encoding",(char *)x->doc->encoding));
      if(x->doc->standalone) _attr.insert(attributes_type::value_type("standalone","yes"));
      else                   _attr.insert(attributes_type::value_type("standalone", "no"));     
    } 
#endif

};

const XMLStartDeclaration::name_type XMLStartDeclaration::_name = "#declaration";
const XMLStartDeclaration::data_type XMLStartDeclaration::_data = "xml";


//{unsecret}
//{group:XML Elements}
//Summary: XML declaration
class XMLDeclaration : public XMLElemBase
{
  protected:
	  static const attributes_type _attr;
    data_type _data;
    name_type _name;
    
  public:
		typedef XMLDeclaration self;

		typedef string name_type;
		typedef string data_type;
		typedef pair<name_type, data_type> value_type;
		typedef string attribute_second_type;
		typedef pair<attribute_first_type,attribute_second_type> attribute_type;
		typedef map <attribute_first_type,attribute_second_type> attributes_type;

  public:
		XMLDeclaration() {}
		explicit XMLDeclaration(const char      *n) : _name(n), _data() {}
		explicit XMLDeclaration(const name_type &n) : _name(n), _data() {}
		XMLDeclaration(const char      *n, const char      *d) : _name(n), _data(d) {}
		XMLDeclaration(const name_type &n, const data_type &d) : _name(n), _data(d) {}

		template<class A, class B> XMLDeclaration(const A &n, const B &d) : _name(name_type(n)), _data(data_type(d)) {}

	  virtual self *copy() const { return new self(*this); }

		name_type       &name()       { return _name; }
		const name_type &name() const { return _name; }

		data_type       &data()       { return _data; }
		const data_type &data() const { return _data; }

		const attributes_type &attributes() const { return _attr; }
		
		virtual string to_string() const { string s; s+="<?"+name() + " " + data() + "?>"; return s; }
	  virtual string to_end_string() const {return "";}
  
  public:
#ifdef MSXML
		explicit XMLDeclaration(MSXML2::IXMLDOMNode *x) {	name() = x->nodeName; data() = x->text; }     
#elif defined(LIBXML)
    explicit XMLDeclaration(libxml2::xmlNode *x) { name() = (char *)x->name; data() = (char *)x->content; } 
#endif

};

const XMLDeclaration::attributes_type XMLDeclaration::_attr;


//{unsecret}
//{group:XML Elements}
//Summary: XML element
class XMLElem : public XMLElemBase
{
public:
	typedef XMLElem self;
	typedef XMLElemBase base;

  protected:
    value_type _val;
    attributes_type _attr;

  public:
    XMLElem() {}

    explicit XMLElem(const char      *n) : _val(value_type(name_type(n),data_type())), _attr() {}
    explicit XMLElem(const name_type &n) : _val(value_type(name_type(n),data_type())), _attr() {}
    
    explicit XMLElem(const XMLComment &x) : _val(value_type(x.name(),x.data())), _attr() {}

    template<class A,        class A1,class B1> XMLElem(const A &a,             const pair<A1,B1> &p1                                              ) : _val(value_type(name_type(a),data_type( ))),  _attr() { attributes()[attribute_first_type(p1.first)]=attribute_second_type(p1.second); }
    template<class A,        class A1,class B1> XMLElem(const A &a,             const pair<A1,B1> &p1, const pair<A1,B1> &p2                       ) : _val(value_type(name_type(a),data_type( ))),  _attr() { attributes()[attribute_first_type(p1.first)]=p1.second; attributes()[attribute_first_type(p2.first)]=p2.second; }
    template<class A,        class A1,class B1> XMLElem(const A &a,             const pair<A1,B1> &p1, const pair<A1,B1> &p2, const pair<A1,B1> &p3) : _val(value_type(name_type(a),data_type( ))),  _attr() { attributes()[attribute_first_type(p1.first)]=p1.second; attributes()[attribute_first_type(p2.first)]=p2.second; attributes()[attribute_first_type(p3.first)]=p3.second; }

    template<class A,class B>                   XMLElem(const A &a, const B &b                                                                     ) : _val(value_type(name_type(a),data_type(b))), _attr() {}
    template<class A,class B,class A1,class B1> XMLElem(const A &a, const B &b, const pair<A1,B1> &p1                                              ) : _val(value_type(name_type(a),data_type(b))), _attr() { attributes()[attribute_first_type(p1.first)]=attribute_second_type(p1.second); }
    template<class A,class B,class A1,class B1> XMLElem(const A &a, const B &b, const pair<A1,B1> &p1, const pair<A1,B1> &p2                       ) : _val(value_type(name_type(a),data_type(b))), _attr() { attributes()[attribute_first_type(p1.first)]=p1.second; attributes()[attribute_first_type(p2.first)]=p2.second; }
    template<class A,class B,class A1,class B1> XMLElem(const A &a, const B &b, const pair<A1,B1> &p1, const pair<A1,B1> &p2, const pair<A1,B1> &p3) : _val(value_type(name_type(a),data_type(b))), _attr() { attributes()[attribute_first_type(p1.first)]=p1.second; attributes()[attribute_first_type(p2.first)]=p2.second; attributes()[attribute_first_type(p3.first)]=p3.second; }

    self  &operator=(const XMLComment &x) { name() = x.name(); data() = x.data(); _attr.clear(); return *this; }
	
	  virtual self *copy() const { return new self(*this); }

    XMLElem(const self &x) : _val(x._val), _attr(x._attr) {}

    name_type       &name()       { return _val.first; }
    const name_type &name() const { return _val.first; }

    data_type       &data()       { return _val.second; }
    const data_type &data() const { return _val.second; }

    attributes_type       &attributes()       { return _attr; }
    const attributes_type &attributes() const { return _attr; }

	  virtual string to_string() const
	  {
	    if (name()=="#text") return data();
	    if (name()==XMLComment::_name) return "<!--" + data() + "-->"; 
		  string s;
		  s+="<"+name();
		  for (attributes_type::const_iterator i=attributes().begin(); i!=attributes().end(); ++i)
		  	s+=" "+(*i).first+"="+"\""+(*i).second+"\"";
		  s+=">";
		  s+=data();
		  return s;
	  }

	  virtual string to_end_string() const 
	  {
	    if (name()=="#text") return "";
	    if (name()==XMLComment::_name) return "";
	    string s; s+="</" + name() + ">"; return s; 
	  }

#ifdef MSXML
		explicit XMLElem(MSXML2::IXMLDOMNode *x)
		{
			name()=x->GetnodeName();
			if (x->nodeType==MSXML2::NODE_TEXT) data()=x->Gettext();
			MSXML2::IXMLDOMNamedNodeMapPtr pAttr = x->Getattributes();
			if (pAttr) 
				for (int i=0; i<pAttr->Getlength(); ++i) 
					_attr.insert(attribute_type((const char *)pAttr->Getitem(i)->GetnodeName(),(const char *)pAttr->Getitem(i)->Gettext()));
			}
#elif defined(LIBXML)
	  explicit XMLElem(libxml2::xmlNode *x)
	  {
      name()=(char *)x->name;
      if(libxml2::xmlNodeIsText(x) && !libxml2::xmlIsBlankNode(x))
      {
        string s((char *)x->content);
        erase_first_last(s," \n\t");
        data() = s;
      }
      for (libxml2::xmlAttrPtr pAttr=x->properties; pAttr; pAttr=pAttr->next)
        _attr.insert(attribute_type((char *)pAttr->name,(char *)pAttr->children->content));
	  } 
#endif

};	
	
inline XMLElem *copy(const XMLElem *p) { return p->copy(); }
string to_string    (const XMLElem &x) { return x.to_string    (); }
string to_end_string(const XMLElem &x) { return x.to_end_string(); }

template<class E,class T> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os, const XMLElem &x) { return os<<to_string(x); }


template<>
struct value_traits<XMLElemBase>
{
  typedef XMLElemBase       value_type;
  typedef value_type       &reference;
  typedef const value_type &const_reference;
  typedef value_type       *pointer;
  typedef const value_type *const_pointer;

  typedef value_type::name_type name_type;
  typedef value_type::data_type data_type;
  typedef value_type::attribute_first_type attribute_first_type;
  typedef value_type::attribute_second_type attribute_second_type;

  typedef value_type::attribute_type attribute_type;
  typedef value_type::attributes_type attributes_type;
};

template<>
struct value_traits<const XMLElemBase>
{
  typedef const XMLElemBase value_type;
  typedef const value_type &reference;
  typedef const value_type &const_reference;
  typedef const value_type *pointer;
  typedef const value_type *const_pointer;

  typedef const value_type::name_type name_type;
  typedef const value_type::data_type data_type;
  typedef const value_type::attribute_first_type attribute_first_type;
  typedef const value_type::attribute_second_type attribute_second_type;

  typedef const value_type::attribute_type attribute_type;
  typedef const value_type::attributes_type attributes_type;
};


template<class V> class xml_tree;

template<class V> xml_tree<V> to_tree(const xml_tree<V> &x) { return x; }

//{unsecret}
//{group:XML Classes}
//Summary: Pseudo-reference of a XML node
//Remark:
//  This class is given back by non-const node access on xml_tree.
//  Use it like a xml_node class.
//Example:
//  XMLTree xml;
//  xml.push_back("root");
//  XMLTree::reference x = xml.front();
//  x.push_back("first");
//  x.push_back("second");
//See: ^xml_node^
template<class NodeType>
class xml_reference : public bidirectional_tree_node_reference<NodeType>
{
  public:
    typedef xml_reference self;
    typedef bidirectional_tree_node_reference<NodeType> base;

    typedef typename base::node_type       node_type;
    typedef typename base::children_type   children_type;
    typedef typename base::value_type      value_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::iterator        iterator;
    typedef typename base::const_iterator  const_iterator;
    typedef typename base::size_type       size_type;

    typedef typename value_type::name_type       name_type;
    typedef typename value_type::data_type       data_type;
    typedef typename value_type::attribute_type  attribute_type;
    typedef typename value_type::attributes_type attributes_type;

    using base::node;
    using base::value;
    using base::begin;
    using base::end;

  public:
    explicit xml_reference(const iterator &i) : base(i) {}
    xml_reference(const self &x) : base(x) {}

    self &operator=(const node_type &x) { node()=x; return *this; }
    self &operator=(const string    &s) { data()=s; return *this; }
    self &operator=(const char *     s) { data()=s; return *this; }
    self &operator=(const value_type &x) { node()=x; return *this; }
    template<class V> self &operator=(const V &x) { node()=value_type(x); return *this; }
    template<class V> self &operator=(const xml_tree<V> &x) { static_cast<base &>(*this)=x; return *this; }
    
    //{secret}  
    typename value_traits<value_type>::reference       value2()       { return typename value_traits<value_type>::reference      (value()); }  
    //{secret}  
    typename value_traits<value_type>::const_reference value2() const { return typename value_traits<value_type>::const_reference(value()); }  
      
    //{noAutoLink}
    //Summary: Name access
    name_type       &name()       { return value2().name(); }
    const name_type &name() const { return value2().name(); }

    //Summary: Data access
    data_type       &data()       { return value2().data(); }
    const data_type &data() const { return value2().data(); }

    //Summary: Attributes access
    attributes_type       &attributes()       { return value2().attributes(); }
    const attributes_type &attributes() const { return value2().attributes(); }

    //{noAutoLink}
    //Summary: Searches the first node
    iterator       find(const name_type &n)       { return node().find(n); }
    const_iterator find(const name_type &n) const { return node().find(n); }
    template<class A1,class B1> iterator       find(const name_type &n, const pair<A1,B1> &a)       { return node().find(n,a); }
    template<class A1,class B1> const_iterator find(const name_type &n, const pair<A1,B1> &a) const { return node().find(n,a); }

    //Summary: Searches the first node which name satifies the predicate
    template<class Pred> iterator       find_name_if(Pred pr)       { return node().find_name_if(pr); }
    template<class Pred> const_iterator find_name_if(Pred pr) const { return node().find_name_if(pr); }

    //Summary: Node access    
    reference       operator[](const size_type &n)       { return static_cast<      base &>(*this)[n]; }
    const_reference operator[](const size_type &n) const { return static_cast<const base &>(*this)[n]; }
    reference       operator[](const name_type &n)       { iterator       it=find(n); if (it==end()) { push_back(value_type(XMLElem(n))); it=--end(); } return *it; }
    const_reference operator[](const name_type &n) const { const_iterator it=find(n); assert(it!=end()); return *it; }
    
    //Summary: Node access
    template<class A1,class B1> reference operator()(const name_type &n, const pair<A1,B1> &a) { iterator it=find(n,a); if (it==end()) { push_back(value_type(XMLElem(n,a))); it=--end(); } return *it; }

    //Summary: Inserts after the last node
    void push_back(const value_type      &v) { base::push_back(v); }
    void push_back(const node_type       &x) { base::push_back(x); }
    void push_back(const tree<node_type> &x) { base::push_back(x.root()); }
    void push_back(const name_type       &n) { base::push_back(XMLElem(n)); }
    template<class A1,class B1> void push_back(const name_type &n, const pair<A1,B1> &a) { push_back(XMLElem(n,a)); }
    template<class V2> void push_back(const V2 &x) { push_back(value_type(x)); }
};


template<class NodeType, class T>
class xml_cast_reference : public xml_reference<NodeType>
{
  public:
    typedef xml_cast_reference self;
    typedef xml_reference<NodeType> base;

    typedef T cast_type;
    
    typedef typename base::value_type value_type;
    typedef typename base::reference reference;

    using base::attributes;
    using base::node;
    using base::data;
    using base::nodes;
    
  public:
    xml_cast_reference(const reference &x) : base(x) {}
   
    xml_cast_reference(const self &r) : base(static_cast<const base &>(r)) {} 

    template<class V> self &operator=(const xml_tree<V> &x) { static_cast<base &>(*this)=x; return *this; }

    self &operator=(const cast_type &x) { xml_tree<value_type> tr=to_tree<value_type>(x); attributes().insert(tr.attributes().begin(),tr.attributes().end()); data()=tr.data(); nodes()=tr.nodes(); return *this; }

    cast_type cast() const { cast_type v; xml_tree<value_type> tr(node()); tr>>v; return v; }

    operator cast_type() const { return cast(); }
};

template<class E,class Tr,class Node,class T> inline basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os, const xml_cast_reference<Node,T> &X) { return os << X.cast(); }


template<class NodeType, class T>
class xml_basic_reference : public xml_cast_reference<NodeType,T>
{
  public:
    typedef xml_basic_reference self;
    typedef xml_cast_reference<NodeType,T> base;

    typedef typename base::value_type value_type;
    typedef typename base::reference reference;
    typedef typename base::cast_type cast_type;

    using base::attributes;
    using base::node;
    using base::data;
    using base::nodes;

  public:
    xml_basic_reference(const reference &x) : base(x) {}

    xml_basic_reference(const self &r) : base(static_cast<const base &>(r)) {}

    template<class U> self &operator=(const U &x) { data()=to_string(x); return *this; }

    operator cast_type() const { cast_type v; data() >> v; return v; }
};


template<class NodeType, class Cont>
class xml_list_reference : public xml_cast_reference<NodeType,Cont>
{
  public:
    typedef xml_list_reference self;
    typedef xml_cast_reference<NodeType,Cont> base;
    
    typedef typename base::reference      reference;
    typedef typename base::cast_type      cast_type;
    typedef typename base::const_iterator const_iterator;

    using base::attributes;
    using base::node;
    using base::data;
    using base::nodes;

  public:
    xml_list_reference(const reference &x) : base(x) {}
    xml_list_reference(const self &r) : base(static_cast<const base &>(r)) {}

    template<class Cont2> self &operator=(const Cont2 &X) { data()=to_string(X.begin(), X.end()); return *this; }
    operator cast_type() const
    {
      //pas testé
      list<typename Cont::value_type> l;
      istringstream iss(data()); 
      typename cast_type::value_type v;
      while (iss>>v, !iss.eof()) l.append(v);
      return l;
    }
};

template<class V>
class name_attribute_pred
{
  public:
    typedef V value_type;
    typedef typename value_type::name_type       name_type;
    typedef typename value_type::attribute_type  attribute_type;
    typedef typename value_type::attributes_type attributes_type;

  public:
    name_type nam;
    attribute_type attr;

  public:
    template<class A1,class B1> name_attribute_pred(const name_type &n, const pair<A1,B1> &a) : nam(n), attr(a) {}

    bool operator()(const value_type &x) 
    { 
      if (x.name()!=nam) return false;

      typename attributes_type::const_iterator it = x.attributes().find(attr.first);
      if (it==x.attributes().end()) return false;

      return (*it).second==attr.second;
    }
};

template<class InIt,class V   > inline InIt std_find   (InIt begin, InIt end, const V &val) { return find   (begin,end,val ); }
template<class InIt,class Pred> inline InIt std_find_if(InIt begin, InIt end, Pred    pred) { return find_if(begin,end,pred); }

//{unsecret}
//{group:XML Classes}
//Summary: Node of a XML Tree
template<class V>
class xml_node
{
  public:
    typedef xml_node self;

    typedef V value_type;
    
    typedef typename value_type::name_type       name_type;
    typedef typename value_type::data_type       data_type;
    typedef typename value_type::attribute_type  attribute_type;
    typedef typename value_type::attributes_type attributes_type;

    typedef self node_type;
    typedef xml_reference<node_type> reference;
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
    V val;
    iterator par; 
    children_type children;

  public:
    xml_node() : val(), par(), children() {}
    explicit xml_node(const iterator &p) : val(), par(p), children() {}
    xml_node(const value_type &v, const iterator &p) : val(v), par(p), children() {}
    xml_node(const char       *s, const iterator &p) : val(value_type(s)), par(p), children() {}
#ifdef MSXML    
    xml_node(MSXML2::IXMLDOMNode *pNode);
#elif defined(LIBXML)
    xml_node(libxml2::xmlNode *pNode);
#endif  

    xml_node(const self &r) : val(r.val), par(r.par), children(r.children) {}
    xml_node(const self &r, const iterator &p) : val(r.val), par(p), children(r.children) {}
   
    self &operator=(const value_type &v) { val=v; children.resize(0); return *this; }

    //{noAutoLink}
    //Summary: Value access
    value_type &      value()       { return val; }
    const value_type &value() const { return val; } 
    
    //{secret}  
    typename value_traits<value_type>::reference       value2()       { return typename value_traits<value_type>::reference      (value()); }  
    //{secret}  
    typename value_traits<value_type>::const_reference value2() const { return typename value_traits<value_type>::const_reference(value()); }  

    //Summary: Parent access
    iterator       &parent()       { return par; }
    const_iterator &parent() const { return par; }

    //{noAutoLink}
    //Summary: Nodes access
    children_type       &nodes()       { return children; }
    const children_type &nodes() const { return children; }

    //Summary: Number of nodes
    size_type size() const { return nodes().size(); };
    
    //Summary: True if no nodes
    bool empty() const { return nodes().empty(); }

    //Summary: Iterator that points at the first node
    iterator       begin()       { return iterator      (child_iterator      (nodes().begin())); }
    const_iterator begin() const { return const_iterator(const_child_iterator(nodes().begin())); }

    //Summary: Iterator that points just beyond the last node
    iterator       end  ()       { return iterator      (child_iterator      (nodes().end  ())); }
    const_iterator end  () const { return const_iterator(const_child_iterator(nodes().end  ())); }

    //Summary: Node access
    reference       operator[](const size_type &n)       { iterator       it(begin()); advance(it,n); return *it; }
    const_reference operator[](const size_type &n) const { const_iterator it(begin()); advance(it,n); return *it; }
    reference       operator[](const name_type &n)       { iterator       it=find(n); if (it==end()) { push_back(value_type(n)); it=--end(); } return *it; }
    const_reference operator[](const name_type &n) const { const_iterator it=find(n); if (it==end()) { throw error("\""+n+"\" not found"); } return *it; }

    //Summary: First node access
    reference       front()       { return *begin(); }
    const_reference front() const { return *begin(); }
    
    //Summary: Last node access
    reference       back ()       { return *--end(); }
    const_reference back () const { return *--end(); }

    //{secret}
    void update(const iterator &p) { iterator e=end(); for (iterator i=begin(); i!=e; ++i) i->parent()=p; } // inutile si les iterateurs sont permanents

    //{noAutoLink}
    //Summary: Name access
    name_type       &name()       { return value2().name(); }
    const name_type &name() const { return value2().name(); }

    //Summary: Data access
    data_type       &data()       { return value2().data(); }
    const data_type &data() const { return value2().data(); }

    //Summary: Attributes access
    attributes_type       &attributes()       { return value2().attributes(); }
    const attributes_type &attributes() const { return value2().attributes(); }

    //Summary: Searches the first node
    iterator       find(const name_type &s)       { return std_find(name_it(node_value(begin())), name_it(node_value(end())), s).base().base(); }
    const_iterator find(const name_type &s) const { return std_find(name_it(node_value(begin())), name_it(node_value(end())), s).base().base(); }
    template<class A1,class B1> iterator       find(const name_type &n, const pair<A1,B1> &a)       { return find_value_if(name_attribute_pred<value_type>(n,a)); }
    template<class A1,class B1> const_iterator find(const name_type &n, const pair<A1,B1> &a) const { return find_value_if(name_attribute_pred<value_type>(n,a)); }

    //Summary: Searches the first node which name satisfies the predicate 
    template<class Pred> iterator       find_name_if(Pred pr)       { return find_if(name_it(node_value(begin())), name_it(node_value(end())), pr).base().base(); }
    template<class Pred> const_iterator find_name_if(Pred pr) const { return find_if(name_it(node_value(begin())), name_it(node_value(end())), pr).base().base(); }

    //Summary: Searches the first node which value satisfies the predicate 
    template<class Pred> iterator       find_value_if(Pred pr)       { return std_find_if(node_value(begin()), node_value(end()), pr).base(); }
    template<class Pred> const_iterator find_value_if(Pred pr) const { return std_find_if(node_value(begin()), node_value(end()), pr).base(); }

    //Summary: Inserts after the last node
    void push_back(const node_type &x) { nodes().push_back(child_type(x));             }
    void push_back(const name_type &n) { nodes().push_back(child_type(value_type(n))); }
    
    //Summary: Removes the first node
    void pop_front() { nodes().pop_front(); }    
};
#ifdef MSXML    
template<class V> xml_node<V>::xml_node(MSXML2::IXMLDOMNode *pNode)
{
  if (pNode->nodeType==MSXML2::NODE_COMMENT)  value() = XMLComment(pNode);
  else  value() = XMLElem(pNode);
    
  for (MSXML2::IXMLDOMNodePtr pChild=pNode->GetfirstChild(); pChild!=NULL; pChild=pChild->GetnextSibling())
  {
    if (pChild->nodeType==MSXML2::NODE_TEXT) data()=pChild->Gettext(); 
    else push_back(node_type(pChild));
  }
}
#elif defined(LIBXML)
template<class V> xml_node<V>::xml_node(libxml2::xmlNode *pNode)
{
  if(pNode->type==libxml2::XML_COMMENT_NODE)  value() = XMLComment(pNode);
  else  value() = XMLElem(pNode);
  
  for (libxml2::xmlNodePtr pChild=pNode->children; pChild; pChild = pChild->next)
  {
    if(!libxml2::xmlNodeIsText(pChild)) push_back(node_type(pChild));
    else
    if (!libxml2::xmlIsBlankNode (pChild))
    {
      string s((char *)pChild->content);
      erase_first_last(s," \n\t");
      data() = s;
    }
  }
}
#endif   

template<class E,class Tr,class V> inline basic_ostream<E,Tr> &operator<<(basic_ostream<E,Tr> &os, const xml_node<V> &x)
{
  os << x.value(); 
  if (x.nodes().size()) os << "->" << x.nodes(); 
  return os; 
}

//{unsecret}
//{group:XML Classes}
//Summary: XML Container
//Example:
//  xml_tree<XMLElem> xml("root");
//  xml.front().push_back(XMLElem("first","1")); // root->first
//  xml.front().push_back("second"); // root->second
//  xml.front().front().push_back("third"); // root->first->third
//  xml.front().front().push_back("forth"); // root->first->forth
//  xml_tree<XMLElem>::reference x=xml.front(); //root->first
//  x["second"].data() = "2";
template<class V>
class xml_tree : public tree<xml_node<V> >
{
  public:
    typedef xml_tree self;
    typedef tree<xml_node<V> > base;

    typedef typename base::value_type      value_type;
    typedef typename base::node_type       node_type;
    typedef typename base::reference       reference;
    typedef typename base::const_reference const_reference;
    typedef typename base::iterator        iterator;
    typedef typename base::const_iterator  const_iterator;
    typedef typename base::size_type       size_type;

    typedef typename value_type::name_type       name_type;
    typedef typename value_type::data_type       data_type;
    typedef typename value_type::attribute_type  attribute_type;
    typedef typename value_type::attributes_type attributes_type;

    using base::root;
    
  private:
    XMLStartDeclaration  startdec;
    list<XMLDeclaration> dec;

  public:
    xml_tree() : base(value_type(XMLElem("Problem"))) {} 
    explicit xml_tree(const char      *n) : base(value_type(XMLElem(n))) {}
    explicit xml_tree(const node_type &x) : base(x) {}
    explicit xml_tree(const name_type &n) : base(value_type(XMLElem(n))) {}
    xml_tree(const name_type &n, const data_type &d) : base(value_type(XMLElem(n,d))) {}
    template<class A1,class B1                  > xml_tree(const name_type &a,                     const pair<A1,B1> &p1                       ) : base(value_type(XMLElem(a,  p1   ))) {}
    template<class A1,class B1,class A2,class B2> xml_tree(const name_type &a,                     const pair<A1,B1> &p1, const pair<A2,B2> &p2) : base(value_type(XMLElem(a,  p1,p2))) {}
    template<class A1,class B1                  > xml_tree(const name_type &a, const data_type &b, const pair<A1,B1> &p1                       ) : base(value_type(XMLElem(a,b,p1   ))) {}
    template<class A1,class B1,class A2,class B2> xml_tree(const name_type &a, const data_type &b, const pair<A1,B1> &p1, const pair<A2,B2> &p2) : base(value_type(XMLElem(a,b,p1,p2))) {}

    //template<class V2> xml_tree(const V2 &x) : base(value_type(x)) {}

    xml_tree(const value_type &x) : base(x) {}

    self &operator=(const value_type &v) { static_cast<base &>(*this)=v; return *this; }

    xml_tree(const self &r) : base(r) {}

    //{noAutoLink}
    //Summary: Name access
    name_type       &name()       { return root().name(); }
    const name_type &name() const { return root().name(); }

    //Summary: Data access
    data_type       &data()       { return root().data(); }
    const data_type &data() const { return root().data(); }

    //Summary: Attributes access
    attributes_type       &attributes()       { return root().attributes(); }
    const attributes_type &attributes() const { return root().attributes(); }

    //Summary: Start declaration access
    XMLStartDeclaration       &start_declaration()       { return startdec; }
    const XMLStartDeclaration &start_declaration() const { return startdec; }

    //{noAutoLink}
    //Summary: Declarations access
    list<XMLDeclaration>       &declarations()       { return dec; }
    const list<XMLDeclaration> &declarations() const { return dec; }

    //{noAutoLink}
    //Summary: Searches the first node
    iterator       find(const name_type &s)       { return root().find(s); }
    const_iterator find(const name_type &s) const { return root().find(s); }
    template<class A1,class B1> iterator       find(const name_type &s, const pair<A1,B1> &a)       { return root().find(s,a); }
    template<class A1,class B1> const_iterator find(const name_type &s, const pair<A1,B1> &a) const { return root().find(s,a); }

    //Summary: Searches the first node which name satisfies the predicate 
    template<class Pred> iterator       find_name_if(Pred pr)       { return root().find_name_if(pr); }
    template<class Pred> const_iterator find_name_if(Pred pr) const { return root().find_name_if(pr); }

    //Summary: Node access
    reference       operator[](const size_type &n)       { return root()[n]; }
    const_reference operator[](const size_type &n) const { return root()[n]; }
    reference       operator[](const name_type &s)       { return root()[s]; }
    const_reference operator[](const name_type &s) const { return root()[s]; }
    
    //Summary: Node access
    template<class A1,class B1> reference operator()(const name_type &s, const pair<A1,B1> &a) { return root()(s,a); }

    //Summary: Inserts after the last node
    void push_back(const value_type &v) { root().push_back(v); }
    void push_back(const name_type  &s) { root().push_back(s); }
    void push_back(const char *s) { push_back(name_type(s)); }
    void push_back(const tree<node_type> &x) { root().push_back(x); }
    void push_back(const xml_tree<value_type> &x) { push_back(static_cast<const tree<node_type> &>(x)); }
    template<class A1,class B1> void push_back(const name_type &s, const pair<A1,B1> &a) { root().push_back(s,a); }
    template<class V2> void push_back(const V2 &v) { push_back(value_type(v)); }

    //template<class A> void push_descriptor(const A &a) { push_back(to_tree<value_type>(a)); }

    //Summary: Loads a file
    void load(const char *s)
    {    
      start_declaration() = XMLStartDeclaration();
      declarations().clear();
#ifdef MSXML
      MSXML2::IXMLDOMDocumentPtr pDoc;
      check_hr(pDoc.CreateInstance(__uuidof(MSXML2::DOMDocument))); //CLSID_DOMDocument
      if (!(bool)pDoc->load(s)) throw error("xml load failed");

      if (pDoc->GetfirstChild())
      {
        for (MSXML2::IXMLDOMNodePtr pChild=pDoc->GetfirstChild(); pChild && pChild->nodeType==MSXML2::NODE_PROCESSING_INSTRUCTION; pChild=pChild->GetnextSibling())
          if(_stricmp(pChild->GetnodeName(),"xml")==0)  
            start_declaration() = XMLStartDeclaration(pChild);
          else  
            declarations().push_back(XMLDeclaration((MSXML2::IXMLDOMNode *)pChild));
        root() = node_type(pDoc->GetdocumentElement());
      }
#elif defined(LIBXML)
      libxml2::xmlDocPtr pDoc;
      pDoc = libxml2::xmlParseFile(s);
      if (!pDoc) throw error("xml load failed");
      
      start_declaration() = XMLStartDeclaration(pDoc->children);

      if (pDoc->children)
      {
        libxml2::xmlNode *node;
        for (node=pDoc->children; node->type==libxml2::XML_PI_NODE; node=node->next)
          declarations().push_back(XMLDeclaration(node));
        root() = node_type(node);
      }
      libxml2::xmlFreeDoc(pDoc);
#endif
    }
};

//{unsecret}
//{group:XML Type Definitions}
//Summary: XML Container without virtual functions
//Remarks:
//  This definition of xml_tree uses a XMLElem class to save each leaf of the XML tree.
//  The execution speed is quicker than VirtualXMLTree, but uses somewhat more memory.
//See: ^VirtualXMLTree^
typedef xml_tree<XMLElem> XMLTree;

//{unsecret}
//{group:XML Type Definitions}
//Summary: XML Container with virtual functions
//Remarks:
//  This definition of xml_tree uses a XMLElemBase class and its derived classes 
//  (XMLElem, XMLComment, XMLDeclaration, XMLStartDeclaration) to save each leaf of the XML Tree.
//  Datas of the leaves are accessed through virtual functions.
//  This is therefore slower than XMLTree, but uses somewhat less memory.
//See: ^XMLTree^
typedef xml_tree<auto_value<XMLElemBase> > VirtualXMLTree;

template<class E,class T,class Node>
basic_ostream<E,T> &node_to_xml_stream(basic_ostream<E,T> &os,const Node &x,const string &indention=string("  "),int level=0)
{
  for(int i=0;i<level;i++) os << indention;
  
  if(x.empty() && x.data().empty()) {  os << "<" << x.name();  for (typename Node::attributes_type::const_iterator it=x.attributes().begin();it!=x.attributes().end();it++)  os << " " << (*it).first << "=" << "\"" << (*it).second << "\"";  os << " />" << os.widen('\n');  return os;  }
  
  os << to_string(x.value());
    
  if(!x.empty())  
  {
    os << os.widen('\n');
    for (typename Node::const_iterator it=x.begin(); it!=x.end(); ++it) node_to_xml_stream(os,*it,indention,level+1);
    for(int i=0;i<level;i++) os << indention;
  }
  os << to_end_string(x.value()) << os.widen('\n');
  return os;
};

template<class V> string to_string(const xml_tree<V> &x)
{ 
  string s;
  s = x.start_declaration().to_string() + x.start_declaration().to_end_string();
  for(list<XMLDeclaration>::const_iterator it= x.declarations().begin();it != x.declarations().end();++it)
    s+=(*it).to_string()+(*it).to_end_string()+"\n";
  ostringstream oss;
  node_to_xml_stream(oss,x.root());
  return s+oss.str();
}

template<class E,class T,class V> basic_ostream<E,T> &operator<<(basic_ostream<E,T> &os, const xml_tree<V> &x)
{
  string s;
  s = x.start_declaration().to_string() + x.start_declaration().to_end_string();
  for(list<XMLDeclaration>::const_iterator it= x.declarations().begin();it != x.declarations().end();it++)
    s+=(*it).to_string()+(*it).to_end_string()+"\n";
  os << s;
  return node_to_xml_stream(os,x.root());
}

//{unsecret}
//{group:XML Classes}
//Summary: XML file
//Example:
//  XMLFile fin("x.xml");
//  XMLTree xml;
//  fin >> xml; //load
//  ...
//  XMLFile fout("y.xml");
//  fout << xml; // store
class XMLFile : public File
{
  public:
    explicit XMLFile(const string &s) : File(s) {}
};

template<class V> XMLFile &operator>>(XMLFile &f, xml_tree<V> &x) { x.load(f.name().c_str()); return f; }
template<class V> XMLFile &operator<<(XMLFile &f, const xml_tree<V> &x) { ofstream os(f.name().c_str()); os<<x; return f; }


//}

#endif

