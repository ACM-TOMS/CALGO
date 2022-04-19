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

#ifndef IDIRSTREAM_H
#define IDIRSTREAM_H

#if (defined(__ICL) || defined(_MSC_VER)) 
#include <io.h>
#include <direct.h>
#else
#include <sys/stat.h>
#include <glob.h>
#endif

#include "array/iterator.h"

//namespace genial
//{

using namespace std;


#if (defined(__ICL) || defined(_MSC_VER))
typedef _finddata_t file_data;
#else
#define _A_SUBDIR 0x10
#define _A_NORMAL 0x00
#define _A_RDONLY 0x01

struct file_data
{
  char name[256];
  unsigned int attrib;
  unsigned int size;
  long int time_create;
  long int time_access;
  long int time_write;
};
#endif

template<class Ch,class Tr,class A> basic_string<Ch,Tr,A> filename_out_of_pathname(const basic_string<Ch,Tr,A> &s) {  return s.substr(s.rfind("/")+1); };
template<class Ch> basic_string<Ch> filename_out_of_pathname(const Ch *s) { return filename_out_of_pathname(basic_string<Ch>(s)); };

#if (defined(__ICL) || defined(_MSC_VER))
class idirstream : public basic_ios<file_data>
{
  public:
    typedef idirstream self;
  
  private:
    long h;
    string s;
  
  public:
    idirstream() : h(-1), s("*") {}
    explicit idirstream(const char *spec) : h(-1), s(spec) {}
    explicit idirstream(const string &spec) : h(-1), s(spec) {}
 
    ~idirstream() { _findclose(h); }
    self &operator>>(file_data &x) { if (h<0) { h=_findfirst(s.c_str(),&x); if (h<0) setstate(failbit); } else if (_findnext(h,&x)<0) setstate(failbit); return *this; }
    void reset() { _findclose(h); h=-1; }
};
#else
class idirstream : public basic_ios<file_data>
{
  public:
    typedef idirstream self;

  private:
    int i;
    glob_t globbuf;
    struct stat daten;
    string s;

  public:
    idirstream() : i(-1),s("*") {}
    explicit idirstream(const char *spec) : i(-1),s(spec) {}
    explicit idirstream(const string &spec) : i(-1),s(spec) {}

    ~idirstream() { globfree(&globbuf); }
    self &operator>>(file_data &x)
    {
      if (i==-1) { if(glob(s.c_str(),0,0,&globbuf)==0) i=0; else { setstate(failbit); return *this; } }
      if (stat(globbuf.gl_pathv[i],&daten)!=0) { setstate(failbit); return *this; }

      strcpy(x.name,filename_out_of_pathname(globbuf.gl_pathv[i]).c_str());
      x.size=daten.st_size;
      x.time_access=daten.st_atime;
      x.time_create=daten.st_ctime;
      x.time_write=daten.st_mtime;
      x.attrib = _A_NORMAL;
      if (  daten.st_mode & S_IFDIR ) x.attrib |= _A_SUBDIR;
      if (!(daten.st_mode & S_IWUSR)) x.attrib |= _A_RDONLY;
      
      i++;
      return *this;
    }
    void reset() {  i=0; }
};
#endif


char       *name(file_data       &x) { return x.name; }
const char *name(const file_data &x) { return x.name; }

template<>
struct name_function<file_data> : public unary_reference_function<file_data,char *,char *, const char *>
{
  name_function() {}
  char       *operator()(argument_type       &x)       { return name(x); }
  const char *operator()(const argument_type &x) const { return name(x); }
};
template<>
struct name_function<const file_data> : public unary_reference_function<file_data,const char *,const char *, const char *>
{
  name_function() {}
  reference       operator()(argument_type       &x)       { return name(x); }
  const_reference operator()(const argument_type &x) const { return name(x); }
};

struct attributes_function : public unary_reference_function<file_data,unsigned int>
{
  attributes_function() {}
  reference       operator()(argument_type       &x)       { return x.attrib; }
  const_reference operator()(const argument_type &x) const { return x.attrib; }
};

//typedef unary_compose<binder2nd<bitwise_and<const unsigned int> >,attributes_function> has_attribute_pred;
//typedef unary_compose<logical_not<const unsigned int>, has_attribute_pred> has_not_attribute_pred;
//typedef filter_iterator<is_iterator<idirstream>,has_attribute_pred    > subdir_iterator;
//typedef filter_iterator<is_iterator<idirstream>,has_not_attribute_pred> file_iterator;
//
//file_iterator   file_it  ()               { return filter_it(is_iterator<idirstream>(  ), compose1(logical_not<const unsigned int>(),compose1(bind2nd(bitwise_and<const unsigned int>(),_A_SUBDIR),attributes_function()))); }
//file_iterator   file_it  (idirstream &is) { return filter_it(is_iterator<idirstream>(is), compose1(logical_not<const unsigned int>(),compose1(bind2nd(bitwise_and<const unsigned int>(),_A_SUBDIR),attributes_function()))); }
//subdir_iterator subdir_it()               { return filter_it(is_iterator<idirstream>(  ), compose1(bind2nd(bitwise_and<const unsigned int>(),_A_SUBDIR),attributes_function())); }
//subdir_iterator subdir_it(idirstream &is) { return filter_it(is_iterator<idirstream>(is), compose1(bind2nd(bitwise_and<const unsigned int>(),_A_SUBDIR),attributes_function())); }

typedef unary_compose<binder2nd<bitwise_and<unsigned int> >,attributes_function> has_attribute_pred;
typedef unary_compose<logical_not<unsigned int>, has_attribute_pred> has_not_attribute_pred;
typedef filter_iterator<is_iterator<idirstream>,has_attribute_pred    > subdir_iterator;
typedef filter_iterator<is_iterator<idirstream>,has_not_attribute_pred> file_iterator;

file_iterator   file_it  ()               { return filter_it(is_iterator<idirstream>(  ), compose1(logical_not<unsigned int>(),compose1(bind2nd(bitwise_and<unsigned int>(),_A_SUBDIR),attributes_function()))); }
file_iterator   file_it  (idirstream &is) { return filter_it(is_iterator<idirstream>(is), compose1(logical_not<unsigned int>(),compose1(bind2nd(bitwise_and<unsigned int>(),_A_SUBDIR),attributes_function()))); }
subdir_iterator subdir_it()               { return filter_it(is_iterator<idirstream>(  ), compose1(bind2nd(bitwise_and<unsigned int>(),_A_SUBDIR),attributes_function())); }
subdir_iterator subdir_it(idirstream &is) { return filter_it(is_iterator<idirstream>(is), compose1(bind2nd(bitwise_and<unsigned int>(),_A_SUBDIR),attributes_function())); }


template<class Ch,class Tr,class A> basic_string<Ch,Tr,A> prefix(const basic_string<Ch,Tr,A> &s) { return s.substr(0,s.rfind('.')); };
template<class Ch> basic_string<Ch> prefix(const Ch *s) { return prefix(basic_string<Ch>(s)); };

string prefix(const file_data &x) { return prefix(x.name); }


#endif

