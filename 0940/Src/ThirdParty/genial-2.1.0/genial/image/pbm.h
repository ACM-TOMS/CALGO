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

#ifndef PBM_H
#define PBM_H

#include <sstream>

#include "bitstream.h"
#include "imagefile.h"

//namespace genial
//{


template<class Stream, class Image>
Stream &pbm_read(Stream &is, Image &X)
{
  typedef Image image_type;
  typedef typename image_type::int_type int_type;

  string s;
  getline(is, s);
  if (s!="P4") throw error("error: Load PBM (incorrect format)");

  do { getline(is,s); } while (s.size()==0 || s[0]=='#');
  int_type M,N;
  istringstream(s.c_str()) >> N >> M;
  if (M*N==0) throw error("read pbm : null size");
  X.resize(M,N);

  BYTE c;
  for (int_type i=0,imax=X.nrows(); i<imax; ++i)
  {
    typename MatrixRow<typename Image<V>::base>::self Xi=row(X,i);
    for (int_type j=0,jmax=X.ncols(); j<jmax; ++j)
    {
      if (!(j&7)) is>c;
      Xi[j] = (c>>(7-j&7))&1 ? (unsigned char)0 : (unsigned char)-1;
    }
  }

  return is;
}

template<class Stream, class Image>
Stream &pbm_write(Stream &os, const Image &X)
{
  typedef Image image_type;
  typedef typename image_type::int_type int_type;
  typedef typename image_type::index_type index_type;

  os << "P4" << endl;
  os << X.ncols() << " " << X.nrows() << endl;

  BYTE c=0;
  for (int_type i=0,imax=X.nrows(); i<imax; ++i)
  {
    typename MatrixRow<typename Image<V>::base>::self Xi=row(X,i);
    for (int_type j=0,jmax=X.ncols(); j<jmax; ++j)
    {
      if (!Xi[j]) c |= 1<<(7-j&7);
      if ((j&7)==7 || X.rightboundary(index_type(i,j))) { os < c; c=0; }
    }
  }

  return os;
}


//{unsecret}
//{group:Images File Formats}
//Summary: PBM format
class PBMFile : public ImageFile
{
  public:
    PBMFile(const char *s) : ImageFile(s) {}   
};

//{unsecret}
template<class G> void operator>>(const PBMFile &file, Matrix<G> &X)
{
  ifstream is(file.name().c_str(), ios::binary|ios::in);
  if (!is.is_open()) throw error("Error: Load PBM (" + file.name() + " could not be opened)");
  pbm_read(is, X);  
}

//{unsecret}
template<class G> void operator<<(const PBMFile &file, const Matrix<G> &X)
{
  ofstream os(file.name().c_str(), ios::binary|ios::out);
  if (!os.is_open()) throw error("Error: Save PBM (" + file.name() + " could not be opened)");
  pbm_write(os, X);  
}


//}

#endif

