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

#ifndef PPM_H
#define PPM_H

#include <sstream>

#include "bitstream.h"
#include "imagefile.h"


//namespace genial
//{


template<class Stream, class Image>
Stream &ppm_read(Stream &is, Image &X)
{
  typedef Image image_type;
  typedef typename image_type::int_type int_type;

  string s;
  getline(is, s);
  if (s!="P6") throw error("Error: Load PPM (incorrect format)");

  do { getline(is,s); } while (s.size()==0 || s[0]=='#');
  int_type M,N;
  istringstream(s.c_str()) >> N >> M;
  if (M*N==0) throw error("read ppm: null size");
  X.resize(M,N);

  do { getline(is,s); } while (s.size()==0 || s[0]=='#');
  unsigned int pixmax;
  istringstream(s.c_str()) >> pixmax;
  if (pixmax>255) throw error("read ppm: bad format");

  for (typename Image::iterator it=X.begin(),itmax=X.end(); it!=itmax; ++it)
    { BYTE R,G,B; is > R > G > B; *it=RGBA(R,G,B); }

  return is;
}

template<class Stream, class Image>
Stream &ppm_write(Stream &os, const Image &X)
{
  typedef Image image_type;
  typedef typename image_type::int_type int_type;

  os << "P6" << endl;
  os << X.ncols() << " " << X.nrows() << endl;
  os << 255 << endl;

  for (typename Image::const_iterator it=X.begin(),itmax=X.end(); it!=itmax; ++it)
    { RGBA c=*it; os < (BYTE)c.x < (BYTE)c.y < (BYTE)c.z; }

  return os;
}

//{unsecret}
//{group:Images File Formats}
//Summary: PPM format
class PPMFile : public ImageFile
{
  public:
    PPMFile(const string &s) : ImageFile(s) {}   
};

//{unsecret}
template<class G> void operator>>(const PPMFile &file, Matrix<G> &X)
{
  ifstream is(file.name().c_str(), ios::binary|ios::in);
  if (!is.is_open()) throw error("Error: Load PPM (" + file.name() + " could not be opened)");
  ppm_read(is, X);  
}

//{unsecret}
template<class G> void operator<<(const PPMFile &file, const Matrix<G> &X)
{
  ofstream os(file.name().c_str(), ios::binary|ios::out);
  if (!os.is_open()) throw error("Error: Save PPM (" + file.name() + " could not be opened)");
  ppm_write(os, X);  
}

//}

#endif

