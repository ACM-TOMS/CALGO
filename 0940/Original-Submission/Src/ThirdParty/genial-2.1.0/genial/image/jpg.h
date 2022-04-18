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

#ifndef JPG_H
#define JPG_H

#include "imagefile.h"
#include "gdi.h"


//{unsecret}
//{group:Images File Formats}
//Summary: JPG format
class JPGFile : public ImageFile
{
  public:
    JPGFile(const string &s,const unsigned long qual=75) : ImageFile(s) { quality=qual; }
    unsigned long quality;
};

//{unsecret}
template<class G> void operator>>(const JPGFile &file, Matrix<G> &X)
{
  typedef Matrix<G> image_type;
  typedef typename image_type::value_type value_type;
  typedef typename image_type::int_type int_type;
  gdi::Bitmap Y(wstring(file.name().begin(),file.name().end()).c_str());
  if (Y.GetHeight()*Y.GetWidth()==0) throw error("Error: Load JPG");
  X.resize(Y.GetHeight(),Y.GetWidth());
  gdi::Color color(0,0,0,0);
  for(int_type i=0,imax=X.nrows();i<imax;i++)
  {
    typename MatrixRow<image_type>::self Xi=row(X,i);
    for(int_type j=0,jmax=X.ncols();j<jmax;j++)
    {
      Y.GetPixel(j,i,&color);
      Xi[j]=value_type(RGBA(color.GetRed(),color.GetGreen(),color.GetBlue()));
    }
  }
}

//{unsecret}
template<class G> void operator<<(const JPGFile &file, const Matrix<G> &X)
{
  typedef const Matrix<G> image_type;
  typedef typename image_type::int_type  int_type;
  gdi::Bitmap Y(X.ncols(),X.nrows());
  gdi::Color color(0,0,0,0);
  for(int_type i=0,imax=X.nrows();i<imax;i++)
  {
    typename MatrixRow<image_type>::self Xi=row(X,i);
    for(int_type j=0,jmax=X.ncols();j<jmax;j++)
    {
      RGBA c(Xi[j]);
      color.SetValue(gdi::Color::MakeARGB(255,c.x,c.y,c.z));
      Y.SetPixel(j,i,color);
    }
  }
  CLSID Clsid;
  GetEncoderClsid(L"image/jpeg",&Clsid);
  gdi::EncoderParameters encoderParameters;
  encoderParameters.Count = 1;
  encoderParameters.Parameter[0].Guid = gdi::EncoderQuality;
  encoderParameters.Parameter[0].Type = gdi::EncoderParameterValueTypeLong;
  encoderParameters.Parameter[0].NumberOfValues = 1;
  encoderParameters.Parameter[0].Value = (unsigned long *)&file.quality;
  check_gdi(Y.Save(wstring(file.name().begin(),file.name().end()).c_str(),&Clsid,&encoderParameters),"Error: Save JPG");
}

#endif