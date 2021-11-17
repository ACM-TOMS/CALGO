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

#ifndef IOIMAGE_H
#define IOIMAGE_H

#include "image/image.h"
#include "image/imagefile.h"
#include "gdi.h"


//namespace genial
//{

template<class G>
void operator>>(const ImageFile &file, Matrix<G> &X)
{
  //Get File Extension
  string ext=file.name().substr(file.name().rfind('.')+1);

#ifdef PPM_H   
  if(ext=="ppm")  { PPMFile(file.name()) >> X;return;  }
#endif
#ifdef PGM_H
  if(ext=="pgm")  { PGMFile(file.name()) >> X;return;  }
#endif
#ifdef PBM_H 
  if(ext=="pbm")  { PBMFile(file.name()) >> X;return;  }
#endif
#ifdef BMP_H
  if(ext=="bmp")  { BMPFile(file.name()) >> X;return;  }
#endif
#ifdef GIF_H
  if(ext=="gif")  { GIFFile(file.name()) >> X;return;  }
#endif
#ifdef PNG_H
  if(ext=="png")  { PNGFile(file.name()) >> X;return;  }
#endif
#ifdef JPG_H 
  if(ext=="jpg" || ext=="jpeg") { JPGFile(file.name()) >> X;return;  }
#endif
#ifdef TIF_H
  if(ext=="tif" || ext=="tiff") { TIFFile(file.name()) >> X;return;  }
#endif
  
  throw error("Read File: Format \"" + ext + "\" unknown");
}

template<class G>
void operator<<(const ImageFile &file, const Matrix<G> &X)
{
  //Get File Extension
  string ext=file.name().substr(file.name().rfind('.')+1);

#ifdef PPM_H   
  if(ext=="ppm")  { PPMFile(file.name()) << X; return; }
#endif
#ifdef PGM_H
  if(ext=="pgm")  { PGMFile(file.name()) << X; return; }
#endif
#ifdef PBM_H 
  if(ext=="pbm")  { PBMFile(file.name()) << X; return; }
#endif
#ifdef BMP_H
  if(ext=="bmp")  { BMPFile(file.name()) << X; return; }
#endif
#ifdef GIF_H
  if(ext=="gif")  { GIFFile(file.name()) << X; return; }
#endif
#ifdef PNG_H
  if(ext=="png")  { PNGFile(file.name()) << X; return; }
#endif
#ifdef JPG_H 
  if(ext=="jpg" || ext=="jpeg") { JPGFile(file.name()) << X; return; }
#endif
#ifdef TIF_H
  if(ext=="tif" || ext=="tiff") { TIFFile(file.name()) << X; return; }
#endif

  throw error("Write File: Format \"" + ext + "\" unknown");
}


//}

#endif