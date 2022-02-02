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

#ifndef IMAGEFILE_H
#define IMAGEFILE_H


#include <fstream>
#include <sstream>

#include "image/image.h"
#include "image/bitstream.h"

//namespace genial
//{

class ImageFile : public File
{
  public:
    ImageFile(const char   *s) : File(s) {}   
    ImageFile(const string &s) : File(s) {}   
};

//}

#endif








