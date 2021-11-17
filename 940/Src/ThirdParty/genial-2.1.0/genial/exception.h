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

#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <string>
#include <exception>

using namespace std;

class error : public exception
{
  protected:
    string msg;

  public:
    error() throw() : exception() {}
    explicit error(const string &s) throw() : msg(s) {}
    explicit error(const char *s) throw() : msg(s) {}

    error(const error &r) : msg(r.msg) {}

    virtual ~error() throw() {}

    virtual const char *what() const throw() { return msg.c_str(); } 
    
    string       &message()       { return msg; }
    const string &message() const { return msg; }
};

#endif
