///////////////////////////////////////////////////////////////////////////////
//
// File: ConsistentObjectAccess.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: A wrapper around objects or pointers to objects which give
// the same interface for both types (i.e., -> will work for both).
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_CONSISTENT_ACCESS_OBJECT_HPP
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_CONSISTENT_ACCESS_OBJECT_HPP

#include <boost/shared_ptr.hpp>
#include <boost/call_traits.hpp>
#include <boost/type_traits.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

#include <iostream>

namespace Nektar
{    
    template<typename DataType>
    struct ConsistentObjectAccess
    {     
        static DataType& reference(DataType& o) { return o; }
        static const DataType& const_reference(const DataType& o) { return o; }
        static DataType* pointer(DataType& o) { return &o; }
        static const DataType* const_pointer(const DataType& o) { return &o; }
        
        static bool ReferencesObject(const DataType& o) { return true; }
    };
    
    template<typename DataType>
    struct ConsistentObjectAccess<DataType*> 
    {
        static const DataType& const_reference(DataType* o) { ASSERTL1(o != 0, "Can't dereference null pointer."); return *o; }
        static const DataType* const_pointer(DataType* o) { return o; }
        static bool ReferencesObject(DataType* o) { return o != 0; }

        static DataType& reference(DataType* o) { ASSERTL1(o != 0, "Can't dereference null pointer."); return *o; }
        static DataType* pointer(DataType* o) { return o; }
    };    
    

    template<typename DataType>
    struct ConsistentObjectAccess<boost::shared_ptr<DataType> >
    {
        static const DataType& const_reference(const boost::shared_ptr<DataType>& o) { ASSERTL1(o, "Can't dereference null pointer."); return *o; }
        static const DataType* const_pointer(const boost::shared_ptr<DataType>& o) { return o.get(); }
        static DataType& reference(const boost::shared_ptr<DataType>& o) { ASSERTL1(o, "Can't dereference null pointer."); return *o; }
        static DataType* pointer(const boost::shared_ptr<DataType>& o) { return o.get(); }
        static bool ReferencesObject(const boost::shared_ptr<DataType>& o) { return o; }
    };
}
    
#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_CONSISTENT_ACCESS_OBJECT_HPP

/**
    $Log: ConsistentObjectAccess.hpp,v $
    Revision 1.4  2007/07/22 23:03:25  bnelson
    Backed out Nektar::ptr.

    Revision 1.3  2007/07/20 00:39:36  bnelson
    Replaced boost::shared_ptr with Nektar::ptr

    Revision 1.2  2007/01/29 01:35:17  bnelson
    Removed memory manager requirements.

    Revision 1.1  2006/11/06 17:06:20  bnelson
    *** empty log message ***

 **/
 
