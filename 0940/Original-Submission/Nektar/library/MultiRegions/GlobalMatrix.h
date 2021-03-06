///////////////////////////////////////////////////////////////////////////////
//
// File GlobalMatrix.h
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
// Description: GlobalMatrix header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALMATRIX_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALMATRIX_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/GlobalMatrixKey.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /// Represents a matrix of all degrees of freedom.
        class GlobalMatrix
        {
        public:
            typedef std::pair<  int,  int> CoordType;
            typedef std::map< CoordType, NekDouble > COOMatType;

            /// Construct a new matrix.
            MULTI_REGIONS_EXPORT GlobalMatrix(unsigned int rows,
                         unsigned int columns,
                         const COOMatType &cooMat);

            MULTI_REGIONS_EXPORT ~GlobalMatrix() {}

            /// Returns a pointer to the DNekSparseMat matrix.
            const DNekSparseMatSharedPtr GetMatrix(void) const
            {
                return m_matrix;
            }

            /// Perform a matrix-vector multiply.
            MULTI_REGIONS_EXPORT void Multiply(const Array<OneD,const NekDouble> &in,
                                Array<OneD,      NekDouble> &out);

        private:
            /// Pointer to a double-precision Nektar++ sparse matrix.
            DNekSparseMatSharedPtr  m_matrix;
        };

        /// Shared pointer to a GlobalMatrix object.
        typedef boost::shared_ptr<GlobalMatrix> GlobalMatrixSharedPtr;
        /// Mapping from global matrix keys to global matrices.
        typedef map<GlobalMatrixKey,GlobalMatrixSharedPtr> GlobalMatrixMap;
        /// Shared pointer to a global matrix map.
        typedef boost::shared_ptr<GlobalMatrixMap> GlobalMatrixMapShPtr;

    } //end of namespace
} //end of namespace

#endif
