///////////////////////////////////////////////////////////////////////////////
//
// File GlobalMatrixKey.cpp
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
// Description: Definition of GlobalMatrixKey
//
///////////////////////////////////////////////////////////////////////////////


#include <MultiRegions/GlobalMatrixKey.h>

namespace Nektar
{
    namespace MultiRegions
    {
        GlobalMatrixKey::GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                    const AssemblyMapSharedPtr &locToGloMap,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varCoeffs) :
            m_matrixType(matrixType),
            m_expansionType(StdRegions::eNoExpansionType),
            m_locToGloMap(locToGloMap),
            m_constFactors(factors),
            m_varCoeffs(varCoeffs)
        {
        }

        GlobalMatrixKey::GlobalMatrixKey(const GlobalMatrixKey &key,
                        const StdRegions::ExpansionType expType):
            m_matrixType(key.m_matrixType),
            m_expansionType(expType),
            m_locToGloMap(key.m_locToGloMap),
            m_constFactors(key.m_constFactors),
            m_varCoeffs(key.m_varCoeffs)
        {
        }

        GlobalMatrixKey::GlobalMatrixKey(const GlobalMatrixKey &key):
            m_matrixType(key.m_matrixType),
            m_expansionType(key.m_expansionType),
            m_locToGloMap(key.m_locToGloMap),
            m_constFactors(key.m_constFactors),
            m_varCoeffs(key.m_varCoeffs)
        {
        }

        GlobalMatrixKey::~GlobalMatrixKey()
        {
        }

        bool operator<(const GlobalMatrixKey &lhs, const GlobalMatrixKey &rhs)
        {
            if(lhs.m_matrixType < rhs.m_matrixType)
            {
                return true;
            }

            if(lhs.m_matrixType > rhs.m_matrixType)
            {
                return false;
            }


            if(lhs.m_expansionType < rhs.m_expansionType)
            {
                return true;
            }
            

            if(lhs.m_expansionType > rhs.m_expansionType)
            {
                return false;
            }
            
            if(lhs.m_constFactors.size() < rhs.m_constFactors.size())
            {
                return true;
            }
            else if(lhs.m_constFactors.size() > rhs.m_constFactors.size())
            {
                return false;
            }
            else
            {
                StdRegions::ConstFactorMap::const_iterator x, y;
                for(x = lhs.m_constFactors.begin(), y = rhs.m_constFactors.begin();
                    x != lhs.m_constFactors.end(); ++x, ++y)
                {
                    if (x->second < y->second)
                    {
                        return true;
                    }
                    if (x->second > y->second)
                    {
                        return false;
                    }
                }
            }

            if(lhs.m_varCoeffs.size() < rhs.m_varCoeffs.size())
            {
                return true;
            }
            else if(lhs.m_varCoeffs.size() > rhs.m_varCoeffs.size())
            {
                return false;
            }
//            else
//            {
//                StdRegions::VarCoeffMap::const_iterator x, y;
//                for (x = lhs.m_varCoeffs.begin(), y = rhs.m_varCoeffs.begin();
//                     x != lhs.m_varCoeffs.end(); ++x, ++y)
//                {
//                    if (x->second.get() < y->second.get())
//                    {
//                        return true;
//                    }
//                    if (x->second.get() > y->second.get())
//                    {
//                        return false;
//                    }
//                }
//            }

            if(!rhs.m_locToGloMap.get())
            {
                return false;
            }
            else if(!lhs.m_locToGloMap.get() && rhs.m_locToGloMap.get() )
            {
                return true;
            }
            if(lhs.m_locToGloMap->GetHash() < rhs.m_locToGloMap->GetHash())
            {
                return true;
            }

            return false;
        }

        std::ostream& operator<<(std::ostream& os, const GlobalMatrixKey& rhs)
        {
            os << "MatrixType: " << rhs.GetMatrixType() << endl;
            os << "Number of constants: " << rhs.GetNConstFactors() << endl;
            StdRegions::ConstFactorMap::const_iterator x;
            for(x = rhs.GetConstFactors().begin(); x != rhs.GetConstFactors().end(); ++x)
            {
                os << "  Constant " << StdRegions::ConstFactorTypeMap[x->first]
                   << ": " << x->second << endl;
            }
            os << "Number of variable coefficients: " 
               << rhs.GetNVarCoeffs() << endl;

            return os;
        }
    }
}
