///////////////////////////////////////////////////////////////////////////////
//
// File $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/LocalRegions/PointExp.h,v $
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
// Description: Definition of a Point expansion 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef POINTEXP_H
#define POINTEXP_H

#include <LocalRegions/LocalRegions.hpp>

#include <StdRegions/StdPointExp.h>

#include <SpatialDomains/MeshComponents.h>
#include <LocalRegions/LocalRegionsDeclspec.h>
#include <LocalRegions/Expansion0D.h>
#include <StdRegions/StdExpansion.h>

namespace Nektar
{
    namespace LocalRegions
    {
        class PointExp: virtual public StdRegions::StdPointExp, virtual public Expansion0D
        {
        public:
            LOCAL_REGIONS_EXPORT PointExp(const SpatialDomains::VertexComponentSharedPtr &m_geom);
            LOCAL_REGIONS_EXPORT ~PointExp(void);

            inline const Array<OneD, const NekDouble>& GetCoeffs(void) const
            {
                return m_coeffs;
            }

            inline NekDouble  GetCoeffs(int i) const
            {
                ASSERTL1(i == 0,"index out of range");

                return m_coeffs[i];
            }
		
            inline NekDouble  GetPhys(int i) const
            {
                ASSERTL1(i == 0,"index out of range");
		
                return m_phys[i];
            }
            
            inline NekDouble  GetCoeff(int i) const
            {
                ASSERTL1(i == 0,"index out of range");

                return m_coeffs[i];
            }

            inline Array<OneD, NekDouble>& UpdateCoeffs(void)
            {
                return(m_coeffs);
            }

            inline void  SetCoeff(const NekDouble value)
            {
                m_coeffs[0] = value;
            }

            inline const Array<OneD, const NekDouble>& GetPhys(void) const
            {
                return m_phys;
            }
  
            inline Array<OneD, NekDouble>& UpdatePhys(void) 
            {
                return(m_phys);
            }

            inline void  SetPhys(const NekDouble value)
            {
                m_phys[0] = value;
            }

            inline void GetCoords(NekDouble &x, NekDouble &y, NekDouble &z)
            {
                m_geom->GetCoords(x,y,z);
            }

            inline void GetCoords(Array<OneD,NekDouble> &coords)
            {
                m_geom->GetCoords(coords);
            }
            
            inline const SpatialDomains::VertexComponentSharedPtr &GetGeom(void) const
            {
                return m_geom;
            }

            inline const SpatialDomains::VertexComponentSharedPtr &GetVertex(void) const
            {
                return m_geom;
            }
            
        protected:
            Array<OneD, NekDouble > m_coeffs; //!< Array containing expansion coefficients
            Array<OneD, NekDouble > m_phys; //!< Array containing physical point which is likely to be the same as the coefficient but is defined for consistency (It is also used in Robin boundary conditions) 
            SpatialDomains::VertexComponentSharedPtr m_geom;
            
            const SpatialDomains::GeometrySharedPtr v_GetGeom() const
            {
                return m_geom;
            }
        };
        
        // type defines for use of PointExp in a boost vector
        typedef boost::shared_ptr<PointExp> PointExpSharedPtr;
        typedef std::vector<PointExpSharedPtr> PointExpVector;
        typedef std::vector<PointExpSharedPtr>::iterator PointExpVectorIter;
        
        const static Array<OneD, PointExpSharedPtr> NullPointExpSharedPtrArray;
    } //end of namespace
} //end of namespace

#endif // POINTEXP_H
