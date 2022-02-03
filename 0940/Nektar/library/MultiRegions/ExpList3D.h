///////////////////////////////////////////////////////////////////////////////
//
// File ExpList3D.h
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
// Description: Expansion list 3D header definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPLIST3D_H
#define EXPLIST3D_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <vector>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList2D.h>

#include <LocalRegions/HexExp.h>
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/PyrExp.h>
#include <LocalRegions/TetExp.h>
#include <SpatialDomains/MeshGraph3D.h>
#include <SpatialDomains/Conditions.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /// Abstraction of a three-dimensional multi-elemental expansion which
        /// is merely a collection of local expansions.
        class ExpList3D:  public ExpList
        {
        public:
            /// Default constructor.
            MULTI_REGIONS_EXPORT ExpList3D();

            /// Copy constructor.
            MULTI_REGIONS_EXPORT ExpList3D(const ExpList3D &In);

            MULTI_REGIONS_EXPORT ExpList3D(  
                        const LibUtilities::SessionReaderSharedPtr &pSession,
                        const LibUtilities::BasisKey &TBa,
                        const LibUtilities::BasisKey &TBb,
                        const LibUtilities::BasisKey &TBc,
                        const LibUtilities::BasisKey &HBa,
                        const LibUtilities::BasisKey &HBb,
                        const LibUtilities::BasisKey &HBc,
                        const SpatialDomains::MeshGraphSharedPtr &graph3D,
                        const LibUtilities::PointsType TetNb
                                            = LibUtilities::SIZE_PointsType);

            /// Sets up a list of local expansions based on an input mesh.
            MULTI_REGIONS_EXPORT ExpList3D(
                        const LibUtilities::SessionReaderSharedPtr &pSession,
                        const SpatialDomains::MeshGraphSharedPtr &graph3D);

            /// Sets up a list of local expansions based on an expansion vector
            ExpList3D(const SpatialDomains::ExpansionMap &expansions);

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~ExpList3D();

        protected:

            /// Set up the normals on each expansion.
            virtual void v_SetUpPhysNormals();

        private:
            /// Definition of the total number of degrees of freedom and
            /// quadrature points. Sets up the storage for \a m_coeff and \a
            ///  m_phys.
            void SetCoeffPhys(void);

/*            LocalRegions::HexExpVector    m_hex;
            LocalRegions::PrismExpVector  m_prism;
            LocalRegions::PyrExpVector    m_pyr;
            LocalRegions::TetExpVector    m_tet;
*/

            virtual void v_ReadGlobalOptimizationParameters();

            virtual void v_WriteVtkPieceHeader(std::ofstream &outfile, int expansion);

        };

        /// Shared pointer to an ExpList3D object.
        typedef boost::shared_ptr<ExpList3D>      ExpList3DSharedPtr;
        /// Vector of pointers to ExpList3D objects.
        typedef std::vector<ExpList3DSharedPtr>   ExpList3DVector;
        /// Iterator over an ExpList3DVector.
        typedef std::vector<ExpList3DSharedPtr>::iterator ExpList3DVectorIter;

    } //end of namespace
} //end of namespace

#endif
