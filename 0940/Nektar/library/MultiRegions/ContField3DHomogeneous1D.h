///////////////////////////////////////////////////////////////////////////////
//
// File ContField3DHomogeneous1D.h
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
// Description: Field definition in three-dimensions for a continuous
// expansion with a homogeneous direction in 1D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_CONTFIELD3DHOMO1D_H
#define NEKTAR_LIBS_MULTIREGIONS_CONTFIELD3DHOMO1D_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <LibUtilities/Communication/Comm.h>

#include <MultiRegions/DisContField3DHomogeneous1D.h>
#include <MultiRegions/ContField2D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class ContField3DHomogeneous1D: public DisContField3DHomogeneous1D
        {
        public:
            MULTI_REGIONS_EXPORT ContField3DHomogeneous1D();

            MULTI_REGIONS_EXPORT ContField3DHomogeneous1D(
                           const LibUtilities::SessionReaderSharedPtr &pSession,
                           const LibUtilities::BasisKey &HomoBasis,
                           const NekDouble lhom,
						   const bool useFFT,
						   const bool dealiasing,
                           const SpatialDomains::MeshGraphSharedPtr &graph2D,
                           const std::string &variable,
						   const bool CheckIfSingularSystem = false);
            
            /// Copy constructor.
            MULTI_REGIONS_EXPORT ContField3DHomogeneous1D(const ContField3DHomogeneous1D &In);

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~ContField3DHomogeneous1D();
            
        protected:

        private:
            /// Template method virtual forwarded for LocalToGlobal()
            virtual void v_LocalToGlobal();

            /// Template method virtual forwarded for GlobalToLocal()
            virtual void v_GlobalToLocal();

            /// Solves the three-dimensional Helmholtz equation, subject to the
            /// boundary conditions specified.
            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const FlagList &flags,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeff,
                    const Array<OneD, const NekDouble> &dirForcing);
        };

        typedef boost::shared_ptr<ContField3DHomogeneous1D>  
            ContField3DHomogeneous1DSharedPtr;

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_CONTFIELD3DHOMO1D_H
