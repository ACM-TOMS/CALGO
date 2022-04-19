///////////////////////////////////////////////////////////////////////////////
//
// File DisContField1D.h
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
// Description: Field definition in one-dimension for a discontinuous
// LDG-H expansion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD1D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD1D_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList0D.h>
#include <LocalRegions/PointExp.h>
#include <SpatialDomains/MeshGraph1D.h>
#include <SpatialDomains/Conditions.h>
#include <MultiRegions/GlobalLinSys.h>
//#include <MultiRegions/AssemblyMapDG.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>


namespace Nektar
{
    namespace MultiRegions
    {
        /// This class is the abstraction of a global discontinuous two-
        /// dimensional spectral/hp element expansion which approximates the
        /// solution of a set of partial differential equations.
        class DisContField1D: public ExpList1D
        {
        public:
            /// Default constructor.
            MULTI_REGIONS_EXPORT DisContField1D();

            /// Constructs a 1D discontinuous field based on a mesh and boundary
            /// conditions.
            MULTI_REGIONS_EXPORT DisContField1D(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                const std::string &variable);
            
            /// Constructor for a DisContField1D from a List of subdomains
            /// New Constructor for arterial network 
            MULTI_REGIONS_EXPORT DisContField1D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::CompositeMap& domain,
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                const std::string &variable,
                int i);

            /// Constructs a 1D discontinuous field based on an existing field.
            MULTI_REGIONS_EXPORT DisContField1D(const DisContField1D &In);
            
            /// Constructs a 1D discontinuous field based on an existing field.
	    /// (needed in order to use ContField( const ExpList1D &In) constructor
            MULTI_REGIONS_EXPORT DisContField1D(const ExpList1D &In);

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~DisContField1D();
            
            /// For a given key, returns the associated global linear system.
            MULTI_REGIONS_EXPORT GlobalLinSysSharedPtr GetGlobalBndLinSys(
                const GlobalLinSysKey &mkey);

        protected:
            /// The number of boundary segments on which Dirichlet boundary
            /// conditions are imposed.
            int m_numDirBndCondExpansions;

            /// Discretised boundary conditions.
            /**
             * It is an array of size equal to the number of boundary points
             * and consists of entries of the type LocalRegions#PointExp. Every
             * entry corresponds to a point on a single boundary region.
             */
            Array<OneD,MultiRegions::ExpListSharedPtr>         m_bndCondExpansions;

            /// An array which contains the information about the boundary
            /// condition on the different boundary regions.
            Array<OneD,SpatialDomains::BoundaryConditionShPtr> m_bndConditions;

            /// Global boundary matrix.
            GlobalLinSysMapShPtr                               m_globalBndMat;
            
            /// Trace space storage for points between elements.
            ExpListSharedPtr                                   m_trace;
            Array<OneD, ExpListSharedPtr>                      m_traces;
            Array<OneD, NekDouble>                             tmpBndSol;

            /// Local to global DG mapping for trace space.
            AssemblyMapDGSharedPtr                        m_traceMap;

            /// Discretises the boundary conditions.
            void GenerateBoundaryConditionExpansion(
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                SpatialDomains::BoundaryConditions &bcs,
                const std::string variable);
            
            // Discretises the boundary conditions in case of multidomain solver.
            void GenerateMultiDomainBoundaryConditionExpansion(
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                SpatialDomains::BoundaryConditions &bcs,
                const std::string variable,
                int subdomain);
            
            /// Generate a associative map of periodic vertices in a mesh.
            void GetPeriodicVertices(
                                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                                const SpatialDomains::BoundaryConditions &bcs,
                                const std::string variable,
                                      map<int,int>& periodicVertices);

            virtual ExpListSharedPtr &v_GetTrace()
            {
                return m_trace;
            }
            
            virtual ExpListSharedPtr &v_GetTrace(int i)
            {
                return m_traces[i];
            }
            
            virtual AssemblyMapDGSharedPtr &v_GetTraceMap(void)
            {
                return m_traceMap;
            }
            
            virtual void v_AddTraceIntegral(
                const Array<OneD, const NekDouble> &Fn,
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_GetFwdBwdTracePhys(
                      Array<OneD,       NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);
            virtual void v_GetFwdBwdTracePhys(
                const Array<OneD, const NekDouble> &field,
                      Array<OneD,       NekDouble> &Fwd,
                      Array<OneD,       NekDouble> &Bwd);
            virtual void v_ExtractTracePhys(
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_ExtractTracePhys(
                const Array<OneD, const NekDouble> &inarray, 
                      Array<OneD,       NekDouble> &outarray);

            /// Populates the list of boundary condition expansions.
            void SetBoundaryConditionExpansion(
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                const SpatialDomains::BoundaryConditions &bcs,
                const std::string variable,
                Array<OneD, MultiRegions::ExpListSharedPtr>
                    &bndCondExpansions,
                Array<OneD, SpatialDomains
                    ::BoundaryConditionShPtr> &bndConditions);
            
            /// Populates the list of boundary condition expansions in multidomain case.
            void SetMultiDomainBoundaryConditionExpansion(
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                const SpatialDomains::BoundaryConditions &bcs,
                const std::string variable,
                Array<OneD, MultiRegions::ExpListSharedPtr>
                    &bndCondExpansions,
                Array<OneD, SpatialDomains
                    ::BoundaryConditionShPtr> &bndConditions,
                int subdomain);
            
            void GenerateFieldBnd1D(
                SpatialDomains::BoundaryConditions &bcs,
                const std::string variable);
            
            virtual map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo();
            
            virtual const Array<OneD,const MultiRegions::ExpListSharedPtr>
                &v_GetBndCondExpansions()
            {
                return m_bndCondExpansions;
            }

            virtual const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>
                &v_GetBndConditions()
            {
                return m_bndConditions;
            }

            virtual MultiRegions::ExpListSharedPtr &v_UpdateBndCondExpansion(int i)
            {
                return m_bndCondExpansions[i];
            }

            virtual Array<OneD, SpatialDomains::BoundaryConditionShPtr> &v_UpdateBndConditions()
            {
                return m_bndConditions;
            }

            virtual void v_GetBoundaryToElmtMap(
                Array<OneD,int> &ElmtID, Array<OneD,int> &VertID);
			
            /// Evaluate all boundary conditions at a given time..
            virtual void v_EvaluateBoundaryConditions(
                const NekDouble time = 0.0,
                const NekDouble x2_in = NekConstants::kNekUnsetDouble,
                const NekDouble x3_in = NekConstants::kNekUnsetDouble);
            
            /// Solve the Helmholtz equation.
            virtual void v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const FlagList &flags,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeff,
                    const Array<OneD, const NekDouble> &dirForcing);
        };

        typedef boost::shared_ptr<DisContField1D>   DisContField1DSharedPtr;
    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD1D_H
