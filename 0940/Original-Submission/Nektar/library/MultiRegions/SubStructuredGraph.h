///////////////////////////////////////////////////////////////////////////////
//
// File SubStructuredGraph.h
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
// Description: a collection of classes that facilitates the implementation
//              of the multi-level static condensation routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_SUBSTRUCTUREDGRAPH_H
#define MULTIREGIONS_SUBSTRUCTUREDGRAPH_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

namespace Nektar
{
    namespace MultiRegions
    {

        class BottomUpSubStructuredGraph;
        class SubGraph;
        class MultiLevelBisectedGraph;
        class PatchMap;

        typedef boost::shared_ptr<BottomUpSubStructuredGraph> BottomUpSubStructuredGraphSharedPtr;
        typedef boost::shared_ptr<SubGraph>                   SubGraphSharedPtr;
        typedef boost::shared_ptr<MultiLevelBisectedGraph>    MultiLevelBisectedGraphSharedPtr;
        typedef boost::shared_ptr<PatchMap>                   PatchMapSharedPtr;


        class PatchMap
        {
        public:

            MULTI_REGIONS_EXPORT  PatchMap(void);

            MULTI_REGIONS_EXPORT  PatchMap(const int vals);

            MULTI_REGIONS_EXPORT  ~PatchMap(void);

#if 0 
            MULTI_REGIONS_EXPORT void SetPatchMap(const int n, const int patchId, const int dofId,const bool bndPatch,const NekDouble sign);
#else
            MULTI_REGIONS_EXPORT void SetPatchMap(const int n, const int patchId, const int dofId,const unsigned int bndPatch,const NekDouble sign);
#endif
            inline Array<OneD, const int> GetPatchId() const 
            {
                return m_patchId;
            }

            inline Array<OneD, const int>  GetDofId() const
            {
                return m_dofId;
            }
#if 0 
            inline Array<OneD, const bool> IsBndDof() const
            {
                return m_bndPatch;
            }
#else
            inline Array<OneD, const unsigned int> IsBndDof() const
            {
                return m_bndPatch;
            }
#endif
            inline Array<OneD, const NekDouble> GetSign() const
            {
                return m_sign;
            }
            
        protected:
            Array<OneD, int > m_patchId;
            Array<OneD, int > m_dofId;
            //            Array<OneD, bool> m_bndPatch; 
            Array<OneD, unsigned int> m_bndPatch; 
            Array<OneD, NekDouble> m_sign; 
        };


        class SubGraph
        {
        public:
            
            MULTI_REGIONS_EXPORT SubGraph(const int nVerts, const int idOffset = 0):
            m_nVerts(nVerts),
                m_idOffset(idOffset)
                {
                }
            
            MULTI_REGIONS_EXPORT ~SubGraph(void)
            {
            }
            
            inline int GetNverts(void) const
            {
                return m_nVerts;
            }

            inline void SetNverts(const int i) 
            {
                m_nVerts = i;
            }

            inline int GetIdOffset(void) const
            {
                return m_idOffset;
            }

            inline void SetIdOffset(const int i) 
            {
                m_idOffset = i;
            }

        protected:
            int m_nVerts;
            int m_idOffset;
        };

        bool SubGraphWithoutVerts(const SubGraphSharedPtr g);

        class MultiLevelBisectedGraph
        {
        public:
            MULTI_REGIONS_EXPORT MultiLevelBisectedGraph(const Array<OneD, const int> sepTree);
            MULTI_REGIONS_EXPORT MultiLevelBisectedGraph(const int nBndDofs);

            MULTI_REGIONS_EXPORT ~MultiLevelBisectedGraph(void);

            MULTI_REGIONS_EXPORT int  GetTotDofs() const;

            MULTI_REGIONS_EXPORT void SetGlobalNumberingOffset();

            MULTI_REGIONS_EXPORT void DumpNBndDofs(void) const;

            MULTI_REGIONS_EXPORT void CollectLeaves(vector<SubGraphSharedPtr>& leaves) const;

            MULTI_REGIONS_EXPORT inline int  GetNdaughterGraphs() const;

            MULTI_REGIONS_EXPORT int CutLeaves();

            MULTI_REGIONS_EXPORT int CutEmptyLeaves();

            inline const SubGraphSharedPtr GetBndDofsGraph() const
            {
                return m_BndDofs;
            }

        protected:
            SubGraphSharedPtr m_BndDofs;
            MultiLevelBisectedGraphSharedPtr m_leftDaughterGraph;
            MultiLevelBisectedGraphSharedPtr m_rightDaughterGraph;        
        };


        class BottomUpSubStructuredGraph
        {
        public:
            MULTI_REGIONS_EXPORT BottomUpSubStructuredGraph(const Array<OneD, const int> septree);
            MULTI_REGIONS_EXPORT BottomUpSubStructuredGraph(const MultiLevelBisectedGraphSharedPtr& graph);
            MULTI_REGIONS_EXPORT BottomUpSubStructuredGraph(const int nVerts);

            MULTI_REGIONS_EXPORT ~BottomUpSubStructuredGraph(void);

            MULTI_REGIONS_EXPORT int GetTotDofs() const;

            MULTI_REGIONS_EXPORT void UpdateBottomUpReordering(Array<OneD,       int>& perm,  Array<OneD,  int>& iperm) const;

            MULTI_REGIONS_EXPORT void ExpandGraphWithVertexWeights(const Array<OneD, const int>& wgts);

            MULTI_REGIONS_EXPORT void MaskPatches(const int leveltomask, Array<OneD, NekDouble>& maskarray) const;
            
            MULTI_REGIONS_EXPORT int GetNpatchesWithInterior(const int whichlevel) const;

            MULTI_REGIONS_EXPORT void GetNintDofsPerPatch(const int whichlevel, Array<OneD, unsigned int>& outarray) const;
            
            MULTI_REGIONS_EXPORT int GetInteriorOffset(const int whichlevel, const int patch = 0) const;

            MULTI_REGIONS_EXPORT vector<SubGraphSharedPtr> GetInteriorBlocks(const int whichlevel) const;

            MULTI_REGIONS_EXPORT int GetNumGlobalDofs(const int whichlevel) const;

            MULTI_REGIONS_EXPORT int GetNlevels() const;

            MULTI_REGIONS_EXPORT void Dump() const;

        protected:
            vector<SubGraphSharedPtr> m_IntBlocks;
            BottomUpSubStructuredGraphSharedPtr m_daughterGraph;

        private:
            void SetBottomUpReordering(Array<OneD, int>& iperm) const;

            inline BottomUpSubStructuredGraphSharedPtr GetDaughterGraph() const
            {
                return m_daughterGraph;
            }
            inline vector<SubGraphSharedPtr> GetInteriorBlocks() const
            {
                return m_IntBlocks;
            }
        };


        namespace 
        {
            typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
        }

        MULTI_REGIONS_EXPORT void CuthillMckeeReordering(const BoostGraph& graph,
                                    Array<OneD, int>& perm,
                                    Array<OneD, int>& iperm);

        MULTI_REGIONS_EXPORT void MultiLevelBisectionReordering(const BoostGraph& graph,
                                           const Array<OneD, const int>& vwgts,
                                           Array<OneD, int>& perm,
                                           Array<OneD, int>& iperm,
                                           BottomUpSubStructuredGraphSharedPtr& substructgraph,
                                           const int mdswitch = 1);
        // The parameter MDSWITCH.
        // This parameters defines the maximal size of the smallest patches.
        // If at a certain level, a patch bundles less than MDSWITCH graph-vertices,
        // metis is not going to partition this subgraph any further.
        // Some quick and basis test have shown that 30 seems to be good value.
        // However, this optimal value will probably depend on the polynomial order
        // of the expansion and there is still room for optimisation here.

        MULTI_REGIONS_EXPORT void NoReordering(const BoostGraph& graph,
                          Array<OneD, int>& perm,
                          Array<OneD, int>& iperm);



        
    } // end of namespace
} // end of namespace

#endif // MULTIREGIONS_SUBSTRUCTUREDGRAPH_H

/**
* $Log: SubStructuredGraph.h,v $
* Revision 1.4  2009/11/19 11:41:07  pvos
* Fixed various bugs
*
* Revision 1.3  2009/11/09 15:57:11  pvos
* multi-level recursion bug fixes
*
* Revision 1.2  2009/11/02 11:19:44  pvos
* Fixed a bug for reordering a graph without edges
*
* Revision 1.1  2009/10/30 14:02:55  pvos
* Multi-level static condensation updates
*
**/
