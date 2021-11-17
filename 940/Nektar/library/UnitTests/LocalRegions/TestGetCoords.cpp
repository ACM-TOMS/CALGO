///////////////////////////////////////////////////////////////////////////////
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/HexExp.h>
#include <SpatialDomains/MeshGraph.h>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
    namespace HexExpTests
    {
        SpatialDomains::SegGeomSharedPtr CreateSegGeom(unsigned int id, 
            SpatialDomains::VertexComponentSharedPtr v0,
            SpatialDomains::VertexComponentSharedPtr v1)
        {
            SpatialDomains::VertexComponentSharedPtr vertices[] = {v0, v1};
            SpatialDomains::SegGeomSharedPtr result(new SpatialDomains::SegGeom(id, 3, vertices));
            return result;
        }

        SpatialDomains::HexGeomSharedPtr CreateHex(
            SpatialDomains::VertexComponentSharedPtr v0,
            SpatialDomains::VertexComponentSharedPtr v1,
            SpatialDomains::VertexComponentSharedPtr v2,
            SpatialDomains::VertexComponentSharedPtr v3,
            SpatialDomains::VertexComponentSharedPtr v4,
            SpatialDomains::VertexComponentSharedPtr v5,
            SpatialDomains::VertexComponentSharedPtr v6,
            SpatialDomains::VertexComponentSharedPtr v7)
        {
            Nektar::SpatialDomains::SegGeomSharedPtr e0 = CreateSegGeom(0, v0, v1);
            Nektar::SpatialDomains::SegGeomSharedPtr e1 = CreateSegGeom(1, v1, v2);
            Nektar::SpatialDomains::SegGeomSharedPtr e2 = CreateSegGeom(2, v2, v3);
            Nektar::SpatialDomains::SegGeomSharedPtr e3 = CreateSegGeom(3, v3, v0);
            Nektar::SpatialDomains::SegGeomSharedPtr e4 = CreateSegGeom(4, v0, v4);
            Nektar::SpatialDomains::SegGeomSharedPtr e5 = CreateSegGeom(5, v1, v5);
            Nektar::SpatialDomains::SegGeomSharedPtr e6 = CreateSegGeom(6, v2, v6);
            Nektar::SpatialDomains::SegGeomSharedPtr e7 = CreateSegGeom(7, v3, v7);
            Nektar::SpatialDomains::SegGeomSharedPtr e8 = CreateSegGeom(8, v4, v5);
            Nektar::SpatialDomains::SegGeomSharedPtr e9 = CreateSegGeom(9, v5, v6);
            Nektar::SpatialDomains::SegGeomSharedPtr e10 = CreateSegGeom(10, v6, v7);
            Nektar::SpatialDomains::SegGeomSharedPtr e11 = CreateSegGeom(11, v4, v7);

            Nektar::SpatialDomains::SegGeomSharedPtr edgesF0[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e0, e1, e2, e3
            };
            Nektar::SpatialDomains::SegGeomSharedPtr edgesF1[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e0, e5, e8, e4
            };
            Nektar::SpatialDomains::SegGeomSharedPtr edgesF2[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e1, e6, e9, e5
            };
            Nektar::SpatialDomains::SegGeomSharedPtr edgesF3[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e2, e6, e10, e7
            };
            Nektar::SpatialDomains::SegGeomSharedPtr edgesF4[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e3, e7, e11, e4
            };
            Nektar::SpatialDomains::SegGeomSharedPtr edgesF5[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e8, e9, e10, e11
            };

            Nektar::StdRegions::Orientation edgeorient0[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF0[0], *edgesF0[1]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF0[1], *edgesF0[2]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF0[2], *edgesF0[3]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF0[3], *edgesF0[0])
            };
            Nektar::StdRegions::Orientation edgeorient1[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF1[0], *edgesF1[1]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF1[1], *edgesF1[2]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF1[2], *edgesF1[3]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF1[3], *edgesF1[0])
            };
            Nektar::StdRegions::Orientation edgeorient2[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF2[0], *edgesF2[1]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF2[1], *edgesF2[2]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF2[2], *edgesF2[3]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF2[3], *edgesF2[0])
            };
            Nektar::StdRegions::Orientation edgeorient3[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF3[0], *edgesF3[1]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF3[1], *edgesF3[2]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF3[2], *edgesF3[3]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF3[3], *edgesF3[0])
            };
            Nektar::StdRegions::Orientation edgeorient4[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF4[0], *edgesF4[1]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF4[1], *edgesF4[2]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF4[2], *edgesF4[3]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF4[3], *edgesF4[0])
            };
            Nektar::StdRegions::Orientation edgeorient5[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF5[0], *edgesF5[1]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF5[1], *edgesF5[2]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF5[2], *edgesF5[3]),
                Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edgesF5[3], *edgesF5[0])
            };

            Nektar::SpatialDomains::QuadGeomSharedPtr face0(new SpatialDomains::QuadGeom(0, edgesF0, edgeorient0));
            Nektar::SpatialDomains::QuadGeomSharedPtr face1(new SpatialDomains::QuadGeom(1, edgesF1, edgeorient1));
            Nektar::SpatialDomains::QuadGeomSharedPtr face2(new SpatialDomains::QuadGeom(2, edgesF2, edgeorient2));
            Nektar::SpatialDomains::QuadGeomSharedPtr face3(new SpatialDomains::QuadGeom(3, edgesF3, edgeorient3));
            Nektar::SpatialDomains::QuadGeomSharedPtr face4(new SpatialDomains::QuadGeom(4, edgesF4, edgeorient4));
            Nektar::SpatialDomains::QuadGeomSharedPtr face5(new SpatialDomains::QuadGeom(5, edgesF5, edgeorient5));

            Nektar::SpatialDomains::QuadGeomSharedPtr qfaces[] = {face0, face1, face2, face3, face4, face5};
            SpatialDomains::HexGeomSharedPtr hexGeom(new SpatialDomains::HexGeom(qfaces));
            return hexGeom;
        }

        BOOST_AUTO_TEST_CASE(TestHexExpThatIsStdRegion)
        {
            SpatialDomains::VertexComponentSharedPtr v0(new SpatialDomains::VertexComponent(3u, 0u, -1.0, -1.0, -1.0));
            SpatialDomains::VertexComponentSharedPtr v1(new SpatialDomains::VertexComponent(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::VertexComponentSharedPtr v2(new SpatialDomains::VertexComponent(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::VertexComponentSharedPtr v3(new SpatialDomains::VertexComponent(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::VertexComponentSharedPtr v4(new SpatialDomains::VertexComponent(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::VertexComponentSharedPtr v5(new SpatialDomains::VertexComponent(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::VertexComponentSharedPtr v6(new SpatialDomains::VertexComponent(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::VertexComponentSharedPtr v7(new SpatialDomains::VertexComponent(3u, 7u, -1.0, 1.0, 1.0));
            
            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);
            
            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr hexExp = 
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);

            Array<OneD, NekDouble> c0 = Array<OneD, NekDouble>(hexExp->GetTotPoints());
            Array<OneD, NekDouble> c1 = Array<OneD, NekDouble>(hexExp->GetTotPoints());
            Array<OneD, NekDouble> c2 = Array<OneD, NekDouble>(hexExp->GetTotPoints());
            hexExp->GetCoords(c0, c1, c2);
            boost::shared_ptr<StdRegions::StdHexExp> stdHex =
                boost::dynamic_pointer_cast<StdRegions::StdHexExp>(hexExp);
            stdHex->GetCoords(c0, c1, c2);
            double epsilon = 1.0e-8;
            BOOST_CHECK_CLOSE(c0[0], -1.0, epsilon);
            BOOST_CHECK_CLOSE(c0[1], -0.76505532392946474, epsilon);
            BOOST_CHECK_CLOSE(c0[2], -0.28523151648064510, epsilon);
            BOOST_CHECK_CLOSE(c0[3], 0.28523151648064510, epsilon);
            BOOST_CHECK_CLOSE(c0[4], 0.76505532392946474, epsilon);
            BOOST_CHECK_CLOSE(c0[5], 1.0, epsilon);
        }
        
        BOOST_AUTO_TEST_CASE(TestScaledAndTranslatedHexExp)
        {
            SpatialDomains::VertexComponentSharedPtr v0(new SpatialDomains::VertexComponent(3u, 0u, 0.0, 0.0, 0.0));
            SpatialDomains::VertexComponentSharedPtr v1(new SpatialDomains::VertexComponent(3u, 1u, 0.5, 0.0, 0.0));
            SpatialDomains::VertexComponentSharedPtr v2(new SpatialDomains::VertexComponent(3u, 2u, 0.5, 0.5, 0.0));
            SpatialDomains::VertexComponentSharedPtr v3(new SpatialDomains::VertexComponent(3u, 3u, 0.0, 0.5, 0.0));
            SpatialDomains::VertexComponentSharedPtr v4(new SpatialDomains::VertexComponent(3u, 4u, 0.0, 0.0, 0.5));
            SpatialDomains::VertexComponentSharedPtr v5(new SpatialDomains::VertexComponent(3u, 5u, 0.5, 0.0, 0.5));
            SpatialDomains::VertexComponentSharedPtr v6(new SpatialDomains::VertexComponent(3u, 6u, 0.5, 0.5, 0.5));
            SpatialDomains::VertexComponentSharedPtr v7(new SpatialDomains::VertexComponent(3u, 7u, 0.0, 0.5, 0.5));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr hexExp = 
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);

            Array<OneD, NekDouble> c0 = Array<OneD, NekDouble>(hexExp->GetTotPoints());
            Array<OneD, NekDouble> c1 = Array<OneD, NekDouble>(hexExp->GetTotPoints());
            Array<OneD, NekDouble> c2 = Array<OneD, NekDouble>(hexExp->GetTotPoints());
            hexExp->GetCoords(c0, c1, c2);

            double epsilon = 1.0e-8;
            BOOST_CHECK_EQUAL(c0[0], 0.0);
            BOOST_CHECK_CLOSE(c0[1], .05873616902, epsilon);
            BOOST_CHECK_CLOSE(c0[2], .17869212088, epsilon);
            BOOST_CHECK_CLOSE(c0[3], .32130787912, epsilon);
            BOOST_CHECK_CLOSE(c0[4], .44126383098, epsilon);
            BOOST_CHECK_CLOSE(c0[5], .5, epsilon);
        }
    }
}
