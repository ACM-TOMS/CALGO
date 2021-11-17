///////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
#undef NEKTAR_USE_EXPRESSION_TEMPLATES
#endif


#include "VectorAdditionTests.h"
#include "MatrixMultiplication.h"

#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/StandardMatrix.hpp>
#include <LibUtilities/LinearAlgebra/MatrixOperations.hpp>

#define DOT(a,b) (Dot(a,b))
#define PROD(a,b) (a*b)
#include "ConjugateGradient.hpp"

namespace expt_tests
{
    void NektarVectorAdditionWithoutExpressionTemplates::RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
    {
        RunVectorAdditionTest<Nektar::NekVector<double>,0>(iterations, numInputs, size, t);
    }

    std::string NektarVectorAdditionWithoutExpressionTemplates::TestName()
    {
        return "VectorAdditionWithoutExpressionTemplates";
    }

    void MatrixWithoutExpressionTemplates::RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
    {
        RunMatrixMultiplicationTest<Nektar::NekMatrix<double>, 0, CreateDefaultMatrix<Nektar::NekMatrix<double> > >(iterations, numInputs, size, t);
    }

    std::string MatrixWithoutExpressionTemplates::TestName()
    {
        return "MatrixWithoutExpressionTemplates";
    }

    void NektarConjugateGradientWithoutExpressionTemplates::RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
    {
        RunConjugateGradientTest<Nektar::NekMatrix<double>, Nektar::NekVector<double>, 0>(iterations, size, t);
    }

    std::string NektarConjugateGradientWithoutExpressionTemplates::TestName()
    {
        return "NektarConjugateGradientWithoutExpressionTemplates";
    }
}
