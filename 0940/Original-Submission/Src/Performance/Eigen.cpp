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

#define EIGEN_DONT_VECTORIZE

#include "VectorAdditionTests.h"
#include "MatrixMultiplication.h"



#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

//#include "ThirdParty/eigen-eigen-65ee2328342f/Eigen/Dense"
#include "ThirdParty/eigen3.12/Eigen/Dense"
 
#define DOT(a,b) a.dot(b)
#define PROD(a,b) (a*b)

#include "ConjugateGradient.hpp"


namespace expt_tests
{
    void VectorAdditionEigen::RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
    {
        RunVectorAdditionTest<Eigen::VectorXf, 1>(iterations, numInputs, size, t);
    }

    std::string VectorAdditionEigen::TestName()
    {
        return "VectorAdditionEigen";
    }

    void EigenMultiplication::RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
    {
        RunMatrixMultiplicationTest<Eigen::MatrixXf, 1, CreateDefaultMatrix<Eigen::MatrixXf> >(iterations, numInputs, size, t);
    }

    std::string EigenMultiplication::TestName()
    {
        return "EigenMultiplication";
    }

    void EigenConjugateGradient::RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
    {
        RunConjugateGradientTest<Eigen::MatrixXf, Eigen::VectorXf, 1>(iterations, size, t);
    }

    std::string EigenConjugateGradient::TestName()
    {
        return "EigenConjugateGradient";
    }

}
