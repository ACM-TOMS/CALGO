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

#include "ITest.h"
#include "HandCoded.h"
#include "TimingTestRunner.h"
#include "NektarMoveConstructors.h"
#include "MET.h"
#include "VectorAdditionTests.h"
#include "MatrixMultiplication.h"
#include "ConjugateGradient.hpp"
#include "Sets.hpp"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <set>
#include "UBlas.h"

#include <boost/shared_ptr.hpp>



int main(int argc, char** argv)
{
    if( !(argc == 1 || argc == 5 || argc == 6 ) )
    {
        std::cout << "Usage: Performance <TestName> <ProblemSize> <ProblemArgs> <NumTests> <ConfidenceInterval>" << std::endl;
        return 1;
    }

    int problemSize = 100;
    int problemArgs = 2;
    int numTests = 100;
    double confidenceInterval = .9;
    const char* testName = "ConjugateGradient";

    if( argc > 1 )
    {
        problemSize = atoi(argv[2]);
        problemArgs = atoi(argv[3]);
        numTests = atoi(argv[4]);
        confidenceInterval = .9;
        if( argc == 6 )
        {
            confidenceInterval = atof(argv[5]);
        }
        testName = argv[1];
    }

    std::vector<boost::shared_ptr<expt_tests::ITest> > tests;

    if( strcmp(testName, "VectorAddition") == 0 )
    {
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::NektarVectorAdditionWithoutExpressionTemplates()));
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::HandCodedVectorAddition()));
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::NektarVectorAdditionWithExpressionTemplates()));
        
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::NektarWithMoveConstructorsVectorAddition()));
        
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::VectorAdditionBlitz()));

        // Armadillo doesn't work in 32 bit visual studio, out of stack space.
        #ifndef  _MSC_VER
          tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::VectorAdditionArmadillo()));
        #endif
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::VectorAdditionGenial()));
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::VectorAdditionMET()));
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::VectorAdditionUBlas()));
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::VectorAdditionEigen()));
    }
    else if( strcmp(testName, "MatrixMultiplication") == 0 )
    {
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::MatrixWithoutExpressionTemplates()));
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::HandCodedMatrixMultiplication()));
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::MatrixWithExpressionTemplates()));

        // Really slow - calls operator() for each element during multiplication, which is slow in Nektar++.
        //tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::MetWrappedMatrixMultiplication()));
        
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::ArmadilloMultiplication()));
        
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::EigenMultiplication()));

        // From a practical standpoint, these tests don't matter.  They don't use BLAS for the multiplication
        // and therefore the performance is signifcantly worse than any test above.
        //tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::UblasMatrixMultiplication()));
        //tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::MetMatrixMultiplication()));
    }
    else if( strcmp(testName, "ConjugateGradient") == 0 )
    {
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::NektarConjugateGradientWithoutExpressionTemplates()));
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::NektarConjugateGradientWithExpressionTemplates()));
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::UblasConjugateGradient()));
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::EigenConjugateGradient()));
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::ArmadilloConjugateGradient()));
    }
    else if( strcmp(testName, "Sets") == 0 )
    {
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::BaselineSetUnionTest()));
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::ExpTempSetUnionTest()));
    }
    else if( strcmp(testName, "SetExpression") == 0 )
    {
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::BaselineSetExpressionTest()));
        tests.push_back(boost::shared_ptr<expt_tests::ITest>(new expt_tests::ExpTempSetExpressionTest()));
    }
    else
    {
        std::cout << "Unkown tests " << testName << std::endl;
        return 1;
    }

    std::cout << "Running tests." << std::endl;
    unsigned int iterations = 0;
    std::vector<Stat> runResults = expt_tests::RunTest(tests, problemSize, problemArgs, confidenceInterval, numTests);

    std::cout << "Final Results" << std::endl;
    std::cout << "Problem Size: " << problemSize << std::endl;
    std::cout << "Problem Arguments: " << problemArgs << std::endl;
    std::cout << "Confidence Interval: "<< confidenceInterval << std::endl;
    std::cout << "Iterations: " << numTests << std::endl;
    for(unsigned int j = 0; j < runResults.size(); ++j)
    {
        std::cout << tests[j]->TestName() << ": " << runResults[j].Mean/numTests << std::endl;
    }
    for(unsigned int j = 0; j < runResults.size(); ++j)
    {
        double baselineAvg = runResults[0].Mean/numTests;
        double thisAvg = runResults[j].Mean/numTests;
        std::cout << "(" << problemSize << ", " << runResults[j].Mean/numTests << ")" << std::endl;
    }
    std::cout << "Speedup." << std::endl;
    for(unsigned int j = 0; j < runResults.size(); ++j)
    {
        double baselineAvg = runResults[0].Mean/numTests;
        double thisAvg = runResults[j].Mean/numTests;
        std::cout << "(" << problemSize << ", " << baselineAvg/thisAvg << ")" << std::endl;
    }
    return 0;
}
