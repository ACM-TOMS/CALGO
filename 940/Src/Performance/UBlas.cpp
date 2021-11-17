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



#include "VectorAdditionTests.h"
#include "MatrixMultiplication.h"
#include "UBlas.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#define DOT(a,b) (inner_prod(a,b))
#define PROD(a,b) (prod(a,b))

#include "ConjugateGradient.hpp"

namespace expt_tests
{
    void VectorAdditionUBlas::RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
    {
        RunVectorAdditionTest<boost::numeric::ublas::vector<double>,1 >(iterations, numInputs, size, t);
    }

    std::string VectorAdditionUBlas::TestName()
    {
        return "VectorAdditionUBlas";
    }


    void UblasConjugateGradient::RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
    {
        RunConjugateGradientTest<boost::numeric::ublas::matrix<double>, boost::numeric::ublas::vector<double>, 1>(iterations, size, t);
    }

    std::string UblasConjugateGradient::TestName()
    {
        return "UblasConjugateGradient";
    }


    void UblasMatrixMultiplication::RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
    {
        typedef boost::numeric::ublas::matrix<double> MatrixType;
        switch(numInputs)
        {
            case 2:
            {
                MatrixType v0(size, size);
                MatrixType v1(size, size);
                MatrixType result(size, size);

                // warmup
                result = prod(v0,v1);

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    result = prod(v0,v1);
                }
                t.stop();
            }
            break;

            case 3:
            {
                MatrixType v0(size, size);
                MatrixType v1(size, size);
                MatrixType v2(size, size);
                MatrixType result(size, size);

                result = prod(v0, MatrixType(prod(v1,v2)));

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    result =prod(v0, MatrixType(prod(v1,v2)));
                }
                t.stop();
            }
            break;

            case 4:
            {
                MatrixType v0(size, size);
                MatrixType v1(size, size);
                MatrixType v2(size, size);
                MatrixType v3(size, size);
                MatrixType result(size, size);

                result = prod(v0, MatrixType(prod(v1, MatrixType(prod(v2,v3)))));

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    result =prod(v0, MatrixType(prod(v1, MatrixType(prod(v2,v3)))));
                }
                t.stop();
            }
            break;

            case 5:
            {
                MatrixType v0(size, size);
                MatrixType v1(size, size);
                MatrixType v2(size, size);
                MatrixType v3(size, size);
                MatrixType v4(size, size);
                MatrixType result(size, size);

                result = prod(v0, MatrixType(prod(v1, MatrixType(prod(v2, MatrixType(prod(v3,v4)))))));

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    result =prod(v0, MatrixType(prod(v1, MatrixType(prod(v2, MatrixType(prod(v3,v4)))))));
                }
                t.stop();
            }
            break;

            case 8:
            {
//                MatrixType v0(size, size);
//                MatrixType v1(size, size);
//                MatrixType v2(size, size);
//                MatrixType v3(size, size);
//                MatrixType v4(size, size);
//                MatrixType v5(size, size);
//                MatrixType v6(size, size);
//                MatrixType v7(size, size);
//                MatrixType result(size, size);

//                result = v0*v1*v2*v3*v4*v5*v6*v7;

//                t.start();
//                for(unsigned int i = 0; i < iterations; ++i)
//                {
//                    result = v0*v1*v2*v3*v4*v5*v6*v7;
//                }
//                t.stop();
            }
            break;
        }
    }

    std::string UblasMatrixMultiplication::TestName() { return "UblasMatrixMultiplication"; }

}
