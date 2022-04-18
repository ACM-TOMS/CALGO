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

#ifndef EXPT_TESTS_CONJUGATE_GRADIENT_HPP
#define EXPT_TESTS_CONJUGATE_GRADIENT_HPP

#include <iostream>

namespace expt_tests
{
    template<typename MatrixType, typename VectorType, unsigned int Id>
    void RunConjugateGradientTest(unsigned int numIterations, unsigned int problemSize, boost::timer::cpu_timer& timer)
    {
        MatrixType A(problemSize, problemSize);
        VectorType b(problemSize);
        for(unsigned int i = 0; i < problemSize; ++i)
        {
            b[i] = i+1;
            for(unsigned int j = 0; j < problemSize; ++j)
            {
                A(i,j) = i*problemSize + j + 1;
            }
        }

        VectorType r_k(problemSize);
        VectorType r_k1(problemSize);
        VectorType x_k(problemSize);
        VectorType x_k1(problemSize);
        VectorType p_k(problemSize);
        VectorType p_k1(problemSize);

        

        timer.start();

        for(unsigned int j = 0; j < numIterations; ++j)
        {
            //for(unsigned int i = 0; i < problemSize; ++i)
            //{
            //    x_k[i] = 0.0;
            //}
            //r_k = b - A*x_k;
            r_k = b - PROD(A,x_k);
            p_k = r_k;

            for(unsigned int k = 0; k < problemSize; ++k)
            {
                // TODO - restore before testing.
                //double alpha = (r_k * r_k)/(p_k*A*p_k);
                //double num = r_k * r_k;
                double num = DOT(r_k, r_k);

                //double denom = p_k*(A*p_k);
                double denom = DOT(p_k, (PROD(A,p_k)));

                double alpha = num/denom;
                x_k1 = x_k + alpha*p_k;

                r_k1 = r_k - alpha*(PROD(A,p_k));

                double beta = DOT( r_k1, r_k1)/DOT(r_k, r_k);
                p_k1 = r_k1 + beta*p_k;

                // Setup for next iteration.
                r_k = r_k1;
                x_k = x_k1;
                p_k = p_k1;
            }
        }
        timer.stop();
    }

    class NektarConjugateGradientWithoutExpressionTemplates : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };

    class NektarConjugateGradientWithExpressionTemplates : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };

    class UblasConjugateGradient : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };

    class EigenConjugateGradient : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };

    class ArmadilloConjugateGradient : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };


}

#endif
