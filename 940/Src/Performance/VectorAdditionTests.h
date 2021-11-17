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

#ifndef EXPT_TESTS_VECTOR_ADDITION_H
#define EXPT_TESTS_VECTOR_ADDITION_H

#include "ITest.h"

namespace expt_tests
{

    template<typename VectorType, unsigned int Id>
    void RunVectorAdditionTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
    {
        switch(numInputs)
        {
            case 2:
            {
                VectorType v0(size);
                VectorType v1(size);
                VectorType result(size);

                // warmup
                result = v0+v1;

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    result = v0+v1;
                }
                t.stop();
            }
            break;

            case 4:
            {
                VectorType v0(size);
                VectorType v1(size);
                VectorType v2(size);
                VectorType v3(size);
                VectorType result(size);

                result = v0+v1+v2+v3;

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    result = v0+v1+v2+v3;
                }
                t.stop();
            }
            break;

            case 8:
            {
                VectorType v0(size);
                VectorType v1(size);
                VectorType v2(size);
                VectorType v3(size);
                VectorType v4(size);
                VectorType v5(size);
                VectorType v6(size);
                VectorType v7(size);
                VectorType result(size);

                result = v0+v1+v2+v3+v4+v5+v6+v7;

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    result = v0+v1+v2+v3+v4+v5+v6+v7;
                }
                t.stop();
            }
            break;
        }
    }

    class LocalVectorAdditionWithExpressionTemplates : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };

    class LocalVectorAdditionWithoutExpressionTemplates : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };


    class NektarVectorAdditionWithExpressionTemplates : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };

    class NektarVectorAdditionWithoutExpressionTemplates : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };


    class VectorAdditionBlitz : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };


    class VectorAdditionArmadillo : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };



    class VectorAdditionUBlas : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };

    class VectorAdditionGenial : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };

    class VectorAdditionEigen : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };
}

#endif
