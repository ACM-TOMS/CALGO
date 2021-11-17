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

#ifndef EXPT_TESTS_SETS_HPP
#define EXPT_TESTS_SETS_HPP

#ifndef NEKTAR_USE_EXPRESSION_TEMPLATES
#define NEKTAR_USE_EXPRESSION_TEMPLATES
#endif

#include <iostream>
#include <set>
#include <algorithm>
#include <ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>
#include "ITest.h"
#include <boost/timer/timer.hpp>
#include <boost/typeof/typeof.hpp>
#include  BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()

template<typename T>
std::set<T> operator+(const std::set<T>& lhs, const std::set<T>& rhs)
{
    std::set<T> result;
    std::set_union(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), std::insert_iterator<std::set<T> >(result, result.begin()));
    return result;
}

template<typename T>
std::set<T> operator-(const std::set<T>& lhs, const std::set<T>& rhs)
{
    std::set<T> result;
    std::set_intersection(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), std::insert_iterator<std::set<T> >(result, result.begin()));
    return result;
}

template<typename T>
struct Set
{
    Set() {}
    template<typename T1, typename Op, typename T2>
    Set(const expt::Node<T1, Op, T2>& node)
    {
        expt::ExpressionEvaluator::Evaluate(node, *this);
    }

    std::set<T> Value;
};



template<typename T>
void AddEqual(Set<T>& accumulator, const Set<T>& rhs)
{
    std::set<T> temp;
    std::set_difference(rhs.Value.begin(), rhs.Value.end(), accumulator.Value.begin(), accumulator.Value.end(), std::insert_iterator<std::set<T> >(temp, temp.begin()));
    accumulator.Value.insert(temp.begin(), temp.end());
}

template<typename T>
void Add(Set<T>& accumulator, const Set<T>& lhs, const Set<T>& rhs)
{
    std::set_union(lhs.Value.begin(), lhs.Value.end(), rhs.Value.begin(), rhs.Value.end(), 
      std::insert_iterator<std::set<T> >(accumulator.Value, accumulator.Value.begin()));
}

template<typename T>
Set<T> Add(const Set<T>& lhs, const Set<T>& rhs)
{
   Set<T> result;
   Add(result, lhs, rhs);
   return result;
}



template<typename T>
void SubtractEqual(Set<T>& accumulator, const Set<T>& rhs)
{
    std::set<T> temp;
    std::set_intersection(rhs.Value.begin(), rhs.Value.end(), accumulator.Value.begin(), accumulator.Value.end(), std::insert_iterator<std::set<T> >(temp, temp.begin()));
    accumulator=temp;
}

template<typename T>
void Subtract(Set<T>& accumulator, const Set<T>& lhs, const Set<T>& rhs)
{
    std::set_intersection(lhs.Value.begin(), lhs.Value.end(), rhs.Value.begin(), rhs.Value.end(),
      std::insert_iterator<std::set<T> >(accumulator.Value, accumulator.Value.begin()));
}

template<typename T>
Set<T> Subtract(const Set<T>& lhs, const Set<T>& rhs)
{
   Set<T> result;
   Subtract(result, lhs, rhs);
   return result;
}

BOOST_TYPEOF_REGISTER_TEMPLATE(Set, 1);
GENERATE_ADDITION_OPERATOR(Set, 1, Set, 1);
GENERATE_SUBTRACTION_OPERATOR(Set, 1, Set, 1);

namespace expt_tests
{

    template<typename T>
    void populateSet(std::set<T>& dest, unsigned int count, unsigned int start, unsigned int step)
    {
        for(unsigned int i = 0; i < count; ++i)
        {
            dest.insert(start+i*step);
        }
    }

    template<typename T>
    void populateSet(Set<T>& dest, unsigned int count, unsigned int start, unsigned int step)
    {
        populateSet(dest.Value, count, start, step);
    }

    template<typename SetType>
    void RunSetUnionTest(unsigned int numIterations, unsigned int numInputs, unsigned int problemSize, boost::timer::cpu_timer& timer)
    {
        switch(numInputs)
        {
            case 2:
            {
                SetType s0; populateSet(s0, problemSize, 0, 2);
                SetType s1; populateSet(s1, problemSize, 1, 2);

                timer.start();
                for(unsigned int i = 0; i < numIterations; ++i)
                {
                    SetType result = s0+s1;
                }
                timer.stop();

            }
            break;

            case 3:
            {
                SetType s0; populateSet(s0, problemSize, 0, 3);
                SetType s1; populateSet(s1, problemSize, 1, 3);
                SetType s2; populateSet(s2, problemSize, 2, 3);

                timer.start();
                for(unsigned int i = 0; i < numIterations; ++i)
                {
                    SetType result = (s0+s1)+s2;
                }
                timer.stop();
            }
            break;

            case 4:
            {
                SetType s0; populateSet(s0, problemSize, 0, 4);
                SetType s1; populateSet(s1, problemSize, 1, 4);
                SetType s2; populateSet(s2, problemSize, 2, 4);
                SetType s3; populateSet(s3, problemSize, 3, 4);

                timer.start();
                for(unsigned int i = 0; i < numIterations; ++i)
                {
                    SetType result = ((s0+s1)+s2)+s3;
                }
                timer.stop();
            }
            break;
        }


    }

    class BaselineSetUnionTest : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };

    class ExpTempSetUnionTest : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
            virtual std::string TestName();
    };

    template<typename SetType>
    void RunSetExpressionTest(unsigned int numIterations,unsigned int problemSize, boost::timer::cpu_timer& timer)
    {
        SetType s0; populateSet(s0, problemSize, 0, 3);
        SetType s1; populateSet(s1, problemSize, 1, 3);
        SetType s2; populateSet(s2, problemSize, 2, 3);
        SetType s3; populateSet(s3, problemSize, 0, 3);

        timer.start();
        for(unsigned int i = 0; i < numIterations; ++i)
        {
            SetType result = (s0 + (s1+s2)) - s3;
        }
        timer.stop();

    }

    class BaselineSetExpressionTest : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
            {
                RunSetExpressionTest<std::set<unsigned int> >(iterations, size, t);
            }

            virtual std::string TestName()
            {
                return "BaselineSetExpressionTest";
            }
    };

    class ExpTempSetExpressionTest : public ITest
    {
        public:
            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
            {
                RunSetExpressionTest<Set<unsigned int> >(iterations, size, t);
            }

            virtual std::string TestName()
            {
                return "ExpTempSetExpressionTest";
            }
    };


//    class BaselineSetIntersectionTest : public ITest
//    {
//        public:
//            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
//            virtual std::string TestName();
//    };

//    class ExpTempSetIntersectionTest : public ITest
//    {
//        public:
//            virtual void RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t);
//            virtual std::string TestName();
//    };

}

#endif
