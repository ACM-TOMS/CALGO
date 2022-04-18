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

#ifndef EXPT_TESTS_TIMING_TEST_RUNNER_H
#define EXPT_TESTS_TIMING_TEST_RUNNER_H

#include "ITest.h"
#include <boost/shared_ptr.hpp>
#include <vector>
#include "Stat.h"

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

namespace expt_tests
{

    // Ideally, I'll have a vector of function objects that I can call one at a time.
    // Implemented as derived classes?  But that introduces virtual function overhead to
    // my timing tests.
    //
    // A boost fusion vector of objects will be better
    // actually if I push the iterations on the subclass, there is only one call to the virtual function.
    std::vector<Stat> RunTest(std::vector<boost::shared_ptr<ITest> >& tests, int problemSize, int problemArgs, double confidenceLevel, unsigned int iterations)
    {
        std::vector<Stat> result;
        std::vector<std::vector<double> > data;
        for(unsigned int i = 0; i < tests.size(); ++i)
        {
            data.push_back(std::vector<double>());
        }

        // requested time is in seconds, convert to nanoseconds for testing.
//        boost::timer::nanosecond_type timePerRepAsNano = timePerReplication*1e8;

//        iterations = 0;
//        for(unsigned int i = 0; i < tests.size(); ++i)
//        {
//            unsigned int numTrials = 100;
//            unsigned int localIterations = std::numeric_limits<unsigned int>::max();

//            // Sometimes 100 iterations to get a feel for how long the test will take is so short that the
//            // predicted number of iterations is skewed.  This shows up as a very large number of iterations,
//            // so we check for that here.
////            while( localIterations > 100000 )
////            {
//                //std::cout << "Testing " << numTrials << std::endl;
//                boost::timer::cpu_timer t;
//                tests[i]->RunTest(numTrials, problemArgs, problemSize, t);
//                //std::cout << "Elapsed time: " << t.elapsed().wall << std::endl;
//                //std::cout << "Numerator: " << timePerRepAsNano*numTrials << std::endl;
//                //std::cout << "Division: " << (timePerRepAsNano*numTrials)/t.elapsed().wall << std::endl;
//                localIterations = static_cast<unsigned int>((timePerRepAsNano*numTrials)/t.elapsed().wall);
//                //std::cout << "Local iterations: " << localIterations << std::endl;
////            }

//            iterations = std::max(localIterations, iterations);
//            std::cout << "Number of iterations = " << iterations << std::endl;
//        }

        unsigned int minReplications = 10;
        unsigned int maxReplications = 12;
        bool confidenceIntervalsOverlap = true;
        unsigned int curReplication = 0;
        while( curReplication < minReplications || (confidenceIntervalsOverlap && curReplication < maxReplications) )
        {
            std::cout << "Running replication " << curReplication << std::endl;

            for(unsigned int i = 0; i < tests.size(); ++i)
            {
                boost::timer::cpu_timer t;
                tests[i]->RunTest(iterations, problemArgs, problemSize, t);
                data[i].push_back(t.elapsed().wall/1e9);
                std::cout << t.elapsed().wall/1e9 << std::endl;
            }

            if( curReplication > minReplications )
            {
                // Test the confidence intervals to see if we can stop.
                result.clear();
                for(unsigned int i = 0; i < tests.size(); ++i)
                {
                    result.push_back(Stat(data[i], confidenceLevel));
                }

                confidenceIntervalsOverlap = false;
                for(unsigned int i = 0; i < tests.size(); ++i)
                {
                    for(unsigned int j = 0; j < tests.size(); ++j)
                    {
                        if( i == j ) continue;
                        if( result[i].Overlaps(result[j]) )
                        {
                            std::cout << "Implementation " << i << " and " << j << " have overlapping confidence intervals." << std::endl;
                            std::cout << "[" << result[i].Low() << ", " << result[i].High() << "] - ["
                                << result[j].Low() << ", " << result[j].High() << "]" << std::endl;
                            confidenceIntervalsOverlap = true;
                            break;
            }
                    }
                    if( confidenceIntervalsOverlap ) break;
                }
            }
            curReplication += 1;
        }

        std::cout << "Entries = " << result.size() << std::endl;
        return result;
    }
}

#endif

