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

#ifndef EXPT_TESTS_STAT_H
#define EXPT_TESTS_STAT_H

#include <math.h>
#include <boost/math/distributions/students_t.hpp>

class Stat
{
    public:
        Stat(const std::vector<double>& samples, double confidence) :
            Mean(0.0),
            HalfWidth(0.0),
            Confidence(0.0),
            StdDev(0.0)
        {
            Confidence = confidence;
            double sum = 0.0;
            for(unsigned int i = 0; i < samples.size(); ++i)
            {
                sum += samples[i];
            }
            Mean = sum/(double)samples.size();

            StdDev = 0.0;
            for(unsigned int i = 0; i < samples.size(); ++i)
            {
                StdDev += (samples[i] - Mean)*(samples[i]-Mean);
            }
            StdDev = sqrt(StdDev/(samples.size()-1));

            boost::math::students_t dist(samples.size() - 1);
            double T = boost::math::quantile(boost::math::complement(dist, (1.0-Confidence)/2.0));
            HalfWidth = T*StdDev/sqrt(static_cast<double>(samples.size()));
        }

        bool Overlaps(const Stat& other)
        {
            return (Low() >= other.Low() && Low() <= other.High()) ||
                (High() >= other.Low() && High() <= other.High());
        }

        double Low() const { return Mean - HalfWidth; }
        double High() const { return Mean + HalfWidth; }
        double Mean;
        double HalfWidth;
        double Confidence;
        double StdDev;
};

#endif
