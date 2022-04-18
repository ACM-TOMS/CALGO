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

#include "HandCoded.h"
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

namespace expt_tests
{

    void HandCodedVectorAddition::RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
    {
        switch(numInputs)
        {
            case 2:
            {
                double* v0 = new double[size];
                double* v1 = new double[size];
                double* result = new double[size];

                for(unsigned int i =0; i < size; ++i)
                {
                  result[i] = v0[i] + v1[i];
                }

                t.start();
                for(unsigned int j = 0; j < iterations; ++j)
                {
                  for(unsigned int i =0; i < size; ++i)
                  {
                    result[i] = v0[i] + v1[i];
                  }
                }
                t.stop();

                delete[] v0;
                delete[] v1;
                delete[] result;
            }
            break;

            case 4:
            {
                double* v0 = new double[size];
                double* v1 = new double[size];
                double* v2 = new double[size];
                double* v3 = new double[size];
                double* result = new double[size];

                for(unsigned int i =0; i < size; ++i)
                {
                  result[i] = v0[i] + v1[i] + v2[i] + v3[i];
                }

                t.start();
                for(unsigned int j = 0; j < iterations; ++j)
                {
                  for(unsigned int i =0; i < size; ++i)
                  {
                    result[i] = v0[i] + v1[i] + v2[i] + v3[i];
                  }
                }
                t.stop();

                delete[] v0;
                delete[] v1;
                delete[] v2;
                delete[] v3;
                delete[] result;
            }
            break;

            case 8:
            {
                double* v0 = new double[size];
                double* v1 = new double[size];
                double* v2 = new double[size];
                double* v3 = new double[size];
                double* v4 = new double[size];
                double* v5 = new double[size];
                double* v6 = new double[size];
                double* v7 = new double[size];
                double* result = new double[size];

                for(unsigned int i =0; i < size; ++i)
                {
                  result[i] = v0[i] + v1[i] + v2[i] + v3[i] + v4[i] + v5[i] + v6[i] + v7[i];
                }

                t.start();
                for(unsigned int j = 0; j < iterations; ++j)
                {
                  for(unsigned int i =0; i < size; ++i)
                  {
                    result[i] = v0[i] + v1[i] + v2[i] + v3[i] + v4[i] + v5[i] + v6[i] + v7[i];
                  }
                }
                t.stop();

                delete[] v0;
                delete[] v1;
                delete[] v2;
                delete[] v3;
                delete[] v4;
                delete[] v5;
                delete[] v6;
                delete[] v7;
                delete[] result;
            }
            break;
        }
    }

    std::string HandCodedVectorAddition::TestName()
    {
      return "HandCodedVectorAddition";
    }
 
    void HandCodedMatrixMultiplication::RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
    {
        switch(numInputs)
        {
            //case 2:
            //{
            //    MatrixType v0 = CreateFunc::CreateMatrix(size, size);
            //    MatrixType v1 = CreateFunc::CreateMatrix(size, size);
            //    MatrixType result = CreateFunc::CreateMatrix(size, size);

            //    // warmup
            //    result = v0*v1;

            //    t.start();
            //    for(unsigned int i = 0; i < iterations; ++i)
            //    {
            //        result = v0*v1;
            //    }
            //    t.stop();
            //}
            //break;

            case 3:
            {
              NekMatrix<double> m0(size, size);
              NekMatrix<double> m1(size, size);
              NekMatrix<double> m2(size, size);
              NekMatrix<double> temp(size, size);
              NekMatrix<double> result(size, size);

              Multiply(temp, m0, m1);
              Multiply(result, temp, m2);

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                  Multiply(temp, m0, m1);
                  Multiply(result, temp, m2);
                }
                t.stop();
            }
            break;

            case 4:
            {
              NekMatrix<double> m0(size, size);
              NekMatrix<double> m1(size, size);
              NekMatrix<double> m2(size, size);
              NekMatrix<double> m3(size, size);
              NekMatrix<double> temp(size, size);
              NekMatrix<double> result(size, size);

              Multiply(result, m0, m1);
              Multiply(temp, result, m2);
              Multiply(result, temp, m3);

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                  Multiply(result, m0, m1);
                  Multiply(temp, result, m2);
                  Multiply(result, temp, m3);
                }
                t.stop();
            }
            break;

            //case 5:
            //{
            //    MatrixType v0 = CreateFunc::CreateMatrix(size, size);
            //    MatrixType v1 = CreateFunc::CreateMatrix(size, size);
            //    MatrixType v2 = CreateFunc::CreateMatrix(size, size);
            //    MatrixType v3 = CreateFunc::CreateMatrix(size, size);
            //    MatrixType v4 = CreateFunc::CreateMatrix(size, size);
            //    MatrixType result = CreateFunc::CreateMatrix(size, size);

            //    result = v0*v1*v2*v3*v4;

            //    t.start();
            //    for(unsigned int i = 0; i < iterations; ++i)
            //    {
            //        result = v0*v1*v2*v3*v4;
            //    }
            //    t.stop();
            //}
            //break;

            //case 8:
            //{
            //    MatrixType v0 = CreateFunc::CreateMatrix(size, size);
            //    MatrixType v1 = CreateFunc::CreateMatrix(size, size);
            //    MatrixType v2 = CreateFunc::CreateMatrix(size, size);
            //    MatrixType v3 = CreateFunc::CreateMatrix(size, size);
            //    MatrixType v4 = CreateFunc::CreateMatrix(size, size);
            //    MatrixType v5 = CreateFunc::CreateMatrix(size, size);
            //    MatrixType v6 = CreateFunc::CreateMatrix(size, size);
            //    MatrixType v7 = CreateFunc::CreateMatrix(size, size);
            //    MatrixType result = CreateFunc::CreateMatrix(size, size);

            //    result = v0*v1*v2*v3*v4*v5*v6*v7;

            //    t.start();
            //    for(unsigned int i = 0; i < iterations; ++i)
            //    {
            //        result = v0*v1*v2*v3*v4*v5*v6*v7;
            //    }
            //    t.stop();
            //}
            //break;
        }
    }

    std::string HandCodedMatrixMultiplication::TestName()
    {
      return "HandCodedMatrixMultiplication";
    }

}
