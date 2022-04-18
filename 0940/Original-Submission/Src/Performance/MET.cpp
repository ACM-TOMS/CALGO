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
#include "MET.h"
#include <vecmat.h>

namespace expt_tests
{
    void VectorAdditionMET::RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
    {
        RunVectorAdditionTest<Vec<double>,1 >(iterations, numInputs, size, t);
    }

    std::string VectorAdditionMET::TestName()
    {
        return "VectorAdditionMET";
    }

    // MET does not provide a dot product nor the multiplication of two vectors.


    struct CreateMetMatrix
    {
        static Mat<double> CreateMatrix(unsigned int rows, unsigned int cols)
        {
            return Mat<double>(rows);
        }
    };

    void MetMatrixMultiplication::RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
    {
        typedef Mat<double> MatrixType;
        switch(numInputs)
        {
            case 2:
            {
                MatrixType v0 = Mat<double>(size);
                MatrixType v1 = Mat<double>(size);
                MatrixType result = Mat<double>(size);

                // warmup
                result = v0*v1;

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    result = v0*v1;
                }
                t.stop();
            }
            break;

            case 3:
            {
                MatrixType v0 = Mat<double>(size);
                MatrixType v1 = Mat<double>(size);
                MatrixType v2 = Mat<double>(size);
                MatrixType result = Mat<double>(size);

                MatrixType t0 = Mat<double>(size);
                t0 = v0*v1;
                result = t0*v2;

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    t0 = v0*v1;
                    result = t0*v2;
                }
                t.stop();
            }
            break;

            case 4:
            {
                MatrixType v0 = Mat<double>(size);
                MatrixType v1 = Mat<double>(size);
                MatrixType v2 = Mat<double>(size);
                MatrixType v3 = Mat<double>(size);
                MatrixType result = Mat<double>(size);

                MatrixType t0 = Mat<double>(size);
                MatrixType t1 = Mat<double>(size);
                t0 = v0*v1;
                t1 = t0*v2;
                result = t1*v3;

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    t0 = v0*v1;
                    t1 = t0*v2;
                    result = t1*v3;
                }
                t.stop();
            }
            break;

            case 5:
            {
                MatrixType v0 = Mat<double>(size);
                MatrixType v1 = Mat<double>(size);
                MatrixType v2 = Mat<double>(size);
                MatrixType v3 = Mat<double>(size);
                MatrixType v4 = Mat<double>(size);
                MatrixType result = Mat<double>(size);

                MatrixType t0 = Mat<double>(size);
                MatrixType t1 = Mat<double>(size);
                t0 = v0*v1;
                t1 = t0*v2;
                t0 = t1*v3;
                result = t0*v4;

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    t0 = v0*v1;
                    t1 = t0*v2;
                    t0 = t1*v3;
                    result = t0*v4;
                }
                t.stop();
            }
            break;

            case 8:
            {
                MatrixType v0 = Mat<double>(size);
                MatrixType v1 = Mat<double>(size);
                MatrixType v2 = Mat<double>(size);
                MatrixType v3 = Mat<double>(size);
                MatrixType v4 = Mat<double>(size);
                MatrixType v5 = Mat<double>(size);
                MatrixType v6 = Mat<double>(size);
                MatrixType v7 = Mat<double>(size);
                MatrixType result = Mat<double>(size);

                MatrixType t0 = Mat<double>(size);
                MatrixType t1 = Mat<double>(size);
                t0 = v0*v1;
                t1 = t0*v2;
                t0 = t1*v3;
                t1 = t0*v4;
                t0 = t1*v5;
                t1 = t0*v6;
                result = t1*v7;

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    t0 = v0*v1;
                    t1 = t0*v2;
                    t0 = t1*v3;
                    t1 = t0*v4;
                    t0 = t1*v5;
                    t1 = t0*v6;
                    result = t1*v7;
                }
                t.stop();
            }
            break;
        }
    }

    std::string MetMatrixMultiplication::TestName()
    {
        return "MetMatrixMultiplication";
    }

    struct NektarWrappedMatrix : public Dim2<double, NektarWrappedMatrix >
    {
        explicit NektarWrappedMatrix(int r_, int c_) :
                Dim2<double,NektarWrappedMatrix>(),
                Matrix(r_, c_)
            {

            }

          NektarWrappedMatrix(const NektarWrappedMatrix& x) : Dim2<double,NektarWrappedMatrix>(x),
              Matrix(x.Matrix)
          {
            assignFrom(x);
          }

          virtual ~NektarWrappedMatrix() {
          }

          // Assignment from an expression, matrix, scalar
          template <class E>
            NektarWrappedMatrix& operator=(const Xpr2<double,E>& x)
            {
              return assignFrom(x);
            }
//          template <class M>
//            NektarWrappedMatrix& operator=(const Dim2<double,M>& x)
//            {
//              return assignFrom(x);
//            }
//          NektarWrappedMatrix& assignFrom(double x)
//          {
//            return assignFrom(x);
//          }

            int rows() const 
            { 
              return Matrix.GetRows(); 
            }
            
            int cols() const 
            { 
              return Matrix.GetColumns(); 
            }

            double operator()(int i, int j) const 
            { 
              return Matrix(i,j); 
            }
            double& operator()(int i, int j) 
            { 
              return Matrix(i,j); 
            }

        Nektar::NekMatrix<double> Matrix;
    };

    void MetWrappedMatrixMultiplication::RunTest(unsigned int iterations, unsigned int numInputs, unsigned int size, boost::timer::cpu_timer& t)
    {
        typedef NektarWrappedMatrix MatrixType;
        switch(numInputs)
        {
            case 2:
            {
                MatrixType v0 = MatrixType(size, size);
                MatrixType v1 = MatrixType(size, size);
                MatrixType result = MatrixType(size, size);

                // warmup
                result = v0*v1;

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    result = v0*v1;
                }
                t.stop();
            }
            break;

            case 3:
            {
                MatrixType v0 = MatrixType(size, size);
                MatrixType v1 = MatrixType(size, size);
                MatrixType v2 = MatrixType(size, size);
                MatrixType result = MatrixType(size, size);

                MatrixType t0 = MatrixType(size, size);
                t0 = v0*v1;
                result = t0*v2;

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    t0 = v0*v1;
                    result = t0*v2;
                }
                t.stop();
            }
            break;

            case 4:
            {
                MatrixType v0 = MatrixType(size, size);
                MatrixType v1 = MatrixType(size, size);
                MatrixType v2 =MatrixType(size, size);
                MatrixType v3 = MatrixType(size, size);
                MatrixType result = MatrixType(size, size);

                MatrixType t0 = MatrixType(size, size);
                MatrixType t1 = MatrixType(size, size);
                t0 = v0*v1;
                t1 = t0*v2;
                result = t1*v3;

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    t0 = v0*v1;
                    t1 = t0*v2;
                    result = t1*v3;
                }
                t.stop();
            }
            break;

            case 5:
            {
                MatrixType v0 = MatrixType(size, size);
                MatrixType v1 = MatrixType(size, size);
                MatrixType v2 = MatrixType(size, size);
                MatrixType v3 = MatrixType(size, size);
                MatrixType v4 = MatrixType(size, size);
                MatrixType result = MatrixType(size, size);

                MatrixType t0 = MatrixType(size, size);
                MatrixType t1 = MatrixType(size, size);
                t0 = v0*v1;
                t1 = t0*v2;
                t0 = t1*v3;
                result = t0*v4;

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    t0 = v0*v1;
                    t1 = t0*v2;
                    t0 = t1*v3;
                    result = t0*v4;
                }
                t.stop();
            }
            break;

            case 8:
            {
                MatrixType v0 = MatrixType(size, size);
                MatrixType v1 = MatrixType(size, size);
                MatrixType v2 = MatrixType(size, size);
                MatrixType v3 = MatrixType(size, size);
                MatrixType v4 = MatrixType(size, size);
                MatrixType v5 = MatrixType(size, size);
                MatrixType v6 = MatrixType(size, size);
                MatrixType v7 = MatrixType(size, size);
                MatrixType result = MatrixType(size, size);

                MatrixType t0 = MatrixType(size, size);
                MatrixType t1 = MatrixType(size, size);
                t0 = v0*v1;
                t1 = t0*v2;
                t0 = t1*v3;
                t1 = t0*v4;
                t0 = t1*v5;
                t1 = t0*v6;
                result = t1*v7;

                t.start();
                for(unsigned int i = 0; i < iterations; ++i)
                {
                    t0 = v0*v1;
                    t1 = t0*v2;
                    t0 = t1*v3;
                    t1 = t0*v4;
                    t0 = t1*v5;
                    t1 = t0*v6;
                    result = t1*v7;
                }
                t.stop();
            }
            break;
        }
    }

    std::string MetWrappedMatrixMultiplication::TestName()
    {
        return "MetWrappedMatrixMultiplication";
    }
}
