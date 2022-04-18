package jphase;

import static jphase.Utils.*;
//import static java.lang.Math.*;

import java.util.ArrayList;
import java.util.List;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.BiCG;
import no.uib.cipr.matrix.sparse.IterativeSolver;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;

/**
 * Utilities class for the jPhase package
 * @author German Riaño
 * @author Juan F. Pérez
 * @version 1.0
 * 
 */
public class MatrixUtils {

    /**
     * Precision for computations
     */
    static double Epsilon = 1.0E-10;

    /**
     * Stores the result of the function ln(n!), for n =0,...,100
     */
    private static double a[] = new double[101];

    /**
     * Stores the coefficients for the gamma function
     */
    private static double cof[] = { 76.18009172947146, -86.50532032941677,
            24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
            -0.5395239384953e-5 };

    /**
     * Stores the result of the function n!, for n =0,...,32
     */
    private static double fac[] = new double[33];

    /**
     * Keeps the maximum factorial calculated and stored at fac
     */
    private static int topN = 0;

    /**
     * Returns the kronecker sum of two matrices in dense format
     * @param A Matrix
     * @param B Matrix
     * @return Kronecker Sum A + B
     */
    public static Matrix kroneckerSum(Matrix A, Matrix B) {
        int rows1 = A.numRows();
        int rows2 = B.numRows();
        return MatrixUtils.kronecker(A, Matrices.identity(rows2)).add(
                MatrixUtils.kronecker(Matrices.identity(rows1), B));
    }

    /**
     * Returns the kronecker sum of two matrices and stores it in the
     * predefined format
     * @param A Matrix
     * @param B Matrix
     * @param res Result Matrix such that res.numRows = A.numRows *
     *        B.numRows and res.numCols = A.numCols * B.numCols
     * @return Kronecker Sum A + B
     */
    public static Matrix kroneckerSum(Matrix A, Matrix B, Matrix res) {
        int r1 = A.numRows();
        int r2 = B.numRows();
        res = kronecker(A, Matrices.identity(r2), res.copy()).add(
                kronecker(Matrices.identity(r1), B, res.copy()));
        return res;
    }

    /**
     * Returns the kronecker product of two matrices in dense format
     * @param A Matrix
     * @param B Matrix
     * @return Kronecker Product A x B
     */
    public static Matrix kronecker(Matrix A, Matrix B) {
        int r1 = A.numRows();
        int c1 = A.numColumns();
        int r2 = B.numRows();
        int c2 = B.numColumns();
        Matrix result = new DenseMatrix(r1 * r2, c1 * c2);
        for (int i1 = 0; i1 < r1; i1++) {
            for (int j1 = 0; j1 < c1; j1++) {
                for (int i2 = 0; i2 < r2; i2++) {
                    for (int j2 = 0; j2 < c2; j2++) {
                        result.set(i1 * r2 + i2, j1 * c2 + j2, A.get(i1, j1)
                                * B.get(i2, j2));
                    }
                }
            }
        }
        return result;
    }

    /**
     * Returns the Kronecker product of two matrices in the predefined
     * storage format
     * @param A Matrix
     * @param B Matrix
     * @param res Matrix such that res.numCols = A.numCols * B.numCols
     *        and res.numRows = A.numRows * B.numRows
     * @return Kronecker Product A x B
     */
    public static Matrix kronecker(Matrix A, Matrix B, Matrix res) {
        int r1 = A.numRows();
        int c1 = A.numColumns();
        int r2 = B.numRows();
        int c2 = B.numColumns();
        int r3 = res.numRows();
        int c3 = res.numColumns();
        if (r1 * r2 != r3 || c1 * c2 != c3)
            throw new IndexOutOfBoundsException(
                    "Matrices cannot be Kronecker multiplied: (" + r1 + " * "
                            + r2 + " != " + r3 + "  or  " + c1 + " * " + c2
                            + " != " + c3 + ")");

        for (int i1 = 0; i1 < r1; i1++) {
            for (int j1 = 0; j1 < c1; j1++) {
                for (int i2 = 0; i2 < r2; i2++) {
                    for (int j2 = 0; j2 < c2; j2++) {
                        res.set(i1 * r2 + i2, j1 * c2 + j2, A.get(i1, j1)
                                * B.get(i2, j2));

                    }
                }
            }
        }
        return res;
    }

    /**
     * Returns the Kronecker product of one matrix with one vector
     * (Matrix x Vector) in dense format
     * @param A Matrix
     * @param B Vector
     * @param res Matrix such that res.numRows = A.numRows and
     *        res.numCols = A.numCols * B.size
     * @return Kronecker Product A x B
     */
    public static Matrix kronecker(Matrix A, Vector B,
            Matrix res) {
        int r1 = A.numRows();
        int c1 = A.numColumns();
        int n2 = B.size();
        int r3 = res.numRows();
        int c3 = res.numColumns();
        if (r1 * n2 != r3 || c1 != c3)
            throw new IndexOutOfBoundsException(
                    "Matrices cannot be Kronecker multiplied: \n"
                            + "A.numRows * B.size != res.numRows  or  A.numCols  != res.numCols\n ("
                            + r1 + "*" + n2 + " != " + r3 + "   or   " + c1
                            + " != " + c3 + ")");

        for (int i1 = 0; i1 < r1; i1++) {
            for (int j1 = 0; j1 < c1; j1++) {
                for (int i2 = 0; i2 < n2; i2++) {
                    res.set(i1 * n2 + i2, j1, A.get(i1, j1) * B.get(i2));
                }
            }
        }
        return res;
    }

    /**
     * Returns the Kronecker product of one matrix with one vector
     * (Matrix x Vector) in dense format
     * @param A Matrix
     * @param B Vector
     * @param res Matrix such that res.numRows = A.numRows and
     *        res.numCols = A.numCols * B.size
     * @return Kronecker Product A x B
     */
    public static Matrix kroneckerMxRowVector(Matrix A, Vector B,
            Matrix res) {
        int r1 = A.numRows();
        int c1 = A.numColumns();
        int n2 = B.size();
        int r3 = res.numRows();
        int c3 = res.numColumns();
        if (r1 != r3 || c1 * n2  != c3)
            throw new IndexOutOfBoundsException(
                    "Matrices cannot be Kronecker multiplied: \n"
                            + "A.numRows != res.numRows  or  A.numCols* B.size  != res.numCols\n ("
                            + r1 + " != " + r3 + "   or   " + c1 + "*" + n2
                            + " != " + c3 + ")");

        for (int i1 = 0; i1 < r1; i1++) {
            for (int j1 = 0; j1 < c1; j1++) {
                for (int i2 = 0; i2 < n2; i2++) {
                    res.set( i1, (j1 * c1) + i2, A.get(i1, j1) * B.get(i2));
                }
            }
        }
        return res;
    }

    /**
     * Returns the Kronecker product of one vector with one matrix
     * (Vector x Matrix) in dense format
     * @param A Vector
     * @param B Matrix
     * @param res Matrix such that res.numRows = A.size * B.numRows
     *        and res.numCols = B.numCols
     * @return Kronecker Product A x B
     */
    public static Matrix kronecker(Vector A, Matrix B,
            Matrix res) {
        int n1 = A.size();
        int r2 = B.numRows();
        int c2 = B.numColumns();
        int r3 = res.numRows();
        int c3 = res.numColumns();
        if (n1 * r2 != r3 || c2 != c3)
            throw new IndexOutOfBoundsException(
                    "Matrices cannot be Kronecker multiplied: \n"
                            + "A.size * B.numRows != res.numRows   or   B.numCols != res.numCols ("
                            + n1 + " * " + r2 + " != " + r3 + "   or   " + c2
                            + " != " + c3 + ")");

        for (int i1 = 0; i1 < n1; i1++) {
            for (int i2 = 0; i2 < r2; i2++) {
                for (int j2 = 0; j2 < c2; j2++) {
                    res.set(i1 * r2 + i2, j2, A.get(i1) * B.get(i2, j2));

                }
            }
        }
        return res;
    }

    /**
     * Returns the Kronecker product of two vectors in dense format
     * @param A Vector
     * @param B Vector
     * @return Kronecker Product A x B
     */
    public static DenseVector kroneckerVectors(DenseVector A, DenseVector B) {
        int n1 = A.size();
        int n2 = B.size();
        DenseVector result = new DenseVector(n1 * n2);
        for (int i1 = 0; i1 < n1; i1++) {
            for (int i2 = 0; i2 < n2; i2++) {
                result.set(i1 * n2 + i2, A.get(i1) * B.get(i2));
            }
        }
        return result;
    }

    /**
     * Returns the Kronecker product of two vectors 
     * format
     * @param A Vector
     * @param B Vector
     * @param res Vector such that res.size = A.size * B.size
     * @return Kronecker Product A x B
     */
    public static Vector kroneckerVectors(
            Vector A, Vector B,
            Vector res) {
        int n1 = A.size();
        int n2 = B.size();
        if (n1 * n2 != res.size())
            throw new IndexOutOfBoundsException(
                    "Vectors cannot be Kronecker multiplied: (" + n1 + "*" + n2
                            + " != " + res.size() + ")");

        for (int i1 = 0; i1 < n1; i1++) {
            for (int i2 = 0; i2 < n2; i2++) {
                res.set(i1 * n2 + i2, A.get(i1) * B.get(i2));
            }
        }
        return res;
    }

    /**
     * Returns a one-row matrix in dense format with one in every
     * entry
     * @param m size of the matrix (1, m)
     * @return Row matrix with one in every entry
     */
    public static DenseMatrix OnesRow(int m) {
        DenseMatrix vec = new DenseMatrix(1, m);
        for (int i = 0; i < m; i++)
            vec.set(0, i, 1);
        return vec;
    }

    /**
     * Returns a one-column matrix in dense format with one in every
     * entry
     * @param m size of the matrix (m, 1)
     * @return One column matrix with one in every entry
     */
    public static DenseMatrix OnesCol(int m) {
        DenseMatrix vec = new DenseMatrix(m, 1);
        for (int i = 0; i < m; i++)
            vec.set(i, 0, 1);
        return vec;
    }

    /**
     * Returns a DenseVector with one in every entry
     * @param m size of the DenseVector
     * @return DenseVector with one in every entry
     */
    public static DenseVector OnesVector(int m) {
        DenseVector vec = new DenseVector(m);
        for (int i = 0; i < m; i++)
            vec.set(i, 1);
        return vec;
    }

    /**
     * Returns a Vector with one in every entry in the predefined
     * storage format
     * @param vec Vector to be modified
     * @return Vector with one in every entry
     */
    public static Vector OnesVector(
            Vector vec) {
        for (int i = 0; i < vec.size(); i++)
            vec.set(i, 1);
        return vec;
    }

    /**
     * Calculates x^n
     * @param x base
     * @param n exponent
     * @return x^n
     */
    public static double pow(double x, int n) {
        double mult = 1.0;
        if (n < 0) {
            n = -n;
            x = 1.0 / x;
        }
        for (int i = 1; i <= n; i++) {
            mult *= x;
        }
        return mult;
    }

    /**
     * Calculates the distance between two arrays, defined as the
     * maximum Euclidean distance between every entry in those arrays
     * @param v1 array
     * @param v2 array
     * @return Distance between two arrays
     */
    public static double distance(double v1[], double v2[]) {
        int n = v1.length;
        if (n != v2.length)
            return -1.0;
        double maxi = -1.0, del, del2;
        for (int i = 0; i < n; i++) {
            del2 = (v1[i] - v2[i]);
            del = (v1[i] > 0) ? del2 / v1[i] : del2;
            maxi = Math.max(maxi, del * del);
        }
        return Math.sqrt(maxi);
    }


    /**
     * Returns exp(A x) * Ones, for the value x. It uses the
     * uniformization algorithm as described in page 60 of Latouche
     * and Ramaswami
     * @param A Matrix
     * @param x evaluation point
     * @return exp(A x) * Ones
     */
    public static Matrix expTimesOnes(Matrix A, double x) {
        int n = A.numColumns();
        return exp(A, x, Matrices.identity(n), OnesRow(n));
    }

    /**
     * Returns leftMat * exp(A x) * Ones, for the value x. It uses the
     * uniformization algorithm as described in page 60 of Latouche
     * and Ramaswami
     * @param A Matrix
     * @param x evaluation point
     * @param leftMat Matrix
     * @return leftMat * exp(A x) * rightMat
     */
    public static Matrix expTimesOnes(Matrix A, double x, Matrix leftMat) {
        int n = A.numColumns();
        return exp(A, x, leftMat, OnesCol(n));
    }

    /**
     * Returns leftVec * exp(A x) * OnesVector, for the value x. It
     * uses the uniformization algorithm as described in page 60 of
     * Latouche and Ramaswami
     * @param A Matrix
     * @param x evaluation point
     * @param leftVec Vector
     * @return leftVec * exp(A x) * OnesVector
     */
    public static double expTimesOnes(Matrix A, double x,
            Vector leftVec) {
        return exp(A, x, leftVec, OnesVector(leftVec.copy()));
    }

    /**
     * Returns leftMat * exp(A x) * OnesCol, for all values x = 0 +
     * i*delta, i=0,...,n. It uses the uniformization algorithm as
     * described in page 60 of Latouche and Ramaswami
     * @param A Matrix
     * @param n number of evaluation points
     * @param delta separation between evaluation points
     * @param leftMat Matrix
     * @return leftMat * exp(A x) * OnesCol, for x = 0 + i*delta,
     *         i=0,...,n
     */
    public static Matrix[] expTimesOnes(Matrix A, int n, double delta,
            Matrix leftMat) {
        int m = A.numColumns();
        return exp(A, n, delta, leftMat, OnesCol(m), true);
    }

    /**
     * Returns leftVec * exp(A x) * OnesVector, for all values x = 0 +
     * i*delta, i=0,...,n. It uses the uniformization algorithm as
     * described in page 60 of Latouche and Ramaswami
     * @param A Matrix
     * @param n number of evaluation points
     * @param delta separation between evaluation points
     * @param leftVec Vector
     * @return leftVec * exp(A x) * OnesVector, for x = 0 + i*delta,
     *         i=0,...,n
     */
    public static double[] expTimesOnes(Matrix A, int n, double delta,
            Vector leftVec) {
        return exp(A, n, delta, leftVec, OnesVector(leftVec.copy()), true);
    }

    /**
     * Returns the value of the position (0,0) in the matrix, if its
     * number of columns is equal to one (1)
     * @param A matrix
     * @return Value of the position (0,0) in the matrix
     */
    public static double scalar(Matrix A) {
        if (A.numColumns() == 1) {
            return A.get(0, 0);
        } else {
            throw new NumberFormatException(
                    "Max iterations reached computing exp()");
        }
    }

    /**
     * Returns leftMat * exp(A x) * rightMat, for the value x. It uses
     * the uniformization algorithm as described in page 60 of
     * Latouche and Ramaswami
     * @param A Matrix
     * @param x evaluation point
     * @param leftMat Matrix
     * @param rightMat Matrix
     * @return leftMat * exp(A x) * rightMat
     */
    public static Matrix exp(Matrix A, double x, Matrix leftMat, Matrix rightMat) {
        Matrix result[] = exp(A, 2, x, leftMat, rightMat, true);
        return result[1];
    }

    /**
     * Returns leftVec * exp(A x) * rightVec, for the value x. It uses
     * the uniformization algorithm as described in page 60 of
     * Latouche and Ramaswami
     * @param A Matrix
     * @param x evaluation point
     * @param leftVec Vector
     * @param rightVec Vector
     * @return leftVec * exp(A x) * rightVec
     */
    public static double exp(Matrix A, double x,
            Vector leftVec,
            Vector rightVec) {
        double result[] = exp(A, 2, x, leftVec, rightVec, true);
        return result[1];
    }

    /**
     * Returns exp(A x), for the value x. It uses the uniformization
     * algorithm as described in page 60 of Latouche and Ramaswami
     * @param A Matrix
     * @param x evaluation point
     * @return exp(A x)
     */
    public static Matrix exp(Matrix A, double x) {
        int n = A.numColumns();
        Matrix result[] = exp(A, 2, x, Matrices.identity(n), Matrices
                .identity(n), true);
        return result[1];
    }

    /**
     * Returns leftMat * exp(A x) * rightMat, for the value x. It uses
     * the uniformization algorithm or the RungeKutta method
     * @param A Matrix
     * @param x evaluation point
     * @param leftMat Matrix
     * @param rightMat Matrix
     * @param useUniformization true if the method to use is
     *        Uniformization, false if it is RungeKutta
     * @return leftMat * exp(A x) * rightMat
     */
    public static Matrix exp(Matrix A, double x, Matrix leftMat,
            Matrix rightMat, boolean useUniformization) {
        Matrix result[] = exp(A, 2, x, leftMat, rightMat, useUniformization);
        return result[1];
    }

    /**
     * Computes leftMat * exp(A x) * rightMat, for all values x = 0 +
     * i*delta, i=0,...,n. It uses the uniformization algorithm or the
     * RungeKutta method
     * @param A Matrix
     * @param n number of evaluation points
     * @param delta separation between evaluation points
     * @param leftMat Matrix
     * @param rightMat Matrix
     * @param useUniformization true if the method to use is
     *        Uniformization, false if it is RungeKutta
     * @return leftMat * exp(A x) * rightMat, for x = 0 + i*delta,
     *         i=0,...,n
     */
    public static Matrix[] exp(Matrix A, int n, double delta, Matrix leftMat,
            Matrix rightMat, boolean useUniformization) {
        if (useUniformization) {
            return expUnif(A, n, delta, leftMat, rightMat);
        } else {
            return expRunge(A, n, delta, leftMat, rightMat);
        }
    }

    /**
     * Returns leftVec * exp(A x) * rightVec, for all values x = 0 +
     * i*delta, i=0,...,n. It uses the uniformization algorithm or the
     * RungeKutta method
     * @param A Matrix
     * @param n number of evaluation points
     * @param delta separation between evaluation points
     * @param leftVec Vector
     * @param rightVec Vector
     * @param useUniformization true if the method to use is
     *        Uniformization, false if it is RungeKutta
     * @return leftVec * exp(A x) * rightVec
     */
    public static double[] exp(Matrix A, int n, double delta,
            Vector leftVec,
            Vector rightVec, boolean useUniformization) {
        if (useUniformization) {
            return expUnif(A, n, delta, leftVec, rightVec);
        } else {
            // return expRunge(A, n, delta, leftMat, rightMat);
            return null;
        }
    }

    /**
     * Computes leftMat * exp(A x) * rightMat, for all values x = 0 +
     * i*delta, i=0,...,n. It uses the RungeKutta method
     * @param A Matrix
     * @param n number of evaluation points
     * @param delta separation between evaluation points
     * @param leftMat Matrix
     * @param rightMat Matrix
     * @return leftMat * exp(A x) * rightMat, for x = 0 + i*delta,
     *         i=0,...,n
     */
    public static Matrix[] expRunge(Matrix A, int n, double delta,
            Matrix leftMat, Matrix rightMat) {
        Matrix result[] = new Matrix[n];
        Matrix y = leftMat.copy();

        // from this point we assume x[0]=0, x[i] = i * delta equally
        // separated.
        double bigStep = delta;
        double ldaMax = -1.0;
        ldaMax = computeLdaMax(A);
        double k = Math.ceil(ldaMax * bigStep);
        double step = bigStep / k;
        double xVal = 0;
        for (int i = 0; i < n; i++) {
            // here xval = x[i]
            // result[i] = new Matrix(y.times(rightMat));
            result[i] = y.mult(rightMat, y.copy());
            for (int j = 0; j < k; j++) {
                y = runge4(A, y, step);
                xVal += step;
            }
        }
        return result;
    }

    /**
     * Executes an iteration of the RungeKutta method of order 4
     * @param A Matrix
     * @param y Matrix
     * @param step step seze for the RungeKutta method of order 4
     * @return evaluation
     */
    private static Matrix runge4(Matrix A, Matrix y, double step) {
        // double h = step / 2.0; /* the midpoint */
        Matrix t1, t2, t3, /* temporary storage arrays */
        k1, k2, k3, k4; /* for Runge-Kutta */

        // initialization
        t1 = t2 = t3 = k1 = k2 = k3 = k4 = null;
        t1 = y.add(k1 = y.multAdd(step, A, k1).scale(0.5));
        t1 = y.add(k2 = t1.multAdd(step, A, k2).scale(0.5));
        t3 = y.add(k3 = t2.mult(step, A, k3));
        k4 = t3.mult(step, A, k4);
        k2 = k2.scale(2.0);
        k1 = k1.scale(2.0);

        y = y.add(1.0 / 6.0, k1.add(k2).add(k3).add(k4));
        return y;
    }


    /**
     * Computes the maximum exit rate of the states in the Markov
     * chain
     * @param A Continuous Markov chain transition matrix
     * @return ldaMax maximum exit rate of the states in the Markov
     *         chain
     */
    private static double computeLdaMax(Matrix A) {
        double ldaMax = 0;
        int n = A.numRows();
        for (int i = 0; i < n; i++) {
            ldaMax = Math.max(ldaMax, -A.get(i, i));
        }
        return ldaMax;
    }

    /**
     * Computes leftMat * exp(A x) * rightMat, for the value x. It
     * uses the uniformization algorithm as described in page 60 of
     * Latouche and Ramaswami
     * @param A Matrix
     * @param x evaluation point
     * @param leftMat Matrix
     * @param rightMat Matrix
     * @return leftMat * exp(A x) * rightMat
     */
    public static Matrix[] expUnif(Matrix A, double x, Matrix leftMat,
            Matrix rightMat) {
        return expUnif(A, new double[] { x }, leftMat, rightMat,
                Integer.MAX_VALUE);
    }


    /**
     * Computes leftMat * exp(A x) * rightMat, for all values x = 0 +
     * i*delta, i=0,...,n. It uses the uniformization algorithm as
     * described in page 60 of Latouche and Ramaswami
     * @param A Matrix
     * @param n number of evaluation points
     * @param delta separation between evaluation points
     * @param leftMat Matrix
     * @param rightMat Matrix
     * @return leftMat * exp(A x) * rightMat
     */
    public static Matrix[] expUnif(Matrix A, int n, double delta,
            Matrix leftMat, Matrix rightMat) {
        return expUnif(A, n, delta, leftMat, rightMat, Integer.MAX_VALUE);
    }

    /**
     * Computes leftVec * exp(A x) * rightVec, for all values x = 0 +
     * i*delta, i=0,...,n. It uses the uniformization algorithm as
     * described in page 60 of Latouche and Ramaswami
     * @param A Matrix
     * @param n number of evaluation points
     * @param delta separation between evaluation points
     * @param leftVec Vector
     * @param rightVec Vector
     * @return leftVec * exp(A x) * rightVec
     */
    public static double[] expUnif(Matrix A, int n, double delta,
            Vector leftVec,
            Vector rightVec) {
        return expUnif(A, n, delta, leftVec, rightVec, Integer.MAX_VALUE);
    }

    /**
     * Computes leftMat * exp(A x) * rightMat, for all values x in
     * times. It uses the uniformization algorithm
     * @param A Matrix
     * @param times evaluation points
     * @param leftMat Matrix
     * @param rightMat Matrix
     * @return leftMat * exp(A x) * rightMat
     */
    public static Matrix[] expUnif(Matrix A, double times[], Matrix leftMat,
            Matrix rightMat) {
        return expUnif(A, times, leftMat, rightMat, Integer.MAX_VALUE);
    }

    /**
     * Computes leftMat * exp(A x) * rightMat, for all values x = 0 +
     * i*delta, i=0,...,n. It uses the uniformization algorithm
     * @param A Matrix
     * @param n num of evaluation points
     * @param delta evaluation points separation
     * @param leftMat Matrix
     * @param rightMat Matrix
     * @param truncate upper bound for iterations
     * @return leftMat * exp(A x) * rightMat
     */
    public static Matrix[] expUnif(Matrix A, int n, double delta,
            Matrix leftMat, Matrix rightMat, int truncate) {
        double times[] = new double[n];
        for (int i = 0; i < n; i++)
            times[i] = delta * i;
        return expUnif(A, times, leftMat, rightMat, truncate);
    }

    /**
     * Computes leftVec * exp(A x) * rightVec, for all values x = 0 +
     * i*delta, i=0,...,n. It uses the uniformization algorithm
     * @param A Matrix
     * @param n num of evaluation points
     * @param delta evaluation points separation
     * @param leftVec Vector
     * @param rightVec Vector
     * @param truncate upper bound for iterations
     * @return leftVec * exp(A x) * rightVec
     */
    public static double[] expUnif(Matrix A, int n, double delta,
            Vector leftVec,
            Vector rightVec, int truncate) {
        double times[] = new double[n];
        for (int i = 0; i < n; i++)
            times[i] = delta * i;
        return expUnif(A, times, leftVec, rightVec, truncate);
    }


    /**
     * Computes leftMat * exp(A x) * rightMat, for all values x in
     * times. It uses the uniformization algorithm
     * @param A Matrix
     * @param times evaluation points
     * @param leftMat Matrix
     * @param rightMat Matrix
     * @param truncate upper bound for iterations
     * @return leftMat * exp(A x) * rightMat
     */
    public static Matrix[] expUnif(Matrix A, double times[], Matrix leftMat,
            Matrix rightMat, int truncate) {

        List<Matrix> Avec = new ArrayList<Matrix>();
        double dif;
        double ldaMax = computeLdaMax(A);
        Matrix matP = getNormalized(A, ldaMax);
        int n = times.length;
        long MaxIterations = Math.max(200, (int) (4 * ldaMax * times[n - 1]));
        // TODO: change for a stochastic P

        Matrix M = new DenseMatrix(leftMat.numRows(), rightMat.numColumns());
        Matrix M1 = new DenseMatrix(leftMat.numRows(), A.numColumns());
        // Note: M = L (I-P)^-1 R
        Matrix temp = (Matrices.identity(A.numRows()).add(-1, matP));
        Matrix tempInv = new DenseMatrix(A.numRows(), A.numColumns());
        tempInv = temp.solve(Matrices.identity(temp.numRows()), tempInv);
        M1 = leftMat.mult(tempInv, M1);
        M = M1.mult(rightMat, M);

        Matrix Ak = new DenseMatrix(leftMat.numRows(), rightMat.numColumns());
        Ak = leftMat.mult(rightMat, Ak);
        Avec.add(Ak.copy());
        DenseMatrix sumA = (DenseMatrix) Ak.copy();
        int k = 0;

        Matrix V = new DenseMatrix(matP.numRows(), rightMat.numColumns());
        V = matP.mult(rightMat, V);
        do {
            k++;
            Ak = leftMat.mult(V, Ak); // Ak = leftMat
            Avec.add(Ak.copy());
            sumA = (DenseMatrix) sumA.add(Ak);
            V = matP.mult(V, new DenseMatrix(matP.numRows(), rightMat
                    .numColumns())); // Vk = P * RightMat
            Matrix difMat = new DenseMatrix(sumA.numRows(), sumA.numColumns());
            difMat = sumA.copy();
            difMat = difMat.add(-1, M);
            dif = Math.abs(difMat.norm(Matrix.Norm.Infinity));
        } while ((dif > Epsilon) && (k < MaxIterations) && (k < truncate));

        if ((k == MaxIterations) && (truncate == Integer.MAX_VALUE)) {
            System.err
                    .println("Max iterations reached (" + MaxIterations + ")");
        }

        int maxK = k;
        Matrix As[] = new Matrix[maxK + 1];
        for (k = 0; k <= maxK; k++) {
            As[k] = Avec.get(k);
        }
        Matrix[] result = new Matrix[n];
        for (int i = 0; i < n; i++) {
            double ldaX = ldaMax * times[i];
            double pk = Math.exp(-ldaX);

            if (pk < Double.MIN_VALUE) {
                // if (true){
                System.out.println("pk menor que Double.MIN_VALUE");
                result[i] = resultFromMedian(As, ldaX);
            } else {
                double sumPk = pk;
                result[i] = As[0].copy();
                result[i].scale(pk);

                for (k = 1; k <= maxK; k++) {
                    pk = pk * (ldaX / k);
                    sumPk += pk;
                    result[i] = result[i].add(pk, As[k]);
                }// next k
            }// else
        }// next i
        Matrix GrmResult[] = new Matrix[n];
        for (int i = 0; i < n; i++)
            GrmResult[i] = result[i].copy();
        return GrmResult;
    }

    /**
     * Computes leftVec * exp(A x) * rightVec, for all values x in
     * times. It uses the uniformization algorithm
     * @param A Matrix
     * @param times evaluation points
     * @param leftVec Vector
     * @param rightVec Vector
     * @param truncate upper bound for iterations
     * @return leftVec * exp(A x) * rightVec
     */
    public static double[] expUnif(Matrix A, double times[],
            Vector leftVec,
            Vector rightVec, int truncate) {

        List<Double> Avec = new ArrayList<Double>();

        double dif;
        double ldaMax = computeLdaMax(A);
        Matrix matP = getNormalized(A, ldaMax);
        int n = times.length;
        long MaxIterations = Math.max(200, (int) (4 * ldaMax * times[n - 1]));
        // TODO: change for a stochastic P 
        double m = 0;
        Vector M1 = rightVec.copy();
        //TODO: extend for disperse matrices 
        // Note: M = L (I-P)^-1 R
        Matrix temp = matP.copy().scale(-1).add(Matrices.identity(A.numRows()));
        // M1 = (I-P)^-1 * rightVec
        IterativeSolver solver = new BiCG(M1);
        try {
            solver.solve(temp,rightVec, M1);
        } catch (IterativeSolverNotConvergedException e) {
            e.printStackTrace();
        }
        m = leftVec.dot(M1);

        double ak = 0;
        ak = leftVec.dot(rightVec);
        Avec.add(new Double(ak));

        double sumA = ak;
        int k = 0;
        Vector V = leftVec.copy().zero();
        V = matP.mult(rightVec, V);
        do {
            k++;
            ak = leftVec.dot(V); // Ak = leftMat *(P^k)*RightMat
            Avec.add(new Double(ak));
            sumA += ak;
            V = matP.mult(V, V.copy()); // Vk = P * RightMat

            double difMat = sumA;
            difMat = difMat - m;
            dif = Math.abs(difMat);

        } while ((dif > Epsilon) && (k < MaxIterations) && (k < truncate));

        if ((k == MaxIterations) && (truncate == Integer.MAX_VALUE)) {
            //TODO 
           // throw new RuntimeException("Max iterations reached (" + MaxIterations + ")");
            // return expRunge(n, delta, leftMat, rightMat);
        }

        int maxK = k;
        double As[] = new double[maxK + 1];
        for (k = 0; k <= maxK; k++) {
            As[k] = Avec.get(k).doubleValue();
        }
        double[] result = new double[n];
        for (int i = 0; i < n; i++) {
            double ldaX = ldaMax * times[i];
            double pk = Math.exp(-ldaX);

            if (pk < Double.MIN_VALUE) {
                System.out.println("pk menor que Double.MIN_VALUE");
                result[i] = resultFromMedian(As, ldaX);
            } else {
                double sumPk = pk;
                result[i] = As[0];
                result[i] *= pk;

                for (k = 1; k <= maxK; k++) {
                    pk = pk * (ldaX / k);
                    sumPk += pk;
                    result[i] = result[i] + pk * As[k];
                }// next k
            }// else
        }// next i
        double GrmResult[] = new double[n];
        for (int i = 0; i < n; i++)
            GrmResult[i] = result[i];
        return GrmResult;
    }

    /**
     * Computes the value of exp(A[i] x) when it is to small
     * @param A Matrix
     * @param ldaX
     * @return exp(A[i] x)
     */
    private static Matrix resultFromMedian(Matrix A[], double ldaX) {
        if (ldaX == 0)
            return A[0];
        int median = (int) Math.floor(ldaX);
        double p_plus = Math.exp(-ldaX + median * Math.log(ldaX)
                - lnFactorial(median));
        double p_minus = p_plus;
        double sumPk = p_plus;
        Matrix result = (median < A.length) ? A[median].scale(p_plus)
                : A[median].zero();
        int k_minus = median;
        int k_plus = median;
        int maxK = Math.max(median, A.length - median + 1);
        for (int k = 1; ((k < maxK) && (sumPk < 0.999)); k++) {
            k_minus--;
            k_plus++;
            p_plus = p_plus * (ldaX / k_plus);
            p_minus = p_minus * (k_minus + 1) / ldaX;
            sumPk += p_plus + p_minus;
            if (k_plus < A.length)
                result = result.add(p_plus, A[k_plus]);
            if ((k_minus < A.length) && (k_minus >= 0))
                result = result.add(p_minus, A[k_minus]);
        }// next k
        return result;
    }


    /**
     * Computes the value of exp(A[i] x) when it is to small
     * @param A double
     * @param ldaX
     * @return exp(A[i] x)
     */
    private static double resultFromMedian(double A[], double ldaX) {
        if (ldaX == 0)
            return A[0];
        int median = (int) Math.floor(ldaX);
        double p_plus = Math.exp(-ldaX + median * Math.log(ldaX)
                - lnFactorial(median));
        double p_minus = p_plus;
        double sumPk = p_plus;
        double result = (median < A.length) ? A[median] * p_plus : 0;
        int k_minus = median;
        int k_plus = median;
        int maxK = Math.max(median, A.length - median + 1);
        for (int k = 1; ((k < maxK) && (sumPk < 0.999)); k++) {
            k_minus--;
            k_plus++;
            p_plus = p_plus * (ldaX / k_plus);
            p_minus = p_minus * (k_minus + 1) / ldaX;
            sumPk += p_plus + p_minus;
            if (k_plus < A.length)
                result = result + p_plus * A[k_plus];
            if ((k_minus < A.length) && (k_minus >= 0))
                result = result + p_minus * A[k_minus];
        }// next k
        return result;
    }

    /**
     * Computes the normalized matrix for the uniformization algorithm
     * @param A Matrix
     * @return Normalized matrix
     */
    private static Matrix getNormalized(Matrix A) {
        Matrix normal = null;
        double ldaMax = computeLdaMax(A);
        normal = getNormalized(A, ldaMax);
        return normal;
    }

    /**
     * Computes the normalized matrix for the uniformization
     * algorithm, using the given value of Lamda.
     * @param A Matrix
     * @param ldaMax Maximum lambda
     * @return Normalized matrix
     */
    private static Matrix getNormalized(Matrix A, double ldaMax) {
        Matrix normal = null;
        int n = A.numRows();
        normal = A.copy().scale(1.0/ldaMax).add(Matrices.identity(n));
        return normal;
    }

    /**
     * Computes data average
     * @param datos
     * @return Data Average
     */
    public static double average(double[] datos) {
        double media = 0;
        for (int i = 0; i < datos.length; i++)
            media += datos[i];
        return media / datos.length;
    }

    /**
     * Computes the second moment of the data
     * @param data data trace
     * @return Second data trace moment
     */
    public static double average2(double[] data) {
        double media2 = 0;
        for (int i = 0; i < data.length; i++)
            media2 += data[i] * data[i];
        return media2 / data.length;
    }

    /**
     * Return the variance of the data
     * @param data data trace
     * @return Data Variance
     */
    public static double variance(double[] data) {
        double media = average(data);
        return average2(data) - media * media;
    }

    /**
     * Return the Coefficient of Variation of the data trace
     * @param data data trace
     * @return Data Coefficient of Variation
     */
    public static double CV(double[] data) {
        double media = average(data);
        return variance(data) / (media * media);
    }

    /**
     * Concatenates the columns of the matrices, keeping the same
     * number of rows in dense format
     * @param A DenseMatrix
     * @param B DenseMatrix
     * @return Resulting matrix from cancatenation
     */
    public static DenseMatrix concatCols(DenseMatrix A, DenseMatrix B) {
        int r1 = A.numRows();
        int c1 = A.numColumns();
        int r2 = B.numRows();
        int c2 = B.numColumns();
        if (r1 != r2)
            throw new IndexOutOfBoundsException(
                    "Matrices cannot be concatenated: (" + r1 + "," + c1
                            + ")|(" + r2 + "," + c2 + ")");
        DenseMatrix C = new DenseMatrix(r1, c1 + c2);

        for (int i = 0; i < r1; i++) {
            for (int j = 0; j < c1; j++) {
                C.set(i, j, A.get(i, j));
            }
        }
        // C.setMatrix(0, r1 - 1, 0, c1 - 1, A);
        for (int i = 0; i < r1; i++) {
            for (int j = c1; j < c1 + c2; j++) {
                C.set(i, j, B.get(i, j - c1));
            }
        }
        // C.setMatrix(0, r1 - 1, c1, c1 + c2 - 1, B);
        return C;
    }

    /**
     * Concatenates the columns of the matrices, keeping the same
     * number of rows in the predefined format
     * @param A Matrix
     * @param B Matrix
     * @param res Resulting matrix
     * @return Resulting matrix from cancatenation
     */
    public static Matrix concatCols(Matrix A, Matrix B, Matrix res) {
        int r1 = A.numRows();
        int c1 = A.numColumns();
        int r2 = B.numRows();
        int c2 = B.numColumns();
        int r3 = res.numRows();
        int c3 = res.numColumns();

        if (r1 != r2 || r1 != r3 || c1 + c2 != c3)
            throw new IndexOutOfBoundsException(
                    "Matrices cannot be concatenated: (" + r1 + "," + c1
                            + ")|(" + r2 + "," + c2 + ") in (" + r3 + "," + c3
                            + ")");

        for (int i = 0; i < r1; i++) {
            for (int j = 0; j < c1; j++) {
                res.set(i, j, A.get(i, j));
            }
        }
        for (int i = 0; i < r1; i++) {
            for (int j = c1; j < c1 + c2; j++) {
                res.set(i, j, B.get(i, j - c1));
            }
        }
        return res;
    }

    /**
     * Concatenates the rows of the matrices, keeping the same number
     * of colums in dense format
     * @param A DenseMatrix
     * @param B DenseMatrix
     * @return Resulting matrix from cancatenation
     */
    public static DenseMatrix concatRows(DenseMatrix A, DenseMatrix B) {
        int r1 = A.numRows();
        int c1 = A.numColumns();
        int r2 = B.numRows();
        int c2 = B.numColumns();
        if (c1 != c2)
            throw new IndexOutOfBoundsException(
                    "Matrices cannot be concatenated: (" + r1 + "," + c1
                            + ")\n--\n(" + r2 + "," + c2 + ")");
        DenseMatrix C = new DenseMatrix(r1 + r2, c1);
        for (int i = 0; i < r1; i++) {
            for (int j = 0; j < c1; j++) {
                C.set(i, j, A.get(i, j));
            }
        }
        for (int i = r1; i < r1 + r2; i++) {
            for (int j = 0; j < c1; j++) {
                C.set(i, j, B.get(i - r1, j));
            }
        }
        return C;
    }

    /**
     * Concatenates the rows of the matrices, keeping the same number
     * of colums in the predefined format
     * @param A Matrix
     * @param B Matrix
     * @param res Resulting matrix
     * @return Resulting matrix from cancatenation
     */
    public static Matrix concatRows(Matrix A, Matrix B, Matrix res) {
        int r1 = A.numRows();
        int c1 = A.numColumns();
        int r2 = B.numRows();
        int c2 = B.numColumns();
        int r3 = res.numRows();
        int c3 = res.numColumns();
        if (c1 != c2 || c1 != c3 || r1 + r2 != r3)
            throw new IndexOutOfBoundsException(
                    "Matrices cannot be concatenated: (" + r1 + "," + c1
                            + ")\n--\n(" + r2 + "," + c2 + ")in (" + r3 + ","
                            + c3 + ")");

        for (int i = 0; i < r1; i++) {
            for (int j = 0; j < c1; j++) {
                res.set(i, j, A.get(i, j));
            }
        }
        for (int i = r1; i < r1 + r2; i++) {
            for (int j = 0; j < c1; j++) {
                res.set(i, j, B.get(i - r1, j));
            }
        }
        return res;
    }

    /**
     * Concatenates the colums of the left and right upper matrices
     * and the result is concatenated by rows with the concatenation
     * of left and right lower matrices
     * @param leftUp Left upper Matrix
     * @param rightUp Right upper Matrix
     * @param leftDown Left lower Matrix
     * @param rightDown Right lower matrix
     * @param res Resulting Matrix
     * @return resulting matrix from concatenation
     */
    public static Matrix concatQuad(Matrix leftUp, Matrix rightUp,
            Matrix leftDown, Matrix rightDown, Matrix res) {

        int r1 = leftUp.numRows();
        int c1 = leftUp.numColumns();
        int r2 = rightUp.numRows();
        int c2 = rightUp.numColumns();
        int r3 = leftDown.numRows();
        int c3 = leftDown.numColumns();
        int r4 = rightDown.numRows();
        int c4 = rightDown.numColumns();
        int r5 = res.numRows();
        int c5 = res.numColumns();

        if (r1 != r2 || r3 != r4 || r1 + r3 != r5 || c1 != c3 || c2 != c4
                || c1 + c2 != c5)
            throw new IndexOutOfBoundsException(
                    "Matrices cannot be concatenated: (" + r1 + "," + c1
                            + ")|(" + r2 + "," + c2 + ")\n--\n" + "" + r3 + ","
                            + c3 + ")|(" + r4 + "," + c4 + ") in (" + r5 + ","
                            + c5 + ")");

        for (int i = 0; i < r1; i++) {
            for (int j = 0; j < c1; j++)
                res.set(i, j, leftUp.get(i, j));
            for (int j = 0; j < c2; j++)
                res.set(i, c1 + j, rightUp.get(i, j));

        }

        for (int i = r1; i < r1 + r3; i++) {
            for (int j = 0; j < c1; j++)
                res.set(i, j, leftDown.get(i - r1, j));
            for (int j = 0; j < c2; j++)
                res.set(i, c1 + j, rightDown.get(i - r1, j));

        }
        return res;
    }

    /**
     * Concatenates the vectors in dense format
     * @param A DenseVector
     * @param B DenseVector
     * @return Resulting DenseVector from cancatenation
     */
    public static DenseVector concatVectors(DenseVector A, DenseVector B) {
        int c1 = A.size();
        int c2 = B.size();

        DenseVector C = new DenseVector(c1 + c2);

        for (int i = 0; i < c1; i++) {
            C.set(i, A.get(i));
        }

        for (int i = c1; i < c1 + c2; i++) {
            C.set(i, B.get(i - c1));
        }
        return C;
    }

    /**
     * Concatenates the vectors in the predefined format
     * @param A Vector
     * @param B vector
     * @param res Resulting vector
     * @return Resulting vectormatrix from cancatenation
     */
    public static Vector concatVectors(
            Vector A, Vector B,
            Vector res) {
        int n1 = A.size();
        int n2 = B.size();

        if (res.size() != n1 + n2)
            throw new IndexOutOfBoundsException(
                    "res.size() != A.size() + B.size() (" + res.size() + " != "
                            + (A.size() + B.size()) + ")");

        for (int i = 0; i < n1; i++) {
            res.set(i, A.get(i));
        }

        for (int i = n1; i < n1 + n2; i++) {
            res.set(i, B.get(i - n1));
        }

        return res;
    }

    /**
     * Computes the producto of two vectors A x B^T
     * @param A Vector
     * @param B Vector
     * @param res Vector to store the resulting matrix
     * @return res = A x B^T
     */
    public static Matrix multVector(Vector A,
            Vector B, Matrix res) {
        int n = res.numRows();
        int m = res.numColumns();

        if (A.size() != n)
            throw new IndexOutOfBoundsException("A.size() != res.numRows() ("
                    + A.size() + " != " + n + ")");
        if (B.size() != m)
            throw new IndexOutOfBoundsException(
                    "B.size() != res.numColumns() (" + B.size() + " != " + m
                            + ")");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                res.set(i, j, A.get(i) * B.get(j));
            }
        }
        return res;
    }

    /**
     * Computes k power of the matrix A
     * @param A matrix base
     * @param k exponent
     * @return A^k
     */
    public static Matrix matPower(Matrix A, int k) {
        int n = A.numColumns();
        if (k < 0)
            throw new IllegalArgumentException(
                    "The moments are defined for k >= 0");
        else if (k == 0)
            return Matrices.identity(n);
        else {
            Matrix result = A.copy();
            for (int i = 1; i < k; i++) {
                result = result.mult(A, A.copy());
            }
            return result;
        }
    }

    /**
     * Computes the kth power of the matrix A premultiplied by leftVec
     * and postmultiplied by rightVec
     * @param A Matrix base
     * @param k exponent
     * @param leftVec Vector
     * @param rightVec Vector
     * @return leftVec * A^k * rightVec
     */
    public static double matPower(Matrix A, int k,
            Vector leftVec,
            Vector rightVec) {
        if (k == 0)
            return leftVec.dot(rightVec);
        else if (k == 1)
            return leftVec.dot(A.mult(rightVec, rightVec.copy().zero()));
        else
            return leftVec.dot((matPower(A, k)).mult(rightVec, rightVec.copy()
                    .zero()));

    }

    /**
     * Computes the sum of the first k terms of the succesion T^(j-1),
     * from j = 1
     * @param A Matrix base
     * @param k maximum exponent
     * @param leftVec
     * @param rightVec
     * @return leftVec * sum_(j=1)^k T^(j-1) * rightVec
     */
    public static double sumMatPower(Matrix A, int k,
            Vector leftVec,
            Vector rightVec) {
        if (k < 1) {
            System.out.println("Sum of Matrix Power less than one");
            return 0;
        } else {
            Matrix temp = Matrices.identity(A.numRows());
            Matrix sum = temp.copy();
            for (int i = 1; i < k; i++) {
                temp.mult(A, temp.copy());
                sum.add(temp);
            }
            return leftVec.dot((sum.mult(rightVec, rightVec.copy().zero())));
        }

    }
    
    /**
     * Checks whether a vector is sub-stochastic or not 
     * @param a Vector to check
     * @return true if the vector is sub-stochastic, false otherwise
     */
    public static boolean checkSubStochasticVector(Vector a){
    	boolean res = false; 
    	boolean nonneg = true; 
    	double cumsum = a.get(0); 
    	for (int i = 1; i < a.size(); i++){
    		if (a.get(i) < -Epsilon){
    			nonneg = false;
    			break;
    		}
    		cumsum += a.get(i);
    	}
    	if (cumsum <= 1 + Epsilon && nonneg) res = true;  
    	return res;
    }
    
    /**
     * Checks whether a matrix is a sub-generator or not 
     * @param a Vector to check
     * @return true if the vector is sub-stochastic, false otherwise
     */
    public static boolean checkSubGeneratorMatrix(Matrix A){
    	boolean res = true; 
    	double n = A.numRows();
    	double m = A.numColumns();
    	if (n == m){
	    	for (int i = 0; i < n; i++){
	    		double cumsum = A.get(i,0); 
	        	for (int j = 1; j < m; j++){
	        		cumsum += A.get(i,j);
		        	if ((i == j && A.get(i,j) > - Epsilon) || (i != j && A.get(i,j) < - Epsilon)){
		    			res = false;
		    			i += n;
		    			break;
		    		}
	        	}
	        	if (cumsum > Epsilon){
	        		res = false;
	        		break;
	        	}
	    	}
	    }
    	return res;
    }

}
