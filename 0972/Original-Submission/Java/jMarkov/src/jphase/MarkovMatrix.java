package jphase;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;

import Jama.Matrix;



/**
 * This class extends Jama's Matrix to include new operations like the Kronecker product
 * 
 * @author German Riano. Universidad de los Andes. (C) 2007
 *
 */
public class MarkovMatrix extends Matrix {
    private static final long serialVersionUID = 1969;

    static double Epsilon = 1.0E-10;// Precision for computations
    public static boolean useUniformization = true;

    public MarkovMatrix(double[][] mat) {
        super(mat);
    }

    public MarkovMatrix(Matrix mat) {
        super(mat.getArray());
    }

    public static MarkovMatrix readTxt(String stg) {
        String rowStg;
        StringTokenizer tkn = new StringTokenizer(stg, ":;", false);
        // first token are row and col dimensions
        int m = Integer.valueOf(tkn.nextToken()).intValue();
        int n = Integer.valueOf(tkn.nextToken()).intValue();
        double A[][] = new double[m][n];
        int i = 0, j = 0;
        while (tkn.hasMoreTokens()) {
            rowStg = tkn.nextToken();
            StringTokenizer tkn2 = new StringTokenizer(rowStg, " ", false);
            j = 0;
            while (tkn2.hasMoreTokens()) {
                A[i][j] = Double.valueOf(tkn2.nextToken()).doubleValue();
                j++;
            }
            i++;
        }
        return new MarkovMatrix(A);
    }

    static MarkovMatrix Ones(int m, int n) {
        return (MarkovMatrix) (new Matrix(m, n, 1.0));
    }

    static MarkovMatrix Ones(int m) {
        return new MarkovMatrix(new Matrix(m, 1, 1.0));
    }

    public static MarkovMatrix Zeros(int rows, int cols) {
        return new MarkovMatrix(new Matrix(rows, cols));
    }

    // ****************
    // Cover Functions
    // ****************

    @Override
    public Matrix inverse() {
        return (new MarkovMatrix(super.inverse()));
    }

    @Override
    public Matrix uminus() {
        return new MarkovMatrix(super.uminus());
    }

    public static Matrix identity(int n) {
        return new MarkovMatrix(Matrix.identity(n, n));
    }

    @Override
    public Matrix times(double s) {
        return new MarkovMatrix(super.times(s));
    }

    public int size() {
        return super.getColumnDimension();
    }

    @Override
    public Matrix times(Matrix B) {
        return new MarkovMatrix(super.times(B));
    }

    public static MarkovMatrix toMarkovMatrix(Matrix A) {
        return new MarkovMatrix(A);
    }

    /**
     * Solves X*A = B, actually solved as A'*X' = B'
     * @param B right hand side
     * @return solution if A is square, least squares solution
     *         otherwise.
     */

    @Override
    public Matrix solveTranspose(Matrix B) {
        return transpose().solve(B.transpose()).transpose();
    }

    /**
     * Determines if the matrix is stochastic.
     * @return true if the matrix is stochastic.
     */
    public boolean isStochastic() {
        boolean result = true;
        int n = size();
        for (int i = 0; result && i < n; i++) {
            double sum = 0.0;
            for (int j = 0; result && j < n; j++) {
                double val = get(i, j);
                result &= (val >= 0.0);
                sum += get(i, j);
            }
            result &= Math.abs(sum - 1.0) < Epsilon;
        }
        return result;
    }

    /**
     * @return k-th power of A
     */
    public MarkovMatrix power(int k) {
        int n = this.getRowDimension();
        Matrix result = super.identity(n, n);
        for (int i = 1; i <= k; i++) {
            result = result.times(this);
        }
        return (new MarkovMatrix(result));
    }

    /**
     * 
     * @return This matrix times a vector of ones
     */
    public MarkovMatrix timesOne() {
        int m = this.getRowDimension();
        int n = this.getColumnDimension();
        double sum = 0.0;
        double result[][] = new double[m][1];
        for (int i = 0; i < m; i++) {
            sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += this.get(i, j);
            }
            result[i][0] = sum;
        }
        return new MarkovMatrix(result);
    }

    public double scalar() {
        if (size() == 1) {
            return super.get(0, 0);
        } else {
            throw new NumberFormatException(
                    "Max iterations reached computing exp()");
        }
    }

    public MarkovMatrix compExp() {
        int n = getRowDimension();
        int m = getColumnDimension();
        double A[][] = getArrayCopy();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                A[i][j] = Math.exp(A[i][j]);
            }
        }
        return new MarkovMatrix(A);
    }

    /**
     * Computes the entry-wise logarithm of this matrix
     * @return The entry-wise logarithm of this matrix
     */
    public MarkovMatrix compLog() {
        int n = getRowDimension();
        int m = getColumnDimension();
        double A[][] = getArrayCopy();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                A[i][j] = Math.log(Math.abs(A[i][j]));
            }
        }
        return new MarkovMatrix(A);
    }

    public MarkovMatrix plus(double x) {
        int n = getRowDimension();
        int m = getColumnDimension();
        double A[][] = getArrayCopy();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                A[i][j] = (A[i][j] + x);
            }
        }
        return new MarkovMatrix(A);
    }

    @SuppressWarnings("unusedPrivate")
    private MarkovMatrix compTimesExp(double x) {
        int n = getRowDimension();
        int m = getColumnDimension();
        double A[][] = getArrayCopy();
        double sign;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (A[i][j] != 0) {
                    sign = (A[i][j] > 0) ? 1.0 : -1.0;
                    A[i][j] = sign * Math.exp(Math.log(Math.abs(A[i][j])) + x);
                }
            }
        }
        return new MarkovMatrix(A);
    }

    private List<Matrix> powers;
    private int maxPower = 1;

    public Matrix pow(int k) {
        if (powers == null) {
            powers = new ArrayList<Matrix>(100);
            powers.add(identity(getRowDimension()));
            powers.add(this);
        }
        if (k <= maxPower) {
            return powers.get(k);
        } else {
            Matrix newPower = null;
            for (int i = maxPower + 1; i <= k; i++) {
                newPower = powers.get(i - 1);
                newPower = newPower.times(this);
                powers.add(newPower);
            }
            maxPower = k;
            return newPower;
        }
    }

    private MarkovMatrix normal = null;
    private double ldaMax = -1.0;

    MarkovMatrix getNormalized() {
        if (normal == null) {
            int n = this.size();
            computeLdaMax();
            normal = new MarkovMatrix(identity(n, n).plus(
                    this.times(1.0 / ldaMax)));
        }
        return normal;
    }

    private void computeLdaMax() {
        int n = this.size();
        for (int i = 0; i < n; i++) {
            ldaMax = Math.max(ldaMax, -this.get(i, i));
        }
    }

    public MarkovMatrix expTimesOnes(double x) {
        int n = getColumnDimension();
        return exp(x, identity(n), Ones(n));
    }

    public MarkovMatrix expTimesOnes(double x, Matrix leftMatrix) {
        int n = getColumnDimension();
        return exp(x, leftMatrix, Ones(n));
    }

    public MarkovMatrix[] expTimesOnes(int N, double delta, Matrix leftMatrix) {
        int n = getColumnDimension();
        return exp(N, delta, leftMatrix, Ones(n));
    }

    public MarkovMatrix exp(double x, Matrix leftMat, Matrix rightMat) {
        MarkovMatrix result[] = exp(2, x, leftMat, rightMat);
        return result[1];
    }

    public MarkovMatrix exp(double x) {
        int n = this.getRowDimension();
        MarkovMatrix result[] = exp(2, x, identity(n), identity(n));
        return result[1];
    }

    public MarkovMatrix[] exp(int n, double delta, Matrix leftMat,
            Matrix rightMat) {
        if (useUniformization) {
            return expUnif(n, delta, leftMat, rightMat);
        } else {
            return expRunge(n, delta, leftMat, rightMat);
        }
    }
    
    public MarkovMatrix[] expRunge(int n, double delta, Matrix leftMat,
            Matrix rightMat) {
        MarkovMatrix result[] = new MarkovMatrix[n];
        Matrix y = leftMat;

        // from this point we assume x[0]=0, x[i] = i * delta equally
        // separated.
        double bigStep = delta;
        computeLdaMax();
        double k = Math.ceil(ldaMax * bigStep);
        double step = bigStep / k;
        double xVal = 0;
        for (int i = 0; i < n; i++) {
            // here xval = x[i]
            result[i] = new MarkovMatrix(y.times(rightMat));
            for (int j = 0; j < k; j++) {
                y = runge4(y, step);
                xVal += step;
            }
        }
        return result;
    }

    /**
     * Computes leftMat * exp(A x) * rightMat, for the value x. It
     * uses the uniformization algorithm as described in page 60 of
     * Latouche and Ramaswami
     */

    public MarkovMatrix[] expUnif(double x, Matrix leftMat, Matrix rightMat) {
        return expUnif(new double[] { x }, leftMat, rightMat, Integer.MAX_VALUE);
    }

    /**
     * Computes leftMat * exp(A x) * rightMat, for all values x. It
     * uses the uniformization algorithm as described in page 60 of
     * Latouche and Ramaswami
     */

    public MarkovMatrix[] expUnif(int n, double delta, Matrix leftMat,
            Matrix rightMat) {
        return expUnif(n, delta, leftMat, rightMat, Integer.MAX_VALUE);
    }

    /**
     * Computes leftMat * exp(A x) * rightMat, for all values x. It
     * uses the uniformization algorithm as described in page 60 of
     * Latouche and Ramaswami
     */
    public MarkovMatrix[] expUnif(double times[], Matrix leftMat,
            Matrix rightMat) {
        return expUnif(times, leftMat, rightMat, Integer.MAX_VALUE);
    }

    /**
     * Computes leftMat * exp(A x) * rightMat, for all values x= 0,
     * delta, 2delta, 3delta,.... It uses the uniformization algorithm
     * as described in page 60 of Latouche and Ramaswami
     */
    public MarkovMatrix[] expUnif(int n, double delta, Matrix leftMat,
            Matrix rightMat, int truncate) {
        double times[] = new double[n];
        for (int i = 0; i < n; i++)
            times[i] = delta * i;
        return expUnif(times, leftMat, rightMat, truncate);

    }

    /**
     * Computes leftMat * exp(A x) * rightMat, for all values x. It
     * uses the uniformization algorithm as described in page 60 of
     * Latouche and Ramaswami
     * @return leftMat * exp(A x) * rightMat, for all values x
     */
    public MarkovMatrix[] expUnif(double times[], Matrix leftMat,
            Matrix rightMat, int truncate) {

        Vector<Matrix> Avec = new Vector<Matrix>();
        double dif = Double.MAX_VALUE;
        int maxK = Integer.MAX_VALUE;
        boolean isStoch = false;
        Matrix M = null;

        MarkovMatrix matP = getNormalized();
        int n = times.length;
        int MaxIterations = Math.max(200, (int) (3 * ldaMax * times[n - 1]));
        if (!matP.isStochastic()) {
            // Note: M = L (I-P)^-1 R
            Matrix sol = (identity(size()).minus(matP)).solve(rightMat);
            M = leftMat.times(sol);
        } else {
            isStoch = true;
            maxK = MaxIterations - 1;
        }
        Matrix Ak = leftMat.times(rightMat);
        Avec.addElement(Ak);
        Matrix sumA = Ak.copy();
        int k = 0;
        Matrix V = matP.times(rightMat);
        do {
            k++;
            Ak = leftMat.times(V); // Ak = leftMat *(P^k)*RightMat
            Avec.addElement(Ak);
            sumA = sumA.plus(Ak);
            V = matP.times(V); // Vk = P * RightMat
            if (!isStoch)
                dif = Math.abs(sumA.minus(M).normInf());
        } while ((k < maxK) && (dif > Epsilon) && (k < MaxIterations)
                && (k < truncate));
        if ((k == MaxIterations) && (truncate == Integer.MAX_VALUE)) {
            throw new RuntimeException("Max iterations reached ("
                    + MaxIterations + ")");
        }

        if (!isStoch) {
            maxK = k;
        }

        Matrix A[] = new Matrix[maxK + 1];
        for (k = 0; k <= maxK; k++) {
            A[k] = Avec.elementAt(k);
        }
        Matrix[] result = new Matrix[n];
        for (int i = 0; i < n; i++) {
            double ldaX = ldaMax * times[i];
            double pk = Math.exp(-ldaX);

            if (pk < Double.MIN_VALUE) {
                // if (true){
                result[i] = resultFromMedian(A, ldaX);
            } else {
                double sumPk = pk;
                result[i] = A[0].times(pk);
                for (k = 1; k <= maxK; k++) {
                    pk = pk * (ldaX / k);
                    sumPk += pk;
                    result[i] = result[i].plus(A[k].times(pk));
                    if (1.0 - sumPk < Epsilon) break;
                }// next k
            }// else
        }// next i
        MarkovMatrix GrmResult[] = new MarkovMatrix[n];
        for (int i = 0; i < n; i++)
            GrmResult[i] = new MarkovMatrix(result[i]);
        return GrmResult;
    }

    private Matrix resultFromMedian(Matrix A[], double ldaX) {

        if (ldaX == 0)
            return A[0];
        int median = (int) Math.floor(ldaX);
        double p_plus = Math.exp(-ldaX + median * Math.log(ldaX)
                - Utils.lnFactorial(median));
        double p_minus = p_plus;
        double sumPk = p_plus;
        int m = A[0].getColumnDimension();
        Matrix result = (median < A.length) ? A[median].times(p_plus) : Zeros(
                m, m);
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
                result = result.plus(A[k_plus].times(p_plus));
            if ((k_minus < A.length) && (k_minus >= 0))
                result = result.plus(A[k_minus].times(p_minus));
        }// next k

        return result;
    }
    
    /**
     * Computes the kronecker sum of this matrix and B
     * @param B matrix to sum to this
     * @return Kronecker sum of this matrix and B
     */
    public Matrix kroneckerSum(Matrix B) {
        int rows1 = this.getRowDimension();
        int rows2 = B.getRowDimension();
        return this.kronecker(identity(rows2)).plus(
                kronecker(identity(rows1), B));
    }

    /**
     * Computes the kronecker sum of A and B
     * @param A matrix to sum 
     * @param B matrix to sum 
     * @return Kronecker sum of A and B
     */
    public static Matrix kroneckerSum(Matrix A, Matrix B) {
        return (new MarkovMatrix(A)).kroneckerSum(B);
    }

    /**
     * Computes the kronecker product of A and B
     * @param A matrix to multiply
     * @param B matrix to multiply 
     * @return Kronecker product of A and B
     */
    public static Matrix kronecker(Matrix A, Matrix B) {
        int r1 = A.getRowDimension();
        int c1 = A.getColumnDimension();
        int r2 = B.getRowDimension();
        int c2 = B.getColumnDimension();
        Matrix result = new Matrix(r1 * r2, c1 * c2);
        for (int i = 0; i < r1; i++) {
            for (int j = 0; j < c1; j++) {
                result.setMatrix(i * r2, i * r2 + r2 - 1, j * c2, j * c2 + c2
                        - 1, B.times(A.get(i, j)));
            }
        }
        return result;
    }

    /**
     * Computes the kronecker sum of this and B
     * @param B matrix to sum 
     * @return Kronecker sum of this and B
     */
    public Matrix kronecker(Matrix B) {
        return kronecker(this, B);
    }

    /**
     * Concatenates matrices A and B by columns
     * @param A matrix to concatenate
     * @param B matrix to concatenate
     * @return Matrix with the columns of A and B
     */
    public static Matrix concatCols(Matrix A, Matrix B) {
        int r1 = A.getRowDimension();
        int c1 = A.getColumnDimension();
        int r2 = B.getRowDimension();
        int c2 = B.getColumnDimension();
        if (r1 != r2)
            throw new IndexOutOfBoundsException(
                    "Matrices cannot be concatenated: (" + r1 + "," + c1
                            + ")|(" + r2 + "," + c2 + ")");
        Matrix C = new Matrix(r1, c1 + c2);
        C.setMatrix(0, r1 - 1, 0, c1 - 1, A);
        C.setMatrix(0, r1 - 1, c1, c1 + c2 - 1, B);
        return C;
    }

    /**
     * Concatenates matrices A and B by rows
     * @param A matrix to concatenate
     * @param B matrix to concatenate
     * @return Matrix with the rows of A and B
     */
    public static Matrix concatRows(Matrix A, Matrix B) {
        int r1 = A.getRowDimension();
        int c1 = A.getColumnDimension();
        int r2 = B.getRowDimension();
        int c2 = B.getColumnDimension();
        if (c1 != c2)
            throw new IndexOutOfBoundsException(
                    "Matrices cannot be concatenated: (" + r1 + "," + c1
                            + ")--(" + r2 + "," + c2 + ")");
        Matrix C = new Matrix(r1 + r2, c1);
        C.setMatrix(0, r1 - 1, 0, c1 - 1, A);
        C.setMatrix(r1, r1 + r2 - 1, 0, c1 - 1, B);
        return C;
    }

    @Override
    public String toString() {
        ByteArrayOutputStream os = new ByteArrayOutputStream();
        PrintWriter pwt = new PrintWriter(os);
        print(pwt, 5, 2);
        try {
            os.flush();
            pwt.flush();
        } catch (IOException e) {
        }
        ;
        return os.toString();
    }

    public String toTxt() {
        String stg = "";
        int m = this.getRowDimension();
        int n = this.getColumnDimension();
        double A[][] = this.getArray();
        stg += m + ":" + n + ":";
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                stg += A[i][j] + " ";
            }
            stg += ";";
        }
        return stg;
    }

    public String toStringRTF() {
        int r = this.getRowDimension();
        int c = this.getColumnDimension();
        String stg = "";
        for (int i = 0; i < r; i++) {
            stg += new MarkovMatrix(getMatrix(i, i, 0, c - 1)).toString()
                    + "\\par";
        }
        return stg;
    }

    public MarkovMatrix oldExp(double x, Matrix rightMat) {
        double ldaX;
        int maxK = 30000;
        double lOmE = Math.log(1.0 - Epsilon);
        // boolean finished = false;
        MarkovMatrix matP = getNormalized();
        ldaX = ldaMax * x;
        double pk = Math.exp(-ldaX);
        if (pk < Double.MIN_VALUE) {
            pk = Float.MIN_VALUE;
        }
        double logEndLevel = Math.log(pk) + lOmE + ldaX;
        double sumPk = pk;
        double logSumPk;
        // Ck = (P^k)*RightMat
        Matrix matCk = identity(getColumnDimension()); // Ck = (P^k)
        Matrix result = matCk.copy().times(pk);
        int k = 1;
        double alarmLevel = (1e-15) * Float.MAX_VALUE;
        do {
            pk = pk * (ldaX / k);
            sumPk += pk;
            logSumPk = Math.log(sumPk);
            if (sumPk > alarmLevel) {
                pk = pk / sumPk;
                result = result.times(1 / sumPk);
                logEndLevel = logEndLevel - logSumPk;
                sumPk = 1.0;
            }
            matCk = matP.pow(k);
            result = result.plus(matCk.times(pk));
            k++;
        } while ((logSumPk < logEndLevel) && (k < maxK));
        if (k == maxK) {
            throw new NumberFormatException("Max iterations reached");
        }
        return new MarkovMatrix(result.times(1.0 / sumPk).times(rightMat));
    }

    private Matrix runge4(Matrix y, double step) {
        Matrix t1, t2, t3, /* temporary storage arrays */
        k1, k2, k3, k4; /* for Runge-Kutta */

        t1 = y.plus((k1 = y.times(this).times(step)).times(0.5));
        t2 = y.plus((k2 = t1.times(this).times(step)).times(0.5));
        t3 = y.plus(k3 = t2.times(this).times(step));
        k4 = t3.times(this).times(step);

        k2 = k2.times(2.0);
        k1 = k1.times(2.0);

        y = y.plus((k1.plus(k2).plus(k3).plus(k4)).times(1.0 / 6.0));
        return y;
    }

    /**
     * 
     * 
     */
    MarkovMatrix aExpStep = null;

    @SuppressWarnings("unusedPrivate")
    private Matrix runge4B(@SuppressWarnings("unusedArgument")
    double x, Matrix y, double step) {
        if (aExpStep == null) {
            int n = size();
            MarkovMatrix[] y2 = expUnif(2, step, identity(n), identity(n), 4);
            aExpStep = y2[1];
        }
        return y.times(aExpStep);
    }

}
