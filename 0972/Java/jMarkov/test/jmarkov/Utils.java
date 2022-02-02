/*
 * Created on 06-jul-2005
 */
package jmarkov;

import junit.framework.Assert;
import Jama.Matrix;

/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 * 
 */
public class Utils extends Assert {

    /**
     * tests whether the two given array are equal.
     * 
     * @param a1
     *            First array
     * @param a2
     *            Second array
     */
    public static void assertArrayEquals(double[] a1, double[] a2) {
        assertArrayEquals(a1, a2, 1e-8);
    }

    /**
     * Tests whether these two array are equal.
     * 
     * @param a1
     * @param a2
     * @param epsilon
     */
    public static void assertArrayEquals(double[] a1, double[] a2,
            double epsilon) {
        assertArrayEquals("", a1, a2, epsilon);
    }

    /**
     * tests whether the two given array have all values whitin epsilon.
     * 
     * @param msg
     *            Message to write when test fails.
     * 
     * @param a1
     *            First (expected) array.
     * @param a2
     *            Second (actual) array.
     * @param epsilon
     *            Tolerance.
     */
    public static void assertArrayEquals(String msg, double[] a1, double[] a2,
            double epsilon) {
        String formatted = "";
        boolean result = true;
        if (a1 == null)
            result = false;
        else if (a2 == null)
            result = false;
        else if (a1.length != a2.length) {
            result = false;
            formatted = " Diferent array dimensions. Expected: " + a1.length
                    + ", but was :" + a2.length;
        } else {
            int n = a1.length;
            for (int i = 0; i < n && result; i++) {
                result = result && (Math.abs(a1[i] - a2[i]) < epsilon);
                formatted = " Diferent arrays at component " + i
                        + ". Expected: " + a1[i] + ", but was :" + a2[i];
            }
        }
        if (!result) {
            fail(msg + ":" + formatted);
        }
    }

    /**
     * tests whether the two given array have all values whitin epsilon.
     * 
     * @param msg
     *            Message to write when test fails.
     * 
     * @param a1
     *            First (expected) array.
     * @param a2
     *            Second (actual) array.
     */
    public static void assertArrayEquals(String msg, int[] a1, int[] a2) {
        String formatted = "";
        boolean result = true;
        if (a1 == null)
            result = false;
        else if (a2 == null)
            result = false;
        else if (a1.length != a2.length) {
            result = false;
            formatted = " Diferent array dimensions. Expected: " + a1.length
                    + ", but was :" + a2.length;
        } else {
            int n = a1.length;
            for (int i = 0; i < n && result; i++) {
                result = result && (a1[i] == a2[i]);
                formatted = " Diferent arrays at component " + i
                        + ". Expected: " + a1[i] + ", but was :" + a2[i];
            }
        }
        if (!result) {
            fail(msg + ":" + formatted);
        }
    }

    /**
     * Tests whether the two given array are equal.
     * 
     * @param a1
     *            First array[][]
     * @param a2
     *            Second array[][]
     */
    public static void assertArrayEquals(double[][] a1, double[][] a2) {
        assertArrayEquals(a1, a2, 1e-8);
    }

    /**
     * @param a1
     * @param a2
     * @param epsilon
     */
    public static void assertArrayEquals(double[][] a1, double[][] a2,
            double epsilon) {
        assertArrayEquals("", a1, a2, epsilon);
    }

    /**
     * tests whether the two given array have all values whitin epsilon.
     * 
     * @param msg
     *            The message to write when the test fails.
     * 
     * @param a1
     *            First array[][]
     * @param a2
     *            Second array[][]
     * @param epsilon
     *            Tolerance.
     */
    public static void assertArrayEquals(String msg, double[][] a1,
            double[][] a2, double epsilon) {
        boolean result = true;
        String formatted = "";
        if (a1 == null)
            result = false;
        else if (a2 == null)
            result = false;
        else if (a1.length != a2.length)
            result = false;
        else if (a1[1].length != a2[1].length) {
            result = false;
            formatted = " Diferent array dimensions. Expected: " + a1.length
                    + ", but was :" + a2.length;
        } else {
            int n = a1.length;
            int m = a1[1].length;
            for (int i = 0; i < n && result; i++) {
                for (int j = 0; j < m && result; j++) {
                    result = result
                            && (Math.abs(a1[i][j] - a2[i][j]) < epsilon);
                    formatted = " Diferent arrays at component (" + i + ","
                            + "). Expected: " + a1[i][j] + ", but was :"
                            + a2[i][j];

                }
            }
        }
        if (!result)
            fail(msg + ":" + formatted);
    }

    /**
     * Creates a MxM flow line probability.
     * 
     * @param M
     *            size
     * @return A matrx where p(i,i+1) = 1.0
     */
    public static double[][] flowProb(int M) {
        double[][] pr = new double[M][M];
        for (int i = 0; i < M - 1; i++) {
            pr[i][i + 1] = 1.0;
        }
        pr[M - 1][0] = 1.0;
        return pr;
    }

    /**
     * @param N
     *            Number of moving entities
     * @param M
     *            Number of stations
     * @param mu
     *            Process rates.
     * @return troughput rates for the stations.
     */
    public static double[] getMVAthruput(int N, int M, double[] mu) {
        return getMVAthruput(N, M, mu, flowProb(M));
    }

    private static int top_n = 0;
    private static double L[][] = new double[100][];

    /**
     * @param N
     *            Total wip
     * @param M
     *            Number of stations
     * @param mu
     *            Service rate
     * @param alpha
     *            relative importance
     * @return Troughput rates for the stations.
     */
    public static double[] getMVAthruput(int N, int M, double[] mu,
            double[] alpha) {
        double lam[] = new double[M];
        if (N < top_n || top_n == 0) {
            top_n = 0;
            L = new double[100][];
            L[0] = new double[M];
        }
        for (int n = top_n + 1; n <= N; n++) {
            L[n] = new double[M];
            lam = getMVAthruput(n, M, mu, L[n - 1], L[n], alpha);
        }
        top_n = N;
        return lam;
    }

    /**
     * Gets the thruput rate for a system with the given parameters
     * 
     * @param N
     *            Total wip
     * @param M
     *            Number of stations
     * @param mu
     *            Service rate
     * @param prob
     *            Routing probabilities.
     * @return The effective rates as computed using MVA.
     */
    public static double[] getMVAthruput(int N, int M, double[] mu,
            double[][] prob) {
        return getMVAthruput(N, M, mu, getSS(prob));

    }

    /**
     * Resets MV parameters.
     */
    public static void resetMVA() {
        top_n = 0;
        L = new double[100][];
    }

    private static double[] getMVAthruput(int n, int M, double[] mu,
            double[] L_old, double[] L, double alpha[]) {
        double W[] = new double[M];
        double lam[] = new double[M];
        double sum = 0.0;
        for (int m = 0; m < M; m++) {
            W[m] = (1.0 + L_old[m]) * (1 / mu[m]);
            sum += alpha[m] * W[m];
        }
        for (int m = 0; m < M; m++) {
            L[m] = n * W[m] * alpha[m] / sum;
            lam[m] = L[m] / W[m];
        }
        return lam;
    }

    /**
     * Gets a probability matrix
     * 
     * @param M
     *            Stations
     * @return a routing matrix.
     */
    public static double[][] getProb(int M) {
        double prob[][] = new double[M][M];
        for (int i = 0; i < M; i++) {
            double sum = 0.0;
            for (int j = 0; j < M - 1; j++)
                sum += (prob[i][j] = 1.0 / (2 * M));
            prob[i][M - 1] = 1.0 - sum;
        }
        return prob;
    }

    /**
     * Gets a setady state probabilities for the fgieven matrix.
     * 
     * @param prob
     *            Markov probability matrix.
     * @return SS prob vector.
     */
    public static double[] getSS(double[][] prob) {
        int n = prob.length;
        Matrix matA = Matrix.identity(n, n).minus(new Matrix(prob));
        for (int i = 0; i < n; i++) {
            matA.set(i, n - 1, 1.0); // replace last column with ones
        }
        Matrix oneVec = new Matrix(1, prob.length);// unit vector
        oneVec.set(0, n - 1, 1);
        return matA.solveTranspose(oneVec).getColumnPackedCopy();
    }
}
