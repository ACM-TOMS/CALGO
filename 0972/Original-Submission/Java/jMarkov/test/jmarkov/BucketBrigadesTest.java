/**
 * BucketBrigades.java
 * Created: Feb 18, 2006
 */
package jmarkov;

import jmarkov.basic.exceptions.NotUnichainException;
import junit.framework.TestCase;
import examples.jmarkov.BucketBrigades;
import examples.jmarkov.BucketBuffers;

/**
 * @author Germán Riaño. Universidad de los Andes. (C) 2006
 * 
 */
public class BucketBrigadesTest extends TestCase {

    /**
     * @param args
     */
    public static void main(String[] args) {
        //junit.swingui.TestRunner.run(BucketBrigadesTest.class);
        junit.textui.TestRunner.run(BucketBrigadesTest.class);
    }

    /**
     * 
     * test BucketBrigades for identical expo servers
     * @throws NotUnichainException 
     */
    public void testIdenticallExpo() throws NotUnichainException {
        double mu = 1.0;// service rate
        BucketBrigades bb = null;
        for (int M = 1; M <= 5; M += 2) {
            for (int N = 1; N <= 5; N += 2) {// Operators
                double vels[][] = new double[N][M];
                for (int n = 0; n < N; n++) {
                    for (int m = 0; m < M; m++) {
                        vels[n][m] = 1.0;
                    }
                }
                double expTh = N * mu / (M + N - 1);
                bb = new BucketBrigades(N, M, vels);
                bb.setDebugLevel(0);
                bb.setMaxStates(2000);
                double val = bb.getEventRate(M - 1);
                bb.debug(1, "M=" + M + ", N= " + N + ", States = "
                        + bb.getNumStates());
                assertEquals("Error in Thruput ", expTh, val, 1e-10);
            }
        }
    }

    /**
     * Compares the BB system with a CONWIP system results. This requires
     * identical exponential workers (but stations may differ).
     * @throws NotUnichainException 
     */
    public void testIdenticalWorkers() throws NotUnichainException {
        BucketBrigades bb = null;
        Utils.resetMVA();
        for (int M = 1; M <= 7; M += 2) {
            for (int N = 1; N <= 5; N += 2) {// Operators
                double mu[] = new double[M];// basic service rate
                double vels[][] = new double[N][M];
                for (int n = 0; n < N; n++) {
                    for (int m = 0; m < M; m++) {
                        mu[m] = m + 1;
                        vels[n][m] = 1.0 * mu[m];
                    }
                }
                // Jackson jacky = new Jackson(N, M, mu);
                // double expTh = jacky.getEventRate(M - 1);
                double expTh = Utils.getMVAthruput(N, M, mu)[M - 1];
                bb = new BucketBrigades(N, M, vels);
                bb.setDebugLevel(0);
                bb.setMaxStates(2000);
                double val = bb.getEventRate(M - 1);
                bb.debug(1, "M=" + M + ", N= " + N + ", States = "
                        + bb.getNumStates());
                assertEquals("Error in Thruput ", expTh, val, 1e-5);
            }
        }
    }

    /**
     * 
     * test BucketBuffers for identical expo servers
     * @throws NotUnichainException 
     */
    public void testIdenticallExpoBuff() throws NotUnichainException {
        double mu = 1.0;// service rate
        BucketBuffers bbb = null;
        for (int M = 1; M <= 7; M += 2) {
            for (int N = 1; N <= 5; N += 2) {// Operators
                double vels[][] = new double[N][M];
                for (int n = 0; n < N; n++) {
                    for (int m = 0; m < M; m++) {
                        vels[n][m] = 1.0;
                    }
                }
                double expTh = N * mu / (M + N - 1);
                bbb = new BucketBuffers(vels, new int[M]);
                bbb.setDebugLevel(0);
                bbb.setMaxStates(2000);
                double val = bbb.getEventRate(M - 1);
                bbb.debug(1, "M=" + M + ", N= " + N + ", States = "
                        + bbb.getNumStates());
                assertEquals("Error in Thruput ", expTh, val, 1e-10);
            }
        }
    }

    /**
     * Compares the BucketBuffers system with a CONWIP system results. This
     * requires identical exponential workers (but stations may differ).
     * @throws NotUnichainException 
     */
    public void testIdenticalWorkersBuff() throws NotUnichainException {
        BucketBuffers bbb = null;
        Utils.resetMVA();
        for (int M = 1; M <= 7; M += 2) {
            for (int N = 1; N <= 5; N += 2) {// Operators
                double mu[] = new double[M];// basic service rate
                double vels[][] = new double[N][M];
                for (int n = 0; n < N; n++) {
                    for (int m = 0; m < M; m++) {
                        mu[m] = m + 1;
                        vels[n][m] = 1.0 * mu[m];
                    }
                }
                // Jackson jacky = new Jackson(N, M, mu);
                // double expTh = jacky.getEventRate(M - 1);
                double expTh = Utils.getMVAthruput(N, M, mu)[M - 1];
                bbb = new BucketBuffers( vels, new int[M]);
                bbb.setDebugLevel(0);
                bbb.setMaxStates(2000);
                double val = bbb.getEventRate(M - 1);
                bbb.debug(1, "M=" + M + ", N= " + N + ", States = "
                        + bbb.getNumStates());
                assertEquals("Error in Thruput ", expTh, val, 1e-5);
            }
        }
    }
}