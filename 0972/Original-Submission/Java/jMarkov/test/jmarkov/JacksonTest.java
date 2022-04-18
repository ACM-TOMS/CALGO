/*
 * Created on 05-jul-2005
 */
package jmarkov;

import static jmarkov.Utils.assertArrayEquals;
import static jmarkov.Utils.getMVAthruput;
import static jmarkov.Utils.getProb;
import jmarkov.MarkovProcess;
import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.exceptions.NotUnichainException;
import junit.framework.TestCase;
import examples.jmarkov.Jackson;

/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 * 
 */
public class JacksonTest extends TestCase {
    private Jackson jackie = null;

    /**
     * @param args
     */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(JacksonTest.class);
    }

    /**
     * @param name
     */
    public JacksonTest(String name) {
        super(name);
        int wip[] = { 10, 10, 15 };
        int servers[] = { 2, 2, 3 };
        double rates[] = { 10.0, 8.0, 12.0 };
        double probs[][] = { { 0.0, 0.5, 0.5 }, { 0, 0, 1 }, { 1, 0, 0 } };
        jackie = new Jackson(wip, servers, rates, probs);
        jackie.setDebugLevel(0);
        // jackie.setSteadyStateSolver(new JamaSolver(jackie));
    }

    /**
     * 
     */
    public void testGenerate() {
        jackie.generate();
        assertTrue(jackie.getStatus() == MarkovProcess.Status.GENERATED);
    }

    /**
     * 
     */
    public void oldtestStatusBig() {
        Jackson jacko = new Jackson("JacksonFiles/BigWiplevel.txt",
                "JacksonFiles/BigNumservers.txt",
                "JacksonFiles/BigServicesrates.txt",
                "JacksonFiles/BigProbabilities.txt");
        jacko.setMaxStates(Long.MAX_VALUE);
        jacko.generate();
        // double mops[] = jacko.getMOPsAvg();
        assertTrue(jacko.getStatus() == SimpleMarkovProcess.Status.GENERATED);
    }

    /**
     * @throws NotUnichainException 
     * 
     */
    public void testGetSteadyState() throws NotUnichainException {
        double prob[] = jackie.getSteadyState();
        assertEquals("Wrong prob value", 0.03985, prob[665], 1.0E-3);
    }

    /**
     * @throws NotUnichainException 
     * 
     */
    public void testMOPsAvg() throws NotUnichainException {
        double mops[] = jackie.getMOPsAvg();
        assertArrayEquals("MOPS not equal: ", new double[] {
                1.9999978951064057, 1.2500006423494994, 1.6666593265864875,
                30.90737210525042, 2.0512648064879384, 2.041335065754373 },
                mops, 1e-4);
    }

    /**
     * Compares with MVA results.
     * @throws NotUnichainException 
     * 
     */
    public void testMVA() throws NotUnichainException {
        Utils.resetMVA();
        for (int M = 1; M <= 5; M = M + 2) {// Stations
            for (int N = 1; N <= 10; N = N + 2) {// WIP
                double mu[] = new double[M];
                for (int m = 0; m < M; m++)
                    mu[m] = 1 / (m + 1.0);
                double[][] prob = getProb(M);
                double expThr[] = getMVAthruput(N, M, mu, prob);
                Jackson jacky = new Jackson(N, M, mu, prob);
                jacky.setDebugLevel(0);
                double thr[] = jacky.effLambdas();
                assertArrayEquals(
                        "Thruput not equal comparing with MVA, for N = " + N
                                + " and M = " + M + ". ", expThr, thr, 1e-5);
                jacky = null;
            }
        }

    }

}
