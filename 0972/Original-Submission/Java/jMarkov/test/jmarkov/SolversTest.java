/*
 * Created on 06-jul-2005
 */
package jmarkov;

import static jmarkov.Utils.assertArrayEquals;

import java.io.PrintWriter;

import jmarkov.MarkovProcess;
import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.exceptions.NotUnichainException;
import jmarkov.solvers.JamaSolver;
import jmarkov.solvers.MtjSolver;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import Jama.Matrix;
import examples.jmarkov.DriveThru;
import examples.jmarkov.Jackson;
import examples.jmarkov.Kanban;
import examples.jmarkov.QueueMM1N;
import examples.jmarkov.QueueMM2dN;

/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 * 
 */
public class SolversTest extends TestCase {

    /**
     * @param args
     */
    public static void main(String[] args) {
        //junit.swingui.TestRunner.run(AllTests.class);
    	try{
        junit.textui.TestRunner.runAndWait(suite());
    	}catch (Exception e){
    		System.out.println("Exception found");
    		e.printStackTrace();
    	}
    }

    /**
     * @return A suite witl all tests for all solvers.
     * @throws Exception 
     */
    public static Test suite() throws Exception {
        // $JUnit-BEGIN$
        TestSuite suite = new TestSuite("Test for jmarkov.tests");
        final SimpleMarkovProcess<?, ?>[] mps = { new DriveThru(),
                new Jackson(), new Kanban(), new QueueMM1N(), new QueueMM2dN(),
        // new BBPhBuf()
        };

        for (final MarkovProcess<?, ?> mp : mps) {
            mp.setDebugLevel(0);
            mp.generate();
            mp.setDebugLevel(1);
            mp.setSteadyStateSolver(new JamaSolver(mp));
            final double pi0[] = mp.getSteadyState();

            for (final MtjSolver.EnumSolver is : MtjSolver.EnumSolver.values()) {
                String testName = mp.getClass().getSimpleName() + " -- "
                        + is.toString();
                suite.addTest(new TestCase(testName) {
                    private MarkovProcess<?, ?> locMp = null;

                    @Override
                    protected void setUp() throws Exception {
                        locMp = mp;
                        mp.setDebugLevel(0);
                    }

                    @Override
                    protected void runTest() throws Throwable {

                        MtjSolver solver = new MtjSolver(locMp, is, false);
                        try {
                            locMp.setSteadyStateSolver(solver);
                            double pi[] = solver.getSteadyState();
                            assertArrayEquals("Solver test with " + is
                                    + ". Result is not equal to Jama's", pi,
                                    pi0, 1e-2);
                            System.out.printf("Solver %s with " + locMp.label()
                                    + ". Process Time: %d milliseconds\n", is,
                                    (solver.getProcessTime()));
                        } catch (NotUnichainException e) {
                            System.out
                                    .printf(
                                            "Solver %s with "
                                                    + locMp.label()
                                                    + " did not converged!. Process Time: %d milliseconds\n",
                                            is, (solver.getProcessTime()));
                        }
                    }
                });
            }
        }
        // $JUnit-END$
        return suite;
    }

    static void test(SimpleMarkovProcess mp) throws NotUnichainException {
        test(mp, false);
    }

    static void test(SimpleMarkovProcess<?, ?> mp, boolean showRes) throws NotUnichainException {
        mp.setMaxStates(Integer.MAX_VALUE);
        mp.setDebugLevel(0);
        mp.generate();
        // mp.setDebugLevel(2);
        PrintWriter pw = new PrintWriter(System.out, true);
        mp.setSteadyStateSolver(new JamaSolver(mp));
        if (showRes) {
            pw.println("JAMA REWSULTS: ");
            mp.printStates(pw);
        }
        double pi0[] = mp.getSteadyState();
        MtjSolver mtjs = new MtjSolver(mp);
        for (MtjSolver.EnumSolver s : MtjSolver.EnumSolver.values()) {
            mtjs.setCurrentIterSolver(s);
            mp.setSteadyStateSolver(mtjs);
            double pi[] = mp.getSteadyState();
            assertArrayEquals(
                    "SOLVER: " + s + " result is NOT equal to JAMA's", pi, pi0,
                    0.001);
            // assertTrue("SOLVER: " + s + " result is NOT equal to
            // JAMA's",equal);
            // pw.println("SOLVER: " + s + ". Result is "
            // + ((equal) ? "" : "NOT ") + "equal to JAMA's");
            if (showRes) {
                Matrix piM = new Jama.Matrix(new double[][] { pi });
                pw.print("Result = ");
                piM.print(pw, 10, 10);
                // mp.printStates(pw);
            }
        }
    }
}
