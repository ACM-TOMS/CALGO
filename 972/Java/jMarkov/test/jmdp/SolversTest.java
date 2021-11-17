/*
 * Created on 24/12/2005
 */
package jmdp;

import jmarkov.jmdp.solvers.LPBCLAverageSolver;
import jmarkov.jmdp.solvers.LPBCLDiscountedSolver;
import jmarkov.jmdp.solvers.MPSQsOptAverageSolver;
import jmarkov.jmdp.solvers.MPSQsOptDiscountedSolver;
import jmarkov.jmdp.solvers.MPSXpressAverage;
import jmarkov.jmdp.solvers.MPSXpressDiscounted;
import jmarkov.jmdp.solvers.PolicyIterationSolver;
import jmarkov.jmdp.solvers.PolicyIterationSolverAvg;
import jmarkov.jmdp.solvers.RelativeValueIterationSolver;
import jmarkov.jmdp.solvers.Solver;
import jmarkov.jmdp.solvers.ValueIterationSolver;
import junit.framework.TestCase;
import jmarkov.Utils;
import examples.jmdp.ControlProdNonEvents;

// XPRESS (AVG) solution: Cost: -1390.6495
// Optimal Policy:
// In state 'Level 0' (prob= 0.580097) take action 'Order 4'
// In state 'Level 1' (prob= 0.191975) take action 'Order 3'
// In state 'Level 2' (prob= 0.141437) take action 'Order 2'
// In state 'Level 3' (prob= 0.069447) take action 'Order 0'
// In state 'Level 4' (prob= 0.017044) take action 'Order 0'
//
// With disc = 0.9 XPRESS gives: Cost: -14422.4654
// Optimal Policy:
// In state 'Level 0' (val=-13644.2762) take action 'Order 4'
// In state 'Level 1' (val=-13964.2762) take action 'Order 3'
// In state 'Level 2' (val=-14284.2762) take action 'Order 2'
// In state 'Level 3' (val=-14795.2222) take action 'Order 0'
// In state 'Level 4' (val=-15424.2762) take action 'Order 0'

/**
 * @author Germán Riaño. Universidad de los Andes.
 * 
 */
@SuppressWarnings("unchecked")
public class SolversTest extends TestCase {

    String workingDir = "src\\examples\\jmdp\\MPSFolder";
    ControlProdNonEvents prob = null;
    Solver solver;
    double disFactor = 0.9;
    double intRate = (1.0 / disFactor) - 1.0;
    double discOpt[] = { -13644.2762, -13964.2762, -14284.2762, -14795.2222,
            -15424.2762 };
    double avgOpt[] = { -1390.6495, -1390.6495, -1390.6495, -1390.6495,
            -1390.6495 };
    int actOpt[] = { 4, 3, 2, 0, 0 };

    @SuppressWarnings("unchecked")
    /*
     * (non-Javadoc)
     * 
     * @see junit.framework.TestCase#setUp()
     */
    @Override
    protected void setUp() throws Exception {
        super.setUp();
        int M = 4;// Capacity
        double K = 500;// Fixed cost per ShortOrder
        double cost = 400; // variable cost
        double price = 1000; // variable cost
        double holdingCost = 80; // holding cost per computer per period.
        double demandMean = 4.0; // mean of the Poisson demand per stage
        prob = new ControlProdNonEvents(M, K, cost, price, holdingCost,
                intRate, demandMean);

    }

    /*
     * @see junit.framework.TestCase#tearDown()
     */
     @Override
     protected void tearDown() throws Exception {
         super.tearDown();
         prob = null;
         System.gc();
     }

    /**
     * @throws Exception
     */
    /*
     * Discounted tests!!
     */

    public void testValueIteration() throws Exception {
        solver = new ValueIterationSolver(prob, intRate);
        prob.getSolver().setPrintValueFunction(true);
        // prob.printSolution();
        double vf[] = prob.getOptimalValueFunction().get();
        Utils
                .assertArrayEquals("Value Function not equal", discOpt, vf,
                        1.0E-3);
        int ord[] = prob.getOptimalOrderSize();
        Utils.assertArrayEquals("Value Function not equal", actOpt, ord);
    }

    /**
     * @throws Exception
     */
    public void testPolicyIteration() throws Exception {
        solver = new PolicyIterationSolver(prob, intRate);
        prob.setSolver(solver);
        //prob.setDebugLevel(3);
        prob.getSolver().setPrintValueFunction(true);
        // prob.printSolution();
        double vf[] = prob.getOptimalValueFunction().get();
        Utils
                .assertArrayEquals("Value Function not equal", discOpt, vf,
                        1.0E-3);
        int ord[] = prob.getOptimalOrderSize();
        Utils.assertArrayEquals("Value Function not equal", actOpt, ord);
    }
    /**
     * @throws Exception
     */
    public void testModifiedPolicyIteration() throws Exception {
        PolicyIterationSolver solver = new PolicyIterationSolver(prob, intRate);
        prob.setSolver(solver);
        //prob.setDebugLevel(3);
        solver.setModifiedPolicy(true);
        prob.getSolver().setPrintValueFunction(true);
        // prob.printSolution();
        double vf[] = prob.getOptimalValueFunction().get();
        Utils
                .assertArrayEquals("Value Function not equal", discOpt, vf,
                        1.0E-3);
        int ord[] = prob.getOptimalOrderSize();
        Utils.assertArrayEquals("Value Function not equal", actOpt, ord);
    }

    
    //
    // Average Tests
    //
    /**
     * @throws Exception
     */
    public void testRelativeValueIterSolver() throws Exception {
        solver = new RelativeValueIterationSolver(prob);
        prob.setSolver(solver);
        prob.getSolver().setPrintValueFunction(true);
        // prob.printSolution();
        double vf[] = prob.getOptimalValueFunction().get();
        Utils.assertArrayEquals("Value Function not equal", avgOpt, vf, 1.0E-3);
        int ord[] = prob.getOptimalOrderSize();
        Utils.assertArrayEquals("Value Function not equal", actOpt, ord);
    }
    /**
     * @throws Exception
     */
    public void testPolicyIterationSolverAvg() throws Exception {
    	PolicyIterationSolverAvg solver = new PolicyIterationSolverAvg(prob);
        prob.setSolver(solver);
        prob.getSolver().setPrintValueFunction(true);
        // prob.printSolution();
        double vf[] = prob.getOptimalValueFunction().get();
        int ord[] = prob.getOptimalOrderSize();
        Utils.assertEquals("Value Function not equal", avgOpt[0], vf[0], 1.0E-2);
        Utils.assertArrayEquals("Value Function not equal", actOpt, ord);
    }
    //
    // MPS - QSOpt tests
    //
    /**
     * @throws Exception
     */
    public void testMPSQsDiscounted() throws Exception {
        solver = new MPSQsOptDiscountedSolver(prob, intRate, workingDir,
                "MDPQsDisc.mps");
        prob.setSolver(solver);
        prob.getSolver().setPrintValueFunction(true);
        // prob.printSolution();
        double vf[] = prob.getOptimalValueFunction().get();
        Utils
                .assertArrayEquals("Value Function not equal", discOpt, vf,
                        1.0E-2);
        int ord[] = prob.getOptimalOrderSize();
        Utils.assertArrayEquals("Value Function not equal", actOpt, ord);
    }
  

    /**
     * @throws Exception
     */
    public void testMPSQsAverage() throws Exception {
        solver = new MPSQsOptAverageSolver(prob, workingDir, "MDPQsAvg.mps");
        prob.setSolver(solver);
        prob.getSolver().setPrintValueFunction(true);
        // prob.printSolution();
        double vf[] = prob.getOptimalValueFunction().get();
        Utils.assertArrayEquals("Value Function not equal", avgOpt, vf, 1.0E-2);
        int ord[] = prob.getOptimalOrderSize();
        Utils.assertArrayEquals("Value Function not equal", actOpt, ord);
    }

/*  THE TESTS BELOW REQUIRE A PROFESSIONAL VERSION OF XPRESS-MP TO RUN 
 * 	IF YOU WISH TO RUN THESE TESTS INSTALL XPRESS-MP PROFESSIONAL WITH 
 *  THE CORRESPONDING BCL LIBRARIES FOR JAVA, THEN UNCOMMENT BELOW AND RUN AGAIN.
 * //**
     * @throws Exception
     *//*
    public void testLPDiscounted() throws Exception {
        solver = new LPBCLDiscountedSolver(prob, intRate);
        prob.setSolver(solver);
        prob.getSolver().setPrintValueFunction(true);
        // prob.printSolution();
        double vf[] = prob.getOptimalValueFunction().get();
        Utils
                .assertArrayEquals("Value Function not equal", discOpt, vf,
                        1.0E-3);
        int ord[] = prob.getOptimalOrderSize();
        Utils.assertArrayEquals("Value Function not equal", actOpt, ord);
    }

    *//**
     * @throws Exception
     *//*
    
     * AVERAGE SOLVERS
     
    public void testLPAvgSolver() throws Exception {
        solver = new LPBCLAverageSolver(prob);
        prob.setSolver(solver);
        prob.getSolver().setPrintValueFunction(true);
        // prob.printSolution();
        double vf[] = prob.getOptimalValueFunction().get();
        Utils.assertArrayEquals("Value Function not equal", avgOpt, vf, 1.0E-3);
        int ord[] = prob.getOptimalOrderSize();
        Utils.assertArrayEquals("Value Function not equal", actOpt, ord);
    }

    *//**
     * @throws Exception
     *//*
    public void testMPSXpressDisc() throws Exception {
        MPSXpressDiscounted solver = new MPSXpressDiscounted(prob, intRate,
                workingDir, "MDPdisc.mps");
        prob.setSolver(solver);
        // solver.setShowXpressOutput(true);
        prob.getSolver().setPrintValueFunction(true);
        // prob.printSolution();
        double vf[] = prob.getOptimalValueFunction().get();
        Utils
                .assertArrayEquals("Value Function not equal", discOpt, vf,
                        1.0E-2);
    }

    *//**
     * @throws Exception
     *//*
    public void testMPSAvgXpress() throws Exception {
        solver = new MPSXpressAverage(prob, workingDir, "MDPavg.mps");
        prob.setSolver(solver);
        prob.getSolver().setPrintValueFunction(true);
        // prob.printSolution();
        double vf[] = prob.getOptimalValueFunction().get();
        Utils.assertArrayEquals("Value Function not equal", avgOpt, vf, 1.0E-2);
        int ord[] = prob.getOptimalOrderSize();
        Utils.assertArrayEquals("Value Function not equal", actOpt, ord);
    }
    *//**
     * @throws Exception
     *//*
    public void testMPSXpressDiscountedNoPath() throws Exception {
        MPSXpressDiscounted solver = new MPSXpressDiscounted(prob, intRate);
        solver.setShowXpressOutput(true);
        prob.setSolver(solver);
        prob.getSolver().setPrintValueFunction(true);
        // prob.printSolution();
        double vf[] = prob.getOptimalValueFunction().get();
        Utils
                .assertArrayEquals("Value Function not equal", discOpt, vf,
                        1.0E-2);
        int ord[] = prob.getOptimalOrderSize();
        Utils.assertArrayEquals("Value Function not equal", actOpt, ord);
    }*/
    /**
     * @param args
     */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(SolversTest.class);
    }

}
