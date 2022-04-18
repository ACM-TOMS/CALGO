// 2014 Passes all tests that don't need Xpress-MP. Xpress-MP tests commented out
package jmdp;

//XPRESS SOLUTION:
//States ={0,1,2,3,4,5,6,7,8,9}
//Computing parameters...
//Writing Probabilities and Costs to CSV files.
//Solving...
//x has 55 components
//Problem Status    : Optimal found.
//Cost: -2058.0203
//Optimal Policy:
//In state '0' (prob= 0.214870) take action '6'
//In state '1' (prob= 0.156293) take action '5'
//In state '2' (prob= 0.195367) take action '4'
//In state '3' (prob= 0.195367) take action '3'
//In state '4' (prob= 0.146525) take action '2'
//In state '5' (prob= 0.073263) take action '1'
//In state '6' (prob= 0.018316) take action '0'
//
//
//REPORT
//State                    0         1         2         3         4         5         6         7         8         9
//Prob                0.215     0.156     0.195     0.195     0.147     0.073     0.018     0.000     0.000     0.000
//Init Inv.               0         1         2         3         4         5         6         7         8         9
//Buy                     6         5         4         3         2         1         0         0         0         0
//Availabl.               6         6         6         6         6         6         6         7         8         9
//Exp. Sales          3.805     3.805     3.805     3.805     3.805     3.805     3.805     3.915     3.966     3.988
//Lost Sales          0.195     0.195     0.195     0.195     0.195     0.195     0.195     0.085     0.034     0.012
//Sales Income     3804.565  3804.565  3804.565  3804.565  3804.565  3804.565  3804.565  3915.239  3966.373  3987.736
//Truck Cost         50.000    50.000    50.000    50.000    50.000    50.000     0.000     0.000     0.000     0.000
//Order Cost       2400.000  2000.000  1600.000  1200.000   800.000   400.000     0.000     0.000     0.000     0.000
//Holding Cost        0.000    80.000   160.000   240.000   320.000   400.000   480.000   560.000   640.000   720.000
//Backorder Cost      0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000     0.000
//Tot. Oprt. Cost   167.261   247.261   327.261   407.261   487.261   567.261   597.261   610.856   660.176   727.358
//Net Income       1354.565  1674.565  1994.565  2314.565  2634.565  2954.565  3324.565  3355.239  3326.373  3267.736

import jmarkov.jmdp.solvers.LPBCLAverageSolver;
import jmarkov.jmdp.solvers.MPSQsOptAverageSolver;
import jmarkov.jmdp.solvers.MPSXpressAverage;
import jmarkov.jmdp.solvers.RelativeValueIterationSolver;
import jmarkov.jmdp.solvers.Solver;
import junit.framework.TestCase;
import examples.jmdp.ControlProdNonEvents;
import examples.jmdp.ControlProduccion;
import examples.jmdp.InvLevel;
import examples.jmdp.Order;

/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 * 
 */
public class ControlProduccionTest extends TestCase {

    ControlProduccion prob2 = null;
    ControlProdNonEvents prob = null;
    private double intRate = 0.9;

    @Override
    protected void setUp() throws Exception {
        int M = 9;// Capacity
        double K = 50;// Fixed cost per Order
        double cost = 400; // variable cost
        double price = 1000; // variable cost
        double holdingCost = 80; // holding cost per computer per period.
        double demandMean = 4; // mean of the Poisson demand per stage

        prob2 = new ControlProduccion(M, K, cost, price, holdingCost, intRate,
                demandMean);
        prob = new ControlProdNonEvents(M, K, cost, price, holdingCost,
                intRate, demandMean);
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
        prob = null;
        prob2 = null;
    }

    /**
     * @throws Exception
     */
    public void testGetOptimalValueFunctionValueIteration() throws Exception {
        Solver<InvLevel, Order> solv = new RelativeValueIterationSolver<InvLevel, Order>(
                prob);
        prob.setSolver(solv);
        InvLevel initialState = new InvLevel(0);
        double vf_cero = -2058.0203;// -1765.8024
        double vf = prob.getOptimalValueFunction().get(initialState);
        // prob.printMatrices();
        // prob.getSolver().setPrintValueFunction(true);
        // prob.printSolution();
        assertEquals("Value Function not equal", vf_cero, vf, 1.0E-3);
    }

    /**
     * @throws Exception
     */
    public void testGetOptimalNonEventsValueFunctionValueIteration()
            throws Exception {
        Solver<InvLevel, Order> solv = new RelativeValueIterationSolver<InvLevel, Order>(
                prob2);
        prob2.setSolver(solv);
        InvLevel initialState = new InvLevel(0);
        double vf_cero = -2058.0203;// -1765.8024
//        prob2.getSolver().setPrintValueFunction(true);
//        prob2.printSolution();
        double vf = prob.getOptimalValueFunction().get(initialState);
        assertEquals("Value Function not equal", vf_cero, vf, 1.0E-3);
    }

    
/* Xpress - BCL tests have been disables so tests can be ran succesfully without third party software
 * if you want to test the BCL classes, just uncomment this section    
    *//**
     * @throws Exception
     *//*
    public void testGetOptimalNonEventsAverageLPMPSQsOptSolver() throws Exception {
        InvLevel initialState = new InvLevel(0);
        Solver<InvLevel, Order> solv = new MPSQsOptAverageSolver<InvLevel, Order>(
                prob, "examples\\jmdp\\MPSFolder", "ContProd.mps");
        prob.setSolver(solv);
        prob.solve();
        double vf_cero = -2058.0203;// -1765.8024
        double vf = prob.getOptimalValueFunction().get(initialState);
        assertEquals("Value Function not equal.", vf_cero, vf, 1.0E-3);

    }
    *//**
     * @throws Exception
     *//*
    public void testGetOptimalNonEventsAverageLPMPSXpressSolver() throws Exception {
        InvLevel initialState = new InvLevel(0);
        Solver<InvLevel, Order> solv = new MPSXpressAverage<InvLevel, Order>(
                prob, "examples\\jmdp\\MPSFolder", "ContProd.mps");
        prob.setSolver(solv);
        prob.solve();
        double vf_cero = -2058.0203;// -1765.8024
        double vf = prob.getOptimalValueFunction().get(initialState);
        assertEquals("Value Function not equal.", vf_cero, vf, 1.0E-3);

    }

    *//**
     * @throws Exception
     *//*
    public void testGetOptimalNonEventsAverageLPSolver() throws Exception {
        InvLevel initialState = new InvLevel(0);
        Solver<InvLevel, Order> solv = new LPBCLAverageSolver<InvLevel, Order>(
                prob);
        prob.setSolver(solv);
        double vf_cero = -2058.0203;// -1765.8024
        double vf = prob.getOptimalValueFunction().get(initialState);
        assertEquals("Value Function not equal.", vf_cero, vf, 1.0E-3);
    }*/

    /**
     * @param args
     */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(ControlProduccionTest.class);
    }
}
