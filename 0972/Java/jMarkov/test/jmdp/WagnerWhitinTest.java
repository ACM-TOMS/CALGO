// 2014: Passes all test
/**
 * WagnerWhitinTest.java
 * Created: Jul 15, 2005
 */
package jmdp;

import jmarkov.jmdp.solvers.FiniteSolver;
import junit.framework.TestCase;
import examples.jmdp.InvLevel;
import examples.jmdp.Order;
import examples.jmdp.WagnerWhitin;

/**
 * @author Germán Riaño. Universidad de los Andes. (C) 2005
 */
public class WagnerWhitinTest extends TestCase {

    WagnerWhitin prob = null;

    /**
     * @param args
     */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(WagnerWhitinTest.class);
    }

    /**
     * Constructor for WagnerWhitinTest.
     * 
     * @param name
     */
    public WagnerWhitinTest(String name) {
        super(name);
    }

    /*
     * @see TestCase#setUp()
     */
    @Override
    protected void setUp() throws Exception {
        int lastStage = 12;
        int maxInventory = 15;
        int maxBackorders = 5;
        int truckSize = 6;
        double K = 500;
        double b = 2000;
        double p = 22000;
        double c = 20000;
        double h = Math.pow(1.3, 1.0 / 52) - 1.0;
        int[] demand = new int[] { 10, 4, 3, 6, 3, 2, 0, 1, 7, 3, 4, 5 };
        prob = new WagnerWhitin(0, lastStage, maxInventory,
                maxBackorders, truckSize, K, b, p, c, h, demand);
    }

    /*
     * @see TestCase#tearDown()
     */
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
        prob = null;
    }

    /*
     * Test method for 'jmdp.MDP.solve()'
     */
    /**
     * 
     */
    public void testSolve() {
        FiniteSolver<InvLevel, Order> theSolver = new FiniteSolver<InvLevel, Order>(
                prob);
        theSolver.solve();
    }

    /**
     * @throws Exception
     */
    // public void testMain(){
    // WagnerWhitin.main(new String[]{});
    // }

    // public void testgetOptimalSolution(){
    // WagnerWhitin.main(new String[]{});
    // }

    public void testgetOptimalValueFunction() throws Exception {
        InvLevel initial = new InvLevel(0);
        double vf = prob.getOptimalValueFunction().get(initial);
        double vf_cero = 1005.8216;
        assertEquals("Value Function not equal", vf_cero, vf, 1.0E-3);
    }
}
