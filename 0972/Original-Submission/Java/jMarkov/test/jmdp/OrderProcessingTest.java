// 2014 - Tests pass
package jmdp;

import junit.framework.TestCase;
import jmarkov.Utils;
import examples.jmdp.OrderProcessing;

/**
 * @author German Riano, Andres Sarmiento. Universidad de los Andes. (C) 2006
 * 
 */
public class OrderProcessingTest extends TestCase {

    OrderProcessing prob = null;

    @Override
    protected void setUp() throws Exception {
        double setupCost = 60;
        double unfilledOrderCost = 4;
        int maxOrders = 15;
        double theta = 3;
        // double discountFactor = 0.9;
        prob = new OrderProcessing(0, setupCost, unfilledOrderCost, maxOrders,
                theta);
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
        prob = null;
    }

    /**
     * @throws Exception
     * 
     */
    public void testGetOptimalValueFunction() throws Exception {
        prob.solve(1 / 0.9 - 1.0);
        double vf[] = prob.getOptimalValueFunction().get();
        // double vf_cero = 10.4724;
        double expValues[] = { 259.42, 270.60, 280.78, 289.90, 297.92, 304.80,
                310.53, 315.29, 319.42, 319.42, 319.42, 319.42, 319.42, 319.42, 319.42,
                319.42 };
        prob.getSolver().setPrintValueFunction(true);
        prob.printSolution();
        Utils.assertArrayEquals("Value Function not equal", expValues, vf,
                1.0E-2);
    }


    /**
     * @param args
     */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(OrderProcessingTest.class);
    }
}
