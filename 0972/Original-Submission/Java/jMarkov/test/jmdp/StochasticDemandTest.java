// 2014: Passes all test
package jmdp;

import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import junit.framework.TestCase;
import examples.jmdp.InvLevel;
import examples.jmdp.StochasticDemand;

/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 *
 */
public class StochasticDemandTest extends TestCase {

	StochasticDemand prob = null;
	
	@Override
	protected void setUp() throws Exception {
		int lastStage = 12;
		int maxInventory = 15;
		int maxBackorders = 5;
		int truckSize = 6;
		int K = 500;
		double b = 1000;
		double h = Math.pow(1.3, 1 / 52);
		double theta = 4;
		double price = 22000;
		double cost = 20000;
		InvLevel initial = new InvLevel(0);
		States<InvLevel> initSet = new StatesSet<InvLevel>(initial);
		
		prob = new StochasticDemand(initSet, lastStage, maxInventory, 
				maxBackorders, truckSize, K, b, price, cost, h, theta);
	}

	@Override 
	protected void tearDown() throws Exception {
		super.tearDown();
		prob = null;
	}

/**
 * @throws Exception
 */
//	public void testMain() {
//		StochasticDemand.main(new String[]{});
//	}

	public void testGetOptimalValueFunction() throws Exception{
		InvLevel initial = new InvLevel(0);
		double vf = prob.getOptimalValueFunction().get(initial);
		double vf_cero = 126235.38374;
        assertEquals("Value Function not equal",  vf_cero, vf, 1.0E-3);	
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		junit.textui.TestRunner.run(StochasticDemandTest.class);
	}
}
