// 2014: Passes all test
package jmdp;

import jmarkov.jmdp.solvers.FiniteSolver;
import junit.framework.TestCase;
import examples.jmdp.CowHerd;

/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 *
 */
public class CowHerdTest extends TestCase {

	CowHerd prob = null;

	@SuppressWarnings("unchecked")
	@Override
	protected void setUp() throws Exception {
		double price[] = { 0, 10, 15, 25, 7, 10 };
		int mxCows[] = { 50, 75, 112, 168, 252, 378 };
		int initQ = 50;
		// CowQuantity initial = new CowQuantity(50);
		// States<CowQuantity> initSet = new
		// StatesCollection<CowQuantity>(initial);
		prob = new CowHerd(initQ, 5, price, mxCows);
		FiniteSolver<?, ?> theSolver = new FiniteSolver(prob);
		theSolver.solve();
	}

	@Override
	protected void tearDown() throws Exception {
		super.tearDown();
		prob = null;
	}

	/**
	 * @throws Exception
	 */
	public void testMain() throws Exception {
		CowHerd.main(new String[] {});
	}

	/**
	 * @throws Exception 
	 * 
	 */
	public void testGetOptimalValueFunction() throws Exception {
		//CowQuantity initial = new CowQuantity(50);
		//double vf_cero = prob.getOptimalValueFunction().get(initial);
		double vf = prob.getValue(50);
		double vf_cero = -4200.0;
		assertEquals("Value Function not equal", vf_cero, vf, 1.0E-3);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		junit.textui.TestRunner.run(CowHerdTest.class);
	}
}