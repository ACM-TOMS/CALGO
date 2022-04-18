/**
 * TestExamples.java
 * Created: Jul 15, 2005
 */
package jmdp;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import examples.jmdp.ControlProduccion;
import examples.jmdp.CowHerd;
import examples.jmdp.StochasticDemand;
import examples.jmdp.WagnerWhitin;

/**
 * @author Germ�n Ria�o. Universidad de los Andes. (C) 2005
 */
public class TestExamples extends TestCase {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		junit.textui.TestRunner.run(TestExamples.class);
	}

	/**
	 * @return a Test suite wit h all the examples
	 */
	public static Test suite() {
		TestSuite suite = new TestSuite("Test for jmdp.tests");
		final String[] a = new String[] {};

		suite.addTest(new TestCase("WagnerWhitin") {
			@Override
			protected void runTest() throws Throwable {
				WagnerWhitin.main(a);
			}
		});
		suite.addTest(new TestCase("StochasticDemand") {
			@Override
			protected void runTest() throws Throwable {
				StochasticDemand.main(a);
			}
		});
		suite.addTest(new TestCase("CowHerd") {
			@Override
			protected void runTest() throws Throwable {
				CowHerd.main(a);
			}
		});
		suite.addTest(new TestCase("ControlProduccion") {
			@Override
			protected void runTest() throws Throwable {
				ControlProduccion.main(a);
			}
		});
		return suite;
	}

}
