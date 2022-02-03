package jphase;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * This Test includes all the test of the classes included
 * in the jphase package
 * @author Juan F. Pérez
 *
 */
public class AllTests {

	/**
	 * Test that includes all the tests for the jphase package
	 * @return All the tests for the jphase package
	 */
	public static Test suite() {
		TestSuite suite = new TestSuite("Test for jphase.tests");
		//$JUnit-BEGIN$
		suite.addTestSuite(DenseContClosureTest.class);
		suite.addTestSuite(DenseContMomentTest.class);
		suite.addTestSuite(DenseContProbTest.class);
		suite.addTestSuite(SparseContMomentTest.class);
		suite.addTestSuite(SparseContClosureTest.class);
		suite.addTestSuite(SparseContProbTest.class);
		//$JUnit-END$
		return suite;
	}

}
