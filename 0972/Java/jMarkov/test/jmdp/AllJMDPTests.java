/**
 * AllTests.java
 * Created: Jul 17, 2005
 */
package jmdp;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * @author Germán Riaño. Universidad de los Andes. (C) 2005
 *
 */
public class AllJMDPTests {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		junit.textui.TestRunner.runAndWait(suite());
	}

	/**
	 * @return the Test
	 */
	public static Test suite() {
		TestSuite suite = new TestSuite("Test for tests");
//		$JUnit-BEGIN$
        suite.addTestSuite(SolversTest.class);
		suite.addTest(TestExamples.suite());
		//$JUnit-END$
		return suite;
	}

}
