/*
 * Created on 06-jul-2005
 */
package jmarkov;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 *
 */
public class AllTests {

    /**
     * @param args
     */
    public static void main(String[] args) {
        //junit.swingui.TestRunner.run(AllTests.class); junit.swing deprecated
        junit.textui.TestRunner.runAndWait(suite());
    }

    /**
     * @return The suite wit all examples
     */
    public static Test suite() {
        TestSuite suite = new TestSuite("Test for jmarkov.tests");
        //$JUnit-BEGIN$
        suite.addTestSuite(JacksonTest.class);
        suite.addTestSuite(KanbanTest.class);
        suite.addTestSuite(QueueMM1NTest.class);
        suite.addTestSuite(DriveThruTest.class);
        suite.addTestSuite(QueueMH2k1Test.class);
        suite.addTestSuite(TransientTest.class);
        suite.addTestSuite(BucketBrigadesTest.class);
        //$JUnit-END$
        //suite.addTest(SolversTest.suite());
        return suite;
    }

}
