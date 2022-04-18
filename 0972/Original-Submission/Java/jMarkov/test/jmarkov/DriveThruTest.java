/**
 * DriveThruTest.java
 * Created: Jul 6, 2005
 */
package jmarkov;

import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.exceptions.NotUnichainException;
import junit.framework.TestCase;
import examples.jmarkov.DriveThru;

/**
 * @author Germán Riaño. Universidad de los Andes. (C) 2005
 */
public class DriveThruTest extends TestCase {

    private DriveThru theDT = null;

    /**
     * @param args
     */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(DriveThruTest.class);
    }

    /*
     * @see TestCase#setUp()
     */
    @Override
    protected void setUp() throws Exception {
        super.setUp();
        theDT = new DriveThru(80.0, 12.0, 30.0, 4, 2, 1);
        theDT.setDebugLevel(0);
    }

    /*
     * @see TestCase#tearDown()
     */
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
        theDT = null;
    }

    /*
     * Test method for 'jmarkov.SimpleMarkovProcess.generate()'
     */
    /**
     * 
     */
    public final void testGenerate() {
        theDT.generate();
        assertTrue("Drive Thru not generated correctly",
                theDT.getStatus() == SimpleMarkovProcess.Status.GENERATED);
    }

    /*
     * Test method for 'jmarkov.SimpleMarkovProcess.getMOPsAvg()'
     */
    /**
     * @throws NotUnichainException 
     * 
     */
    public final void testGetMOPsAvg() throws NotUnichainException {
        double mops[] = theDT.getMOPsAvg();
        double expMops[] = new double[] { 1.8535673938908643,
                0.9085921200048781, 0.36343683564258694, 1.272028955647465,
                0.0, 0.04537381806445134, 0.04537381806445134,
                1.8989412119553157, 3.17097016760278 };
        Utils.assertArrayEquals("MOPS incorrect in DriveThru", expMops, mops,
                1e-5);

    }

}
