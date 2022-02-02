/*
 * Created on 08-jul-2005
 */
package jmarkov;

import static jmarkov.Utils.assertArrayEquals;
import jmarkov.basic.exceptions.NotUnichainException;
import junit.framework.TestCase;
import examples.jmarkov.Kanban;
/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 *
 */
public class KanbanTest extends TestCase {
    private Kanban kaby;

    /**
     * @param args
     */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(KanbanTest.class);
    }

    /**
     * @param name
     */
    public KanbanTest(String name) {
        super(name);
        kaby = new Kanban();
        kaby.setDebugLevel(0);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    /**
     * 
     */
    /*
     * Test method for 'jmarkov.SimpleMarkovProcess.generate()'
     */
    public final void testGenerate() {
        kaby.generate();
        assertTrue(kaby.isGenerated());

    }

    /**
     * Test method for 'jmarkov.SimpleMarkovProcess.getMOPsAvg()'
     * @throws Exception 
     */
    public final void testGetMOPsAvg() throws Exception {
        double mops[] = kaby.getMOPsAvg();
        double expMops[] = { 0.7003166576899312, 0.7003166613086502,
                0.7003166718329953, 0.6714510596567259, 0.6634469869009703,
                1.6282322852477127, 0.6362363535199593 };
        assertArrayEquals("Error in MOPS", expMops, mops, 1e-4);

    }

    /**
     * Test method for 'jmarkov.SimpleMarkovProcess.getEventRates()'
     * @throws NotUnichainException 
     */
    public final void testGetRates() throws NotUnichainException {
        double ratesAvg[] = kaby.getEventsRates();
        double expRates[] = { 1.4006333153798625, 1.4006333226173004,
                1.4006333436659906 };
        assertArrayEquals("Error in rates.",expRates,  ratesAvg, 1e-4);
    }

}
