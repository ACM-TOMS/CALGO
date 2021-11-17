/**
 * QueueMM1NTest.java
 * Created: Jul 6, 2005
 */
package jmarkov;

import static jmarkov.Utils.assertArrayEquals;
import jmarkov.MarkovProcess.Status;
import junit.framework.TestCase;
import examples.jmarkov.QueueMM1N;

/**
 * @author Germán Riaño. Universidad de los Andes. (C) 2005
 */
public class QueueMM1NTest extends TestCase {

    private QueueMM1N model = new QueueMM1N(4.0, 2.0, 18);
    private double prob[] = null;
    private int N;

    /**
     * @param args
     */
    public static void main(String[] args) {
        //junit.swingui.TestRunner.run(QueueMM1NTest.class);
        junit.textui.TestRunner.run(QueueMM1NTest.class);
    }

    /**
     * Test method for 'jmarkov.SimpleMarkovProcess.generate()'
     */
    public void testGenerate() {
        model.generate();
        assertTrue("QueueMM1N not generates",
                model.getStatus() == Status.GENERATED);
    }

    /**
     * @see junit.framework.TestCase#setUp()
     */
    @Override
    protected void setUp() throws Exception {
        // let's compute probs using B&D
        double lam = model.getLambda();
        double mu = model.getMu();
        N = model.getMax();
        prob = new double[N + 1];
        double sum = prob[0] = 1.0;
        for (int i = 1; i <= N; i++) {
            prob[i] = prob[i - 1] * lam / mu;
            sum += prob[i];
        }
        for (int i = 0; i <= N; i++) {
            prob[i] = prob[i] / sum;
        }
        model.setDebugLevel(0);
    }

    /**
     * @see junit.framework.TestCase#tearDown()
     */
    @Override
    public void tearDown() {
        model = null;
    }

    /*
     * Tests method for 'jmarkov.SimpleMarkovProcess.getSteadyState()'
     */
    /**
     * @throws Exception
     */
    public void testGetSteadyState() throws Exception {

        assertArrayEquals("Steady probabilities not equal in model MM1N", prob,
                model.getSteadyState(), 1e-5);
    }

    /**
     * @throws Exception
     */
    /*
     * Test method for 'jmarkov.SimpleMarkovProcess.getMOPsAvg()'
     */
    public void testGetMOPsAvg() throws Exception {
        double ldaEff = model.getLambda() * (1 - prob[model.getMax()]);
        double sum = 0.0, sumQ = 0.0;
        for (int i = 1; i <= N; i++) {
            sum += i * prob[i];
            sumQ += (i - 1) * prob[i];
        }
        // Number in the system";
        double L = sum;
        // "Number in queue", ;
        double Lq = sumQ;
        // "Server Utilization"
        double rho = ldaEff / model.getMu();

        // boolean equal = equalArrays(mops, new double[] {
        // 7.017612574190854,
        // 6.019569513613666, 0.9980430605771874 }, 1e-4);
        double mops[] = model.getMOPsAvg();

        assertArrayEquals("MOPS not equal: ", new double[] { L, Lq, rho },
                mops, 1e-5);

    }

    /**
     * @throws Exception
     */
    /*
     * Test method for 'jmarkov.SimpleMarkovProcess.getEventsRates()'
     */
    public void testGetEventRates() throws Exception {
        double ldaEff = model.getLambda() * (1 - prob[model.getMax()]);
        assertEquals("Arrival rate not equal: ", ldaEff, model.getEventRate(0),
                1e-5);
    }

}
