/*
 * Created on 12-jul-2005
 */
package jmarkov;

import static jmarkov.Utils.assertArrayEquals;
import jmarkov.basic.exceptions.NotUnichainException;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrix;
import examples.jmarkov.QueueMH2k1;

/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 */
public class QueueMH2k1Test extends TestCase {
    private QueueMH2k1 model = null;
    double[][] matrixR = { { 0.6419753084348551, 0.44444444430105357 },//
            { 0.19753086409800832, 0.44444444437554426 } };

    private double[] arrayA0 = { 4.0, 0.0, 0.0, 4.0 };
    private double[] arrayA1 = { -13.0, 0.0, 9.0, -13.0 };
    private double[] arrayA2 = { 0.0, 9.0, 0.0, 0.0 };
    private double[] arrayB00 = { -7.0 };
    private double[] arrayB01 = { 7.0, 0.0 };
    private double[] arrayB10 = { 0.0, 9.0 };
    private double[] pi0 = { 0.0666858472473757 };
    private double[] pi1 = { 0.074914784737818, 0.0583356184930964 };
    private double[] pi1mod = { 0.4666666665986373, 0.4666666665986373 };
    private double[] mopsAvg = { 6.533333320037249, 0.933314152752624,
            0.466657099381307, 0.466657053371318 };
    private double[] eventRates = { 0.46680093073163, 3.7332566110105,//
            4.19991389443176, 4.19991348034186 };

    /**
     * @param args
     */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(QueueMH2k1Test.class);
    }

    /**
     * @param name
     */
    public QueueMH2k1Test(String name) {
        super(name);

    }

    @Override
    protected void setUp() throws Exception {
        model = new QueueMH2k1(7, 4, 9, 9);
        model.generate();
        model.setDebugLevel(0);
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
        model = null;
    }

    /**
     * Test method for 'jmarkov.GeomProcess.getAMatrices()'
     */

    public void testGetAMatrices() {
        Matrix A[] = model.getAMatrices();
        double A0[] = (new DenseMatrix(A[0])).getData();
        double A1[] = (new DenseMatrix(A[1])).getData();
        double A2[] = (new DenseMatrix(A[2])).getData();
        assertArrayEquals("A0 different", arrayA0, A0, 1e-4);
        assertArrayEquals("A1 different", arrayA1, A1, 1e-4);
        assertArrayEquals("A2 different", arrayA2, A2, 1e-4);

    }

    /**
     * 
     */
    public void testIsStable() {
        boolean stable = model.isStable();
        assertTrue("Stability test failed", stable);
    }

    /**
     * Test method for 'jmarkov.GeomProcess.getRmatrix()'
     * @throws NotUnichainException 
     */
    public void testGetRmatrix() throws NotUnichainException {
        double[][] newR = model.matrixRtoArray();
        assertArrayEquals("matrix R is not equal", matrixR, newR, 1e-1);
    }

    /**
     * Test method for 'jmarkov.GeomProcess.getBMatrices()'
     */
    public void testGetBMatrices() {
        Matrix B[] = model.getBMatrices();
        double B00[] = (new DenseMatrix(B[0])).getData();
        double B01[] = (new DenseMatrix(B[1])).getData();
        double B10[] = (new DenseMatrix(B[2])).getData();
        assertArrayEquals("B00 different", arrayB00, B00, 1e-4);
        assertArrayEquals("B01 different", arrayB01, B01, 1e-4);
        assertArrayEquals("B10 different", arrayB10, B10, 1e-4);
    }

    /**
     * Test method for 'jmarkov.GeomProcess.getVectorPi1()'
     * @throws NotUnichainException 
     */
    public void testGetVectorPi1() throws NotUnichainException {
        double vectorPi1[] = model.getVectorPi1();
        assertArrayEquals("Pi1 vectors are not equal", vectorPi1, pi1, 1e-2);
    }

    /**
     * Test method for 'jmarkov.GeomProcess.getVectorPi1Mod()'
     * @throws NotUnichainException 
     */
    public void testGetVectorPi1Mod() throws NotUnichainException {
        double vectorPi1mod[] = model.getVectorPi1Mod();
        assertArrayEquals("Pi1 vectors are not equal", vectorPi1mod, pi1mod,
                1e-2);

    }

    /**
     * Test method for 'jmarkov.GeomProcess.getVectorPi0()'
     * @throws NotUnichainException 
     */
    public void testGetVectorPi0() throws NotUnichainException {
        double vectorPi0[] = model.getVectorPi0();
        assertArrayEquals("Pi0 vectors are not equal", vectorPi0, pi0, 1e-1);
    }

    /**
     * Test method for 'jmarkov.SimpleMarkovProcess.getMOPsAvg()'
     * @throws Exception
     */
    public void testGetMOPsAvg() throws Exception {
        double averageMOPs[] = model.getMOPsAvg();
        assertArrayEquals("The MOPS average are not equal", mopsAvg,
                averageMOPs, 1e-1);
    }

    /**
     * Test method for 'jmarkov.SimpleMarkovProcess.getEventsRates()'
     * @throws Exception
     */
    public void testGetEventsRates() throws Exception {
        double rates[] = model.getEventsRates();
        assertArrayEquals("The rates are not equal", eventRates, rates, 1e-1);
    }

}
