package jphase.fit;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import jphase.PhaseVar;
import jphase.DenseContPhaseVar;
import jphase.fit.MomentsECCompleteFit;
import junit.framework.TestCase;

/**
 * Set of Tests for the fitting method in the ECCompleteFit class
 * (Erlang Coxian Distribution)
 * @author Juan F. Pérez
 */
public class ECCompleteFitTest extends TestCase {
	
	/**
	 * First Moment
	 */
	double m1;
	
	/**
	 * Second Moment
	 */
	double m2;
	
	/**
	 * Third Moment
	 */
	double m3;
	
	/**
	 * Fitter object
	 */
	MomentsECCompleteFit fitter;
	
	/**
	 * Fitted Variable
	 */
	PhaseVar calcVar;

	/**
	 * 
	 * @param args Not used
	 */
	public static void main(String[] args) {
		junit.textui.TestRunner.run(ECCompleteFitTest.class);
	}
	@Override
	protected void setUp() throws Exception {
		//fitter = new MomentsECCompleteFit();
	}

	@Override
	protected void tearDown() throws Exception {
		super.tearDown();
		
	}
	
	/**
	 * Test for the moment fitter class 
	 */
	public void testParamCaso1() {
		m1 = 1;
		m2 = 2;
		m3 = 8;
		fitter = new MomentsECCompleteFit(m1, m2, m3);
		calcVar = fitter.fit();
		double[] vecRes = {0.875,	0.0, 0.0};
		double[][] matRes ={
				{-1.1667,	1.1667,	0.0000},
				{0.0000,	-3.6904,	0.0093},
				{0.0000,	0.0000,	-0.1717}
		};
		
		DenseContPhaseVar realVar = new DenseContPhaseVar(vecRes, matRes);
		assertTrue("Variable not well fitted from moments (Matrix)", (realVar.getMatrix().add(-1, calcVar.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-3);
		assertTrue("Variable not well fitted from moments (Vector)", (realVar.getVector().add(-1, calcVar.getVector())).norm(Vector.Norm.Infinity)<1.0E-3  );
	}
	
	/**
	 * Test for the moment fitter class 
	 */
	public void testParamCaso2() {
		m1 = 1;
		m2 = 2;
		m3 = 5;
		fitter = new MomentsECCompleteFit(m1, m2, m3);
		calcVar = fitter.fit();
		double[] vecRes = {0.6667,	0.0, 0.0};
		double[][] matRes ={
				{-2,	2,	0},
				{0,	-2,	2},
				{0,	0,	-2}
		};
		
		DenseContPhaseVar realVar = new DenseContPhaseVar(vecRes, matRes);
		assertTrue("Variable not well fitted from moments (Matrix)", (realVar.getMatrix().add(-1, calcVar.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-3);
		assertTrue("Variable not well fitted from moments (Vector)", (realVar.getVector().add(-1, calcVar.getVector())).norm(Vector.Norm.Infinity)<1.0E-3  );
	}
	
	/**
	 * Test for the moment fitter class 
	 */
	public void testParamCaso3() {
		m1 = 1.53;
		m2 = 10.94;
		m3 = 148.57;
		fitter = new MomentsECCompleteFit(m1, m2, m3);
		calcVar = fitter.fit();
		double[] vecRes = {1.0,	0.0};
		double[][] matRes ={
				{-1.845307614,	0.388202144},
				{0.00000,		-0.212909421}
		};

		DenseContPhaseVar realVar = new DenseContPhaseVar(vecRes, matRes);
		assertTrue("Variable non well fitted from moments (Matrix)", (realVar.getMatrix().add(-1, calcVar.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Variable non well fitted from moments (Vector)", (realVar.getVector().add(-1, calcVar.getVector())).norm(Vector.Norm.Infinity)<1.0E-5);
	}
	
	/**
	 * Test for the moment fitter class 
	 */
	public void testParamCaso4() {
		m1 = 1;
		m2 = 2;
		m3 = 4.571358;
		fitter = new MomentsECCompleteFit(m1, m2, m3);
		calcVar = fitter.fit();
		double[] vecRes = {0.583321346,	0,	0,	0,	0,	0};
		double[][] matRes ={
				{-3.500431608,	3.500431608,	0,	0,	0,	0},
				{0,	-3.500431608,	3.500431608,	0,	0,	0},
				{0,	0,	-3.500431608,	3.500431608,	0,	0},
				{0,	0,	0,	-3.500431608,	3.500431608,	0},
				{0,	0,	0,	0,	-3.500431608,	3.500431608},
				{0,	0,	0,	0,	0,	-3.497412587}
		};
		DenseContPhaseVar realVar = new DenseContPhaseVar(vecRes, matRes);
		System.out.println("Fitted variable: "+calcVar.toString());
		assertTrue("Variable non well fitted from moments (Matrix)", (realVar.getMatrix().add(-1, calcVar.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Variable non well fitted from moments (Vector)", (realVar.getVector().add(-1, calcVar.getVector())).norm(Vector.Norm.Infinity)<1.0E-5);
	}
}