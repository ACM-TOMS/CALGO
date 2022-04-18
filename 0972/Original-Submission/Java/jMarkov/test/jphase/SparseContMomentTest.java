package jphase;


import jphase.SparseContPhaseVar;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import junit.framework.TestCase;

/**
 * Set of Tests for the moment methods of the 
 * SparseContPhaseVar class
 * @author Juan F. Pérez
 */
public class SparseContMomentTest extends TestCase {
	
	/**
	 * Sparse Continuous Distribution
	 */
	SparseContPhaseVar var;
	
	/**
	 * Matrix of the Sparse Continuous Distribution
	 */
	double[][] matrix= {
			{-6,	2,	1,	2,	0,	0,	1,	0,	0,	0,	0},
			{ 1,	-5,	1,	0,	2,	0,	1,	0,	0,	0,	0},
			{ 2,	1,	-7,	0,	0,	2,	2,	0,	0,	0,	0},
			{ 2,	0,	0,	-9,	2,	1,	0,	1,	3,	0,	0},
			{ 0,	2,	0,	1,	-8,	1,	0,	1,	0,	3,	0},
			{ 0,	0,	2,	2,	1,-10,	0,	2,	0,	0,	3},
			{ 0,	0,	0,	0,	0,	0,	-2,	2,	0,	0,	0},
			{ 0,	0,	0,	0,	0,	0,	2,	-5,	0,	0,	0},
			{ 0,	0,	0,	0,	0,	0,	0,	0,	-4,	2,	1},
			{ 0,	0,	0,	0,	0,	0,	0,	0,	1,	-3,	1},
			{ 0,	0,	0,	0,	0,	0,	0,	0,	2,	1,	-5}
	};
	
	/**
	 * Vector of the Sparse Continuous Distribution
	 */
	double vector[] = {0.02, 0.04, 0.04, 0.04, 0.08,
			0.08, 0.1, 0.2, 0.04, 0.08, 0.08};

	/**
	 * 
	 * @param args Not used
	 */
	public static void main(String[] args) {
		//junit.swingui.TestRunner.run(SparseContMomentTest.class);
		junit.textui.TestRunner.run(SparseContMomentTest.class);
	}
	
	@Override
	protected void setUp() throws Exception {
		var = new SparseContPhaseVar(vector, matrix); 
	}

	@Override
	protected void tearDown() throws Exception {
		super.tearDown();
		
	}
	
	/**
	 * Test for the method jphase.SparseContPhaseVar.expectedValue
	 */
	public void testExpectedValue() {
		double calcEV = var.expectedValue();
		double realEV = 0.786212;
		assertEquals("Expected Value non equal", realEV, calcEV, 1.0E-5);
		assertTrue("Matrix changed", ((new DenseMatrix(matrix)).add(-1, var.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector)).add(-1, var.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
	}
	
	/**
	 * Test for the method jphase.SparseContPhaseVar.variance
	 */
	public void testVariance() {
		double calcVar = var.variance();
		double realVar = 0.908712;
		assertEquals("Variance non equal", realVar, calcVar, 1.0E-5);
		assertTrue("Matrix changed", ((new DenseMatrix(matrix)).add(-1, var.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector)).add(-1, var.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
	}
	
	/**
	 * Test for the method jphase.SparseContPhaseVar.CV
	 */
	public void testCV() {
		double calcCV = var.CV();
		double realCV = 1.470101;
		assertEquals("Coefficient of Variance non equal", realCV, calcCV, 1.0E-5);
		assertTrue("Matrix changed", ((new DenseMatrix(matrix)).add(-1, var.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector)).add(-1, var.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
	}
	
	/**
	 * Test for the method jphase.SparseContPhaseVar.moment
	 */
	public void testMoment2() {
		double calcM2 = var.moment(2);
		double realM2 = 1.526841;
		assertEquals("Second Moment non equal", realM2, calcM2, 1.0E-5);
		assertTrue("Matrix changed", ((new DenseMatrix(matrix)).add(-1, var.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector)).add(-1, var.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
	}
	
	/**
	 * Test for the method jphase.SparseContPhaseVar.moment
	 */
	public void testMoment3() {
		double calcM3 = var.moment(3);
		double realM3 = 4.382719;
		assertEquals("Third Moment non equal", realM3, calcM3, 1.0E-5);
		assertTrue("Matrix changed", ((new DenseMatrix(matrix)).add(-1, var.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector)).add(-1, var.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
	}
	
	/**
	 * Test for the method jphase.SparseContPhaseVar.moment
	 */
	public void testMoment4() {
		double calcM4 = var.moment(4);
		double realM4 = 16.677425;
		assertEquals("Fourth Moment non equal", realM4, calcM4, 1.0E-5);
		assertTrue("Matrix changed", ((new DenseMatrix(matrix)).add(-1, var.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector)).add(-1, var.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
	}
}