package jphase;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import jphase.SparseContPhaseVar;
import jphase.SparseDiscPhaseVar;
import junit.framework.TestCase;

/**
 * Set of Tests for the closure methods of the 
 * SparseContPhaseVar class
 * @author Juan F. Pérez
 */
public class SparseContClosureTest extends TestCase {

	/**
	 * First Continuous Distribution
	 */
	SparseContPhaseVar var1;
	
	/**
	 * First Continuous Distribution
	 */
	SparseContPhaseVar var2;
	
	/**
	 * Discrete distribution
	 */
	SparseDiscPhaseVar varD;
	
	/**
	 * Matrix for the first continuous distribution
	 */
	double[][] matrix1= {
			{	-2,	2},	
			{	2, -5}};
	/**
	 * Matrix for the second continuous distribution
	 */
	double[][] matrix2= {
			{-4,	2,	1},
			{1,	-3,	1},
			{2,	1,	-5}};
	/**
	 * Matrix for discrete distribution 
	 */
	double[][] matrixD= {
			{0.76,	0.24},
			{	0,	0.76}};
	
	/**
	 * Vector for the first continuous distribution 
	 */
	double[] vector1 = {0.2,	0.4};		
	
	/**
	 * Vector for the second continuous distribution 
	 */
	double[] vector2 = {0.1,	0.2,	0.2};	
	
	/**
	 * Vector for discrete distribution 
	 */
	double[] vectorD = {0.15, 0.85};
	
	/**
	 * Main method
	 * @param args
	 */
	public static void main(String[] args) {
		//junit.swingui.TestRunner.run(SparseContClosureTest.class);
		junit.textui.TestRunner.run(SparseContClosureTest.class);
	}
	
	@Override
	protected void setUp() throws Exception {
		var1 = new SparseContPhaseVar(vector1, matrix1);
		var2 = new SparseContPhaseVar(vector2, matrix2);
		varD = new SparseDiscPhaseVar(vectorD, matrixD);
	}

	@Override
	protected void tearDown() throws Exception {
		super.tearDown();
		
	}

	/**
	 * Test for the method jphase.SparseContPhaseVar.Sum
	 */
	public void testSum() {
        SparseContPhaseVar calcSum = (SparseContPhaseVar)var1.sum(var2);
		double[] vecRes = {0.2,	0.4, 0.04, 0.08, 0.08};
		double[][] matRes ={
				{-2, 2, 0, 0, 0},
				{2, -5, 0.3, 0.6, 0.6},
				{0,	0, -4, 2, 1},
				{0,	0, 1, -3, 1},
				{0,	0, 2, 1, -5}
		};
		SparseContPhaseVar realSum = new SparseContPhaseVar(vecRes, matRes);
		assertTrue("Sum of Variables non equal (Matrix)", (realSum.getMatrix().add(-1, calcSum.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Sum of Variables non equal (Vector)", (realSum.getVector().add(-1, calcSum.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
		assertTrue("Matrix changed", ((new DenseMatrix(matrix1)).add(-1, var1.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector1)).add(-1, var1.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
		assertTrue("Matrix changed", ((new DenseMatrix(matrix2)).add(-1, var2.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector2)).add(-1, var2.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
		
	}

	/**
	 * Test for the method jphase.SparseContPhaseVar.Min
	 */
	public void testMin() {
        SparseContPhaseVar calcSum = (SparseContPhaseVar)var1.min(var2);
		double[] vecRes = {0.02, 0.04, 0.04, 0.04, 0.08, 0.08};
		double[][] matRes ={
				{-6, 2,	1,	2,	0,	0},
				{1,	-5,	1,	0,	2,	0},
				{2,	1,	-7,	0,	0,	2},
				{2,	0,	0,	-9,	2,	1},
				{0,	2,	0,	1,	-8,	1},
				{0,	0,	2,	2,	1,	-10}
		};
		SparseContPhaseVar realSum = new SparseContPhaseVar(vecRes, matRes);
		assertTrue("Min of Variables non equal (Matrix)", (realSum.getMatrix().add(-1, calcSum.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Min of Variables non equal (Vector)", (realSum.getVector().add(-1, calcSum.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
		assertTrue("Matrix changed", ((new DenseMatrix(matrix1)).add(-1, var1.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector1)).add(-1, var1.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
		assertTrue("Matrix changed", ((new DenseMatrix(matrix2)).add(-1, var2.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector2)).add(-1, var2.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
	}
	
	/**
	 * Test for the method jphase.SparseContPhaseVar.Max
	 */
	public void testMax() {
        SparseContPhaseVar calcSum = (SparseContPhaseVar)var1.max(var2);
		double[] vecRes = {0.02, 0.04, 0.04, 0.04, 0.08, 0.08, 0.1, 0.2, 0.04, 0.08, 0.08};

		double[][] matRes ={
				{-6, 2,	1,	2,	0,	0,	1,	0,	0,	0,	0},
				{1,	-5,	1,	0,	2,	0,	1,	0,	0,	0,	0},
				{2,	1,	-7,	0,	0,	2,	2,	0,	0,	0,	0},
				{2,	0,	0, -9,	2,	1,	0,	1,	3,	0,	0},
				{0,	2,	0,	1,	-8,	1,	0,	1,	0,	3,	0},
				{0,	0,	2,	2,	1,	-10,0,	2,	0,	0,	3},
				{0,	0,	0,	0,	0,	0,	-2,	2,	0,	0,	0},
				{0,	0,	0,	0,	0,	0,	2,	-5,	0,	0,	0},
				{0,	0,	0,	0,	0,	0,	0,	0,	-4,	2,	1},
				{0,	0,	0,	0,	0,	0,	0,	0,	1,	-3,	1},
				{0,	0,	0,	0,	0,	0,	0,	0,	2,	1,	-5}

		};
		SparseContPhaseVar realSum = new SparseContPhaseVar(vecRes, matRes);
		assertTrue("Max of Variables non equal (Matrix)", (realSum.getMatrix().add(-1, calcSum.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Max of Variables non equal (Vector)", (realSum.getVector().add(-1, calcSum.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
		assertTrue("Matrix changed", ((new DenseMatrix(matrix1)).add(-1, var1.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector1)).add(-1, var1.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
		assertTrue("Matrix changed", ((new DenseMatrix(matrix2)).add(-1, var2.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector2)).add(-1, var2.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
	}
	
	/**
	 * Test for the method jphase.SparseContPhaseVar.Mix
	 */
	public void testMix() {
        SparseContPhaseVar calc = (SparseContPhaseVar)var1.mix(0.2, var2);
		double[] vecRes = {0.04, 0.08, 0.08, 0.16, 0.16};
		double[][] matRes ={
				{-2, 2, 0, 0, 0},
				{2, -5, 0, 0, 0},
				{0,	0, -4, 2, 1},
				{0,	0, 1, -3, 1},
				{0,	0, 2, 1, -5}
		};
		SparseContPhaseVar real = new SparseContPhaseVar(vecRes, matRes);
		assertTrue("Mix of Variables non equal (Matrix)", (real.getMatrix().add(-1, calc.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Mix of Variables non equal (Vector)", (real.getVector().add(-1, calc.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
		assertTrue("Matrix changed", ((new DenseMatrix(matrix1)).add(-1, var1.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector1)).add(-1, var1.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
		assertTrue("Matrix changed", ((new DenseMatrix(matrix2)).add(-1, var2.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector2)).add(-1, var2.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
	}
	
	/**
	 * Test for the method jphase.SparseContPhaseVar.SumGeom
	 */
	public void testSumGeom() {
        SparseContPhaseVar calc = (SparseContPhaseVar)var2.sumGeom( 0.2);
		double[] vecRes = {0.1, 0.2, 0.2};
		double[][] matRes ={
				{-3.92,	2.16,	1.16},
				{1.08,	-2.84,	1.16},
				{2.16,	1.32,	-4.68}
		};
		
		SparseContPhaseVar real = new SparseContPhaseVar(vecRes, matRes);
		assertTrue("Geometric Sum of Variables non equal (Matrix)", (real.getMatrix().add(-1, calc.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Geometric Sum of Variables non equal (Vector)", (real.getVector().add(-1, calc.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
		assertTrue("Matrix changed", ((new DenseMatrix(matrix2)).add(-1, var2.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector2)).add(-1, var2.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
	}
	
	/**
	 * Test for the method jphase.SparseContPhaseVar.residualTime
	 */
	public void testResidual() {
        SparseContPhaseVar calc = (SparseContPhaseVar)var2.residualTime(2);
		double[] vecRes = {0.318884557,	0.472382406, 0.208733037};
		SparseContPhaseVar real = new SparseContPhaseVar(vecRes, var2.getMatrixArray());
		assertTrue("Residual Time Distribution non equal (Vector)", (real.getVector().add(-1, calc.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
		assertTrue("Matrix changed", ((new DenseMatrix(matrix2)).add(-1, var2.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector2)).add(-1, var2.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
	}
	
	/**
	 * Test for the method jphase.SparseContPhaseVar.eqResidualTime
	 */
	public void testEqResidual() {
		//int n2 = var2.getNumPhases();
		//SparseContPhaseVar calc = new SparseContPhaseVar(n2);
        SparseContPhaseVar calc = (SparseContPhaseVar)var2.eqResidualTime();
		double[] vecRes = {0.3,	0.45,	0.25};
		
		
		SparseContPhaseVar real = new SparseContPhaseVar(vecRes, var2.getMatrixArray());
		//assertTrue("Geometric Sum of Variables non equal (Matrix)", (real.getMatrix().add(-1, calc.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Equilibrium Residual Time Distribution non equal (Vector)", (real.getVector().add(-1, calc.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
		assertTrue("Matrix changed", ((new DenseMatrix(matrix2)).add(-1, var2.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector2)).add(-1, var2.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
	}
	
	/**
	 * Test for the method jphase.SparseContPhaseVar.waitingQ
	 */
	public void testWaitingQ() {
        SparseContPhaseVar calc = (SparseContPhaseVar)var2.waitingQ(0.5);
		double[] vecRes = {0.15, 0.225, 0.125};
		double[][] matRes ={
				{-3.85,	2.225,	1.125},
				{1.15,	-2.775,	1.125},
				{2.3,	1.45,	-4.75}
		};
		SparseContPhaseVar real = new SparseContPhaseVar(vecRes, matRes);
		assertTrue("Waiting in Queue Distribution non equal (Matrix)", (real.getMatrix().add(-1, calc.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-4);
		assertTrue("Waiting in Queue Distribution non equal (Vector)", (real.getVector().add(-1, calc.getVector())).norm(Vector.Norm.Infinity)<1.0E-4  );
		assertTrue("Matrix changed", ((new DenseMatrix(matrix2)).add(-1, var2.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector2)).add(-1, var2.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
	}
	
	/**
	 * Test for the method jphase.SparseContPhaseVar.SumPH
	 */
	public void testSumPH() {
        SparseContPhaseVar calc = (SparseContPhaseVar)var2.sumPH(varD);
		double[] vecRes = {0.0242,	0.1418,	0.0484,	0.2836,	0.0484,	0.2836};
		double[][] matRes ={{-3.9387,	0.0312,	2.1226,	0.0624,	1.1226,	0.0624},
							{0.0000, -3.9387,	0.0000,	2.1226,	0.0000,	1.1226},
							{1.0613,	0.0312,	-2.8774,	0.0624,	1.1226,	0.0624},
							{0.0000,	1.0613,	0.0000,	-2.8774,	0.0000,	1.1226},
							{2.1226,	0.0624,	1.2452,	0.1249,	-4.7548,	0.1249},
							{0.0000,	2.1226,	0.0000,	1.2452,	0.0000,	-4.7548}
							};

		SparseContPhaseVar real = new SparseContPhaseVar(vecRes, matRes);
		assertTrue("Waiting in Queue Distribution non equal (Matrix)", (real.getMatrix().add(-1, calc.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-4);
		assertTrue("Waiting in Queue Distribution non equal (Vector)", (real.getVector().add(-1, calc.getVector())).norm(Vector.Norm.Infinity)<1.0E-4  );
		assertTrue("Matrix changed", ((new DenseMatrix(matrix2)).add(-1, var2.getMatrix())).norm(Matrix.Norm.Maxvalue)<1.0E-5);
		assertTrue("Vector changed", ((new DenseVector(vector2)).add(-1, var2.getVector())).norm(Vector.Norm.Infinity)<1.0E-5  );
	}	
}