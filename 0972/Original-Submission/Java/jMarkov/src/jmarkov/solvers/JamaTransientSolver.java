/*
 * Created on 18/06/2005
 */

package jmarkov.solvers;

import jmarkov.MarkovProcess;
import jmarkov.basic.State;
import jphase.MarkovMatrix;
import Jama.Matrix;

/**
 * This class claculate the transient probabilities.
 * It uses uniformization, and basically it is a wrapper for 
 * the 'expUnif' methods in jphase.MarkovMatrix.
 * 
 * @author Julio Goéz, Germán Riaño. Universidad delos Andes
 * @see  jphase.MarkovMatrix
 * @see jmarkov.MarkovProcess
 */
public class JamaTransientSolver extends TransientSolver {

	private Matrix alpha = null;
	private State lasti0 = null;

	/**
	 * Default constructor
	 * 
	 * @param mp
	 *          the Markov Process
	 */
	public JamaTransientSolver(MarkovProcess mp) {
		super(mp);
	}

	/**
	 * @see jmarkov.solvers.Solver#label()
	 */
	@Override
	public String label() {
		return "Jama Transient";
	}

	/**
	 * Builds a unit row- vector for initial conditions
	 * 
	 * @param i0
	 *          Initial state
	 * @param n
	 *          Number of states
	 * @return vector (0,0,...,1,...,0,0)
	 */
	private Matrix getInitialConditions(State i0, int n) {
		if (i0 == lasti0)
			return alpha;
		Matrix result = new Matrix(1, n);
		result.set(0, i0.getIndex(), 1.0);
		alpha = result;
		return result;
	}

	/**
	 * @see jmarkov.solvers.TransientSolver#getTransientProbs(double,
	 *      jmarkov.basic.State)
	 * @see jphase.MarkovMatrix#expUnif(double, Matrix, Matrix)     
	 */
	@Override
	public double[] getTransientProbs(double time, State i0) {
		return getTransientProbs(new double[] {time},i0)[0];
	}

	/**
	 * @see jmarkov.solvers.TransientSolver#getTransientProbs(double[],
	 *      jmarkov.basic.State)
	 * @see jphase.MarkovMatrix#expUnif(double[], Matrix, Matrix)     
	 */
	@Override
	public double[][] getTransientProbs(double[] times, State i0) {
		MarkovMatrix Q = new MarkovMatrix(mp.getGenerator());
		int n = Q.size();
		Matrix alpha = getInitialConditions(i0, n);
		Matrix I = Matrix.identity(n, n);
		int numTimes = times.length;
		double results[][] = new double[numTimes][];
		MarkovMatrix matrixResult[] = Q.expUnif(times, alpha, I);
		for (int t=0;t<numTimes;t++){
			results[t] = matrixResult[t].getRowPackedCopy(); 
		}
		return results;
	}

	/**
	 * @see jmarkov.solvers.TransientSolver#getTransientProbs(int, double,
	 *      jmarkov.basic.State)
	 * @see jphase.MarkovMatrix#expUnif(int, double, Matrix, Matrix)     
	 */
	@Override
	public double[][] getTransientProbs(int NumberPoints, double delta, State i0) {
		double times[] = new double[NumberPoints];
		for (int t=0;t<NumberPoints;t++) times[t] = delta * t;
		return getTransientProbs(times, i0);
	}

    public String description() {
        return "Transient state solver using JAMA";
    }

}// end class
