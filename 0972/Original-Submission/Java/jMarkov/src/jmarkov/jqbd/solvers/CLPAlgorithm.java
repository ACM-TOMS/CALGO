/**
 * MtjLogRedSolver.java
 * Created: 2/09/2006
 */
package jmarkov.jqbd.solvers;

import jphase.PhaseVar;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;

/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 */
public class CLPAlgorithm extends QBDPhaseSolver {

	/**
	 * 
	 */
	public CLPAlgorithm(PhaseVar arrivals, PhaseVar service) {

		initilicePHPH1(arrivals, service);
	}

	@Override
	public double[][] getRmatrix() {


		int dimen = MA1.numRows();


		Matrix I = new DenseMatrix(Matrices.identity(dimen));

		DenseMatrix mA1I = new DenseMatrix(dimen, dimen);
		G = new DenseMatrix(dimen, dimen);
		U = new DenseMatrix(dimen, dimen);
		N = new DenseMatrix(dimen, dimen);
		DenseMatrix UA = new DenseMatrix(dimen, dimen);
		R = new DenseMatrix(dimen, dimen);

		U.set(MA1);
		

		U.mult(-1, I, UA);// UA = -U
		UA.solve(I, mA1I);// ma1I=(-U)^{-1}
		mA1I.mult(MA2, G);// G=(-U)^{-1}*MA2.

		double compare = 1;
		double epsilon = 1e-1;

		while (compare > epsilon) {
			U.set(MA1);
			MA0.multAdd(G, U); //U = A1 + MA0*G
			
			mA1I = new DenseMatrix(dimen, dimen);
			U.mult(-1, I, UA);// U=A_1+A_0*G
			UA.solve(I, mA1I);// (I-U)^{-1}
			N.set(mA1I);// (I-U)^{-1}
			mA1I.mult(MA2, G);// G = (I-U)^-1*A2

			// Procedimiento para el calculo de la norma infinito
			// ||1-G1||_{\infty}

			double[][] ones = new double[G.numRows()][1];
			for (int i = 0; i < G.numRows(); i++) {
				ones[i][0] = 1;
			}

			DenseMatrix one = new DenseMatrix(ones);

			DenseMatrix check = new DenseMatrix(G.numRows(), 1);

			G.mult(-1, one, check);
			check.add(one);

			compare = Math.abs(check.get(0, 0));

			for (int i = 1; i < check.numRows(); i++) {
				compare = Math.max(compare, Math.abs(check.get(i, 0)));
			}
		}
		/*
		U = new DenseMatrix(dimen, dimen);
		U = (DenseMatrix) MA0.multAdd(G, MA1.copy());*/
		U.set(MA1);
		MA0.multAdd(G, U);
		mA1I = new DenseMatrix(dimen, dimen);
		U.mult(-1, I, UA);// U=A_1+A_0*G
		UA.solve(I, mA1I);// (I-U)^{-1}
		N.set(mA1I);// (I-U)^{-1}
		MA0.mult(mA1I, R);// R=A_0*H1I.

		return Matrices.getArray(R);
	}

	public String label() {
		return "Continuous Linear Progresion algorithm solver.";
	}
}

