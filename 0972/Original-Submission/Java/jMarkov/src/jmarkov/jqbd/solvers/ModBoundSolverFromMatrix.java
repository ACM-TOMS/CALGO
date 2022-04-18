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
public class ModBoundSolverFromMatrix extends QBDPhaseSolver {

	/**
	 * 
	 */
	public ModBoundSolverFromMatrix(PhaseVar arrivals, PhaseVar service) {

		initilicePHPH1(arrivals, service);
	}

	@Override
	public double[][] getRmatrix() {

		double epsilon = 1e-108;

		int dimen = MA1.numRows();

		Matrix I = new DenseMatrix(Matrices.identity(dimen));

		DenseMatrix INV = new DenseMatrix(dimen,dimen);
		Matrix sum = new DenseMatrix(I.copy());
		sum.add(MA0.copy().scale(-1));
		sum.add(MA1.copy().scale(-1));
		sum.solve(I, INV);	//MA1IINV = (I-A1-A0)^-1
		G = new DenseMatrix(dimen, dimen);
		DenseMatrix Gold = new DenseMatrix(dimen, dimen);
		U = new DenseMatrix(dimen, dimen);

		INV.mult(MA2, G);
		double compare = 1;

		while (compare > epsilon) {

			Gold = G.copy();
			MA0.mult(G, U);// U=HL
			U.add(MA1);// U=A1 + A0G

			sum = U.copy().scale(-1);
			sum.solve(I, INV);	//MA1IINV = (-U)^-1
			INV.mult(MA2, G);

			// Procedimiento para el calculo de la norma infinito
			// ||G-Gold||_{\infty}

			Gold.add(G.copy().scale(-1));

			compare = Math.abs(Gold.get(0, 0));

			for (int i = 0; i < Gold.numRows(); i++) {
				for (int j = 0; j < Gold.numColumns(); j++) {
					compare = Math.max(compare, Math.abs(Gold.get(i, j)));
				}
			}

			System.out.println("Compare: " + compare);

		}

		U = new DenseMatrix(dimen, dimen);
		N = new DenseMatrix(dimen, dimen);
		R = new DenseMatrix(dimen, dimen);
		U = (DenseMatrix) MA0.multAdd(G, MA1.copy());
		DenseMatrix mA1I = new DenseMatrix(dimen, dimen);
		DenseMatrix UA = new DenseMatrix(dimen, dimen);

		U.mult(-1, I, UA);// U=A_1+A_0*G
		UA.solve(I, mA1I);// (I-U)^{-1}
		N.set(mA1I);
		MA0.mult(mA1I, R);// R=A_0*H1I.

		return Matrices.getArray(R);
	}

	public String label() {
		return "Modified boundary solver. This solver uses the MTJ p"
		+ "ackage to handle matrices from a phase ditribution.";
	}

}

