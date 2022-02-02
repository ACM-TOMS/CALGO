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
public class MtjLogRedSolverFromMatrix extends QBDPhaseSolver {
    
    /**
     */
    public MtjLogRedSolverFromMatrix(PhaseVar arrivals, PhaseVar service) {
		initilicePHPH1(arrivals, service);
	}
	
	public MtjLogRedSolverFromMatrix (double[][] xMB0, double[][] xMB1, double[][] xMB2, double[][] xMA0, double[][] xMA1, double[][] xMA2){
		MA0 = new DenseMatrix(xMA0);
		MA1 = new DenseMatrix(xMA1);
		MA2 = new DenseMatrix(xMA2);
		MB0 = new DenseMatrix(xMB0);
		MB1 = new DenseMatrix(xMB1);
		MB2 = new DenseMatrix(xMB2);
	}

    @Override
	public double[][] getRmatrix(){

        double epsilon = 1e-8;

        int dimen = MA1.numRows();

        Matrix I = new DenseMatrix(Matrices.identity(dimen));

        DenseMatrix mA1 = new DenseMatrix(dimen, dimen);
        DenseMatrix mA1I = new DenseMatrix(dimen, dimen);
        MA1.mult(-1, I, mA1);
        mA1.solve(I, mA1I);
        DenseMatrix H = new DenseMatrix(dimen, dimen);
        DenseMatrix L = new DenseMatrix(dimen, dimen);
        G = new DenseMatrix(dimen, dimen);
        DenseMatrix T = new DenseMatrix(dimen, dimen);
        U = new DenseMatrix(dimen, dimen);
        N = new DenseMatrix(dimen, dimen);
        DenseMatrix UA = new DenseMatrix(dimen, dimen);
        DenseMatrix M = new DenseMatrix(dimen, dimen);
        DenseMatrix TA = new DenseMatrix(dimen, dimen);
        R = new DenseMatrix(dimen, dimen);

        mA1I.mult(MA0, H);// H=(-A_1)^{-1}A_0

        mA1I.mult(MA2, L);// L=(-A_1)^{-1}A_2
        G.set(L);
        T.set(H);
        double compare = 1;

        while (compare > epsilon) {

            H.mult(L, U);// U=HL
            L.multAdd(H, U);// U=HL+LH
            H.mult(H, M);// M=(H)^2
            H = new DenseMatrix(dimen, dimen);
            mA1I = new DenseMatrix(dimen, dimen);
            // U.multAdd(-1, I, I, H);// H=(I-U)
            // TODO: check this!!
            H = (DenseMatrix) I.copy().add(-1, U);// H=(I-U)
            H.solve(I, mA1I);// H1I=(I-U)^{-1}
            H = new DenseMatrix(dimen, dimen);
            mA1I.mult(M, H);// H=(I-U)^{-1}*M
            M = new DenseMatrix(dimen, dimen);
            L.mult(L, M);// M=(L)^2
            L = new DenseMatrix(dimen, dimen);
            mA1I.mult(M, L);// L=(I-U)^{-1}*M
            T.multAdd(L, G);// G=G+TL
            TA.set(T);
            T = new DenseMatrix(dimen, dimen);
            TA.mult(H, T);// T=TH
            U = new DenseMatrix(dimen, dimen);

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

        U = new DenseMatrix(dimen, dimen);
        // MA0.multAdd(G, MA1, U);
        // TODO: check this!!!
        U = (DenseMatrix) MA0.multAdd(G, MA1.copy());
        mA1I = new DenseMatrix(dimen, dimen);
        U.mult(-1, I, UA);// U=A_1+A_0*G
        UA.solve(I, mA1I);// (I-U)^{-1}
        N.set(mA1I);// N=(I-U)^{-1}
        MA0.mult(mA1I, R);// R=A_0*H1I.

        return Matrices.getArray(R);
    }

    
	public String label() {
        return "MTJ Logarithmic reduction solver. This solver uses the MTJ p"
                + "ackage to handle matrices from a phase ditribution.";
    }
}

