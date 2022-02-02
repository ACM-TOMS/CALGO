/*
 * Created on 13/06/2005
 */
package jmarkov.solvers;

import jmarkov.MarkovProcess;
import jmarkov.basic.exceptions.NotUnichainException;
import Jama.Matrix;

/**
 * Solver implementation for steady state, using JAMA
 * 
 * @author Germán Riaño
 */
public final class JamaSolver extends SteadyStateSolver {

    /**
     * @param mp
     */
    public JamaSolver(MarkovProcess mp) {
        super(mp);
    }

    /**
     * It find the steady state probabilities. If no Solution is found an array
     * of 0's is returned
     */
    @Override
    public double[] getSteadyState() throws NotUnichainException {
        Matrix Q = (new Matrix(mp.getGenerator())).copy();
        mp.debug(1, "Computing  Steady State Probabilities ..");
        int n = Q.getRowDimension();
        // Replace last column for 1's
        for (int i = 0; i < n; i++) {
            Q.set(i, n - 1, 1.0);
        }
        // New row vector with 1 at last component
        Matrix b = new Matrix(1, n);
        b.set(0, n - 1, 1);
        Matrix pi;
        try {
            long initialTime = System.currentTimeMillis();
            pi = Q.solveTranspose(b);
            long processTime = System.currentTimeMillis() - initialTime;
            mp.debug(1, "Solve time : " + processTime + " milliseconds.");
        } catch (Exception e) {
            throw new NotUnichainException(
                    "Exception Solving Steady state probabilities "
                            + "with JAMA solver. ", e);
        }
        return (pi.getRowPackedCopy());
    }

    /**
     * @see jmarkov.solvers.Solver#label()
     */
    @Override
    public String label() {
        return "JAMA solver";
    }

    public String description() {
        return "JAMA SOLVER. This solver uses the JAMA package to handle matrices.";
    }

}
