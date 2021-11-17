/*
 * Created on 13/06/2005
 */
package jmarkov.solvers;

import jmarkov.MarkovProcess;
import jmarkov.basic.State;

/**
 * An abstract class for Transient solvers. At least the method
 * getTransientProbs(double time, State i0) has to be implemented. All others
 * call this method, but the user can provide more efficient implementations.
 * 
 * @see #getTransientProbs(double, State)
 * 
 * @author Germán Riaño. Universidad de los Andes.
 * 
 */
public abstract class TransientSolver extends Solver {

    /**
     * Build a solver with the asssocieted Markov Process.
     * 
     * @param mp
     */
    public TransientSolver(MarkovProcess mp) {
        super(mp);
    }

    /**
     * Computes the steady state probabilities at this given time, assuming the
     * Markov Chain starts in the given state i0.
     * 
     * @param time
     * @param i0
     *            Initial State.
     * @return probabilities array
     */

    public abstract double[] getTransientProbs(double time, State i0);

    /**
     * Computes the steady state probabilities at this given times, assuming the
     * Markov Chain starts in the given state i0.
     * 
     * @param times
     *            An array with the times at which the probabilities are to be
     *            evaluated.
     * @param i0
     *            The initial state (at time t=0).
     * @return probabilities array for each state. The (i,j) entry on the
     *         returned state represents the steady state probability for state
     *         i at time times[j].
     */
    public double[][] getTransientProbs(double times[], State i0) {
        double results[][] = new double[times.length][];
        for (int i = 0; i < times.length; i++) {
            results[i] = getTransientProbs(times[i], i0);
        }
        return results;
    }

    /**
     * Computes the steady state probabilities at times delta, 2delta,
     * 3delta,..., assuming the Markov Chain starts in the given state i0.
     * 
     * @param NumberPoints
     * @param delta
     *            the time gap between measurements.
     * @param i0
     *            Initial state.
     * @return probabilities array for each state. The (i,j) entry on the
     *         returned state represents the steady state probability for state
     *         i at time j * delta.
     */
    public double[][] getTransientProbs(int NumberPoints, double delta, State i0) {
        double times[] = new double[NumberPoints];
        for (int i = 0; i < NumberPoints; i++) {
            times[i] = delta * i;
        }
        return getTransientProbs(times, i0);
    }

}
