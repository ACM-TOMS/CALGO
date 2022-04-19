/**
 * GeometrixSolver.java
 * Created: 2/09/2006
 */
package jmarkov.solvers;

import jmarkov.MarkovProcess;
import jmarkov.basic.State;
import jmarkov.basic.exceptions.NotUnichainException;

/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 */
public abstract class GeometrixSolver extends Solver {

    /**
     * Builds a Geometrix Solver with the given SimpleMarkovProcess.
     * @param mp The Markov Process for which the steady state
     *        probabilities are sought.
     */
    public GeometrixSolver(MarkovProcess mp) {
        super(mp);
    }

    /**
     * This process should be extended in order to compute the R
     * matrix of the QBD. The user can get information of the
     * SimpleMarkovProcess associated with this solver though the
     * methods
     * <code>getRates(), getGenerator, and getRate(State,State)</code>
     * @return a Matrix with the R matrix for the given QBD.
     * @throws NotUnichainException
     * @see jmarkov.SimpleMarkovProcess#getGenerator()
     * @see jmarkov.SimpleMarkovProcess#getRates()
     * @see jmarkov.SimpleMarkovProcess#getRate(State, State)
     */
    public abstract double[][] getRmatrix() throws NotUnichainException;
}
