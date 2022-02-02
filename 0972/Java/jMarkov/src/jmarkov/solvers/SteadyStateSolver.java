/*
 * Created on 13/06/2005
 */
package jmarkov.solvers;

import jmarkov.MarkovProcess;
import jmarkov.basic.State;
import jmarkov.basic.exceptions.NotUnichainException;

/**
 * An abstract class for steady state solver. Any solver for steady state should
 * extend this class.
 * 
 * @author German Riano. universidad de los Andes.
 */
public abstract class SteadyStateSolver extends Solver {

	/**
	 * Builds a Steady State Solver with the given SimpleMarkovProcess.
	 * 
	 * @param mp
	 *            The Markov Process for which the steady state probabilities
	 *            are sought.
	 */
	public SteadyStateSolver(MarkovProcess mp) {
		super(mp);
	}

	/**
	 * This process should be extended in order to compute the steady State
	 * probabilities of the MarkovChain. The user can get information of the
	 * SimpleMarkovProcess associated with this solver though the methods
	 * <code>getRates(), getGenerator, and getRate(State,State)</code>
	 * 
	 * @return an array with the Steady state probabilities for the given
	 *         problem.
	 * @throws NotUnichainException 
	 * @see jmarkov.SimpleMarkovProcess#getGenerator()
	 * @see jmarkov.SimpleMarkovProcess#getRates()
	 * @see jmarkov.SimpleMarkovProcess#getRate(State, State)
	 */
	public abstract double[] getSteadyState() throws NotUnichainException;

}
