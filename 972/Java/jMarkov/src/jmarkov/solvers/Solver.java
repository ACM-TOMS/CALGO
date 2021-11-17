/*
 * Created on 13/06/2005
 */
package jmarkov.solvers;

import jmarkov.MarkovProcess;
import jmarkov.basic.JMarkovElement;

/**
 * This abstract class has to be extended in order to 
 * implement solvers for Steady State and Transient probabilities.
 * Most users do not need to implement this class since a default
 * solver is provided.
 * 
 * @author German Riano. Universidad de los Andes.
 * @see jmarkov.solvers.SteadyStateSolver
 * @see jmarkov.solvers.TransientSolver
 * @see jmarkov.solvers.JamaSolver
 */


public abstract class Solver implements JMarkovElement {
/** The Markovprocess being solved */
  protected MarkovProcess<?,?> mp; 
  
  /**
   * Build a solver for the given SimpleMarkovProcess
 * @param mp Markov Process to be solved.
   */
  public Solver(MarkovProcess mp){
    this.mp = mp;
  }
  
  /**
   * Returns the Markov process currently being solved by this solver.
   * @return the current Markov Process associated with this solver.
   */
  public final MarkovProcess getMP(){
    return mp;
  }
  
  /**
   * The name of this solver. This should be implemented by the extending classes.. 
   */
  public abstract  String label();

  /**
   * Return the name of the Solver.
   * @see #label()
   */
  @Override
	public final String toString(){
	  return label();
  }

}
