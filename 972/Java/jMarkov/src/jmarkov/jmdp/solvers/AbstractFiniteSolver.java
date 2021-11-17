/*
 * AbstractFiniteSolver.java
 * Created: Jul 15, 2005
 */
package jmarkov.jmdp.solvers;

import jmarkov.basic.Action;
import jmarkov.basic.State;
import jmarkov.jmdp.FiniteMDP;

/**
 * Structural class for solvers to extend in order to solve finite horizon problems.
 * 
 * @author Andrés Sarmiento, Germán Riaño. Universidad de los Andes. (C) 2005
 * @param <S> State class
 * @param <A> Action class
 */
public abstract class AbstractFiniteSolver<S extends State, A extends Action>
		extends Solver<S, A> {


	
	/**
	 * @param problem finite horizon problem to be solved
	 */
	protected AbstractFiniteSolver(FiniteMDP<S,A> problem) {
		super(problem);
	}
	
	/**
	 * Returns the problem associated wit this solver.
	 * 
	 * @return the problem associated wit this solver.
	 */
	@Override
	public FiniteMDP<S, A> getProblem() {
		return (FiniteMDP<S, A>)problem;
	}

}
