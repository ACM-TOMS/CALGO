package jmarkov.jmdp;

import jmarkov.basic.Action;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;

/**
 * This class should ONLY be used in FINITE horizon deterministic problems. It
 * must be extended in order to represent the appropriate structure for each
 * FINITE Dynamic Programming problem. The user must implement at least the
 * functions that have been declared abstract. It´s also necessary to create one
 * of the extensions of the class Solver. By default, the program includes the
 * FiniteSolver class to solve finite horizon problems. PolicyIterationSolver
 * and ValueIterationSolver are only for infinite horizon problems. To solve the
 * problem follow the instructions in each of the solvers´ instructions.
 * 
 * @author Andres Sarmiento, Germán Riaño - Universidad de Los
 *         Andes
 * @param <S> States class
 * @param <A> Actions class
 * @see jmarkov.jmdp.solvers.FiniteSolver
 */

public abstract class FiniteDP<S extends State, A extends Action> extends
		FiniteMDP<S, A>

{

	//
	// constructors
	//

	/**
	 * Creates a new FINITE Dynamic Programming (DP) Problem.
	 * 
	 * @param states
	 *            set of all possible states.
	 * @param lastStage
	 *            number of the last stage.
	 */

	/*
	 * public DP(States states, int lastStage) { super(lastStage); }
	 */
	/**
	 * Creates a new FINITE Dynamic Programming (DP) Problem.
	 * @param initial initial set of known states.
	 * 
	 * @param lastStage
	 *            number of the last stage.
	 */

	public FiniteDP(States<S> initial, int lastStage) {
		super(initial, lastStage);
	}

	//
	// Abstract Methods
	//

	/**
	 * Final method. Cannot be implemented by any user.
	 */

	@Override
	public final double prob(S i, S j, A a, int t) {
		/*
		 * If a State j is passed to this function it is necessarily reachable.
		 */
		return 1.0;
	}

	/**
	 * Final function must not be extended by any user.
	 */
	@Override
	public final States<S> reachable(S i, A a, int t) {
		return new StatesSet<S>(destination(i, a, t));
	}

	/**
	 * State where the system will end up if action a is taken from state i at
	 * time t. The user must implement this method.
	 * @param i Current state
	 * @param a Current action
	 * @param t Time stage.
	 * @return Destination states
	 */
	public abstract S destination(S i, A a, int t);


}// class end

