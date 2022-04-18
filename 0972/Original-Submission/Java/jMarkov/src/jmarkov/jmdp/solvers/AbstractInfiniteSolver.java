/*
 * AbstractInfiniteSolver.java
 * Created: Jul 15, 2005
 */
package jmarkov.jmdp.solvers;

import java.io.PrintWriter;

import jmarkov.basic.Action;
import jmarkov.basic.State;
import jmarkov.jmdp.CT2DTConverter;
import jmarkov.jmdp.CTMDP;
import jmarkov.jmdp.DTMDP;
import jmarkov.jmdp.InfiniteMDP;

/**
 * Structural class to be extended by solvers in order to solve infinite horizon
 * problems
 * 
 * @author Germán Riaño. Universidad de los Andes. (C) 2005
 * @param <S>
 *            States class
 * @param <A>
 *            Actions class
 */
public abstract class AbstractInfiniteSolver<S extends State, A extends Action>
        extends Solver<S, A> {

    /**
     * Constructor method for Discrete Time Markov Decision Processes to be
     * solved for discounted cost.
     * 
     * @param problem
     *            Discrete Time Markov Decision Process of type DTMDP
     */
    protected AbstractInfiniteSolver(DTMDP<S, A> problem) {
        super(problem);
    }

    /**
     * Creates a solver for an infinite horizon continuous time problem
     * 
     * @param problem
     *            continuous time problem
     */
    protected AbstractInfiniteSolver(CTMDP<S, A> problem) {
        super(new CT2DTConverter<S, A>(problem));
    }

    /**
     * Returns the problem associated with this solver.
     * 
     * @return the problem associated with this solver.
     */
    @Override
    public InfiniteMDP<S, A> getProblem() {
        return (InfiniteMDP<S, A>) problem;

    }

    /**
     * 
     * @return discrete time problem
     */
    protected DTMDP<S, A> getDiscreteProblem() {
        return (DTMDP<S, A>) getProblem();
    }

    /**
     * @return Returns the iterations in the last solve.
     */
    public abstract long getIterations();

    @Override
    public void printSolution(PrintWriter pw) {
        super.printSolution(pw);
        pw.println(getIterations() + " iterations");

    }

    /**
     * @return Returns the gain in the optimal policy.
     */

}
