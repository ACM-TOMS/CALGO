/**
 * LPSolver.java
 * Created: 7/04/2006
 */
package jmarkov.jmdp.solvers;

import jmarkov.basic.Action;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.exceptions.SolverException;

/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 * This interface defines some methods necesary to implement an LP solver. 
 * @param <S>
 * @param <A>
 */
public interface LPSolver<S extends State, A extends Action> {

    /**
     * Returns the time taken to build and write the MPS file.
     * @return Returns the buildTime.
     */
    public abstract long getBuildTime();

    /**
     * Return the time taken to solve the LP model.
     * @return Returns the lpSolveTime.
     */
    public abstract long getLpSolveTime();

    /**
     * Returns the time needed to build the Solution after the LP was
     * solved.
     * @return Returns the solBuildTime.
     */
    public abstract long getSolBuildTime();

    /**
     * The implementator classes should override this class to solve
     * the problem using the mpsFile that has been created.
     * @throws SolverException 
     */
    public abstract void solveLP() throws SolverException;

    /**
     * The implementator classes should override this class to build
     * the solution after the model has been solved.
     * @return The solution to the problem.
     * @throws SolverException 
     */
    public abstract Solution<S, A> buildSolution() throws SolverException;

}