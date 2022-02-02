package jmarkov.jmdp;

import java.io.PrintWriter;

import jmarkov.DebugReporter;
import jmarkov.basic.Action;
import jmarkov.basic.Policy;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.solvers.Solver;

/*
 * Created on 6/10/2004
 
 */

/**
 * This class is the main framework to build a Dynamic Programming Problem.
 * It was initially created to work over Markov Decision Problems which imply
 * random probabilities but can easily be worked out for deterministic problems
 * if the probabilities are set to one. This class should not me extended
 * directly on problems. The default package has FiniteMDP and InfiniteMDP
 * classes that are intended to be extended on problems. See the examples for a
 * clearer reference.
 * 
 * @author Andres Sarmiento, Germán Riaño - Universidad de Los Andes
 * @param <S>
 *            Current state.
 * @param <A>
 *            Action taken.
 * @see FiniteMDP
 * @see DTMDP
 */
public abstract class MDP<S extends State, A extends Action> {

    /**
     * States whether the problem is a finite horizon problem or not.
     */
    protected boolean finite = false;

    /**
     * Set of initial states.
     */
    protected States<S> initial;

    /**
     * Solver for the problem
     */

    private Solver<S, A> solver = null;

    /**
     * Reporter used for debug information.
     */
    protected DebugReporter reporter = new DebugReporter(1);

    /**
     * Returns the optimal policy. This function should be called only AFTER
     * solve has been called.
     * 
     * @return The optimal policy.
     * @throws SolverException
     *             Exception thrown if a solution cannot be found
     */
    public final Policy<S, A> getOptimalPolicy() throws SolverException {
        return getSolver().getOptimalPolicy();
    }

    /**
     * Indicates if the problems has been solved
     * 
     * @return true if solved
     */
    public final boolean isSolved() {
        return getSolver().isSolved();
    }

    /**
     * @return Returns true if the problem's horizon is finite.
     */
    public final boolean isFinite() {
        return finite;
    }

    /**
     * Solves the problem. If no solver has been defined, this used the default
     * solver.
     * 
     * @throws SolverException
     *             This Exception is thrown if a solution cannot be found.
     * 
     * @see #getDefaultSolver()
     */
    public final void solve() throws SolverException {
        if (solver == null) {
            solver = getDefaultSolver();
        }
        solver.solve();
    }

    /**
     * @return Returns the solver.
     */
    public Solver<S, A> getSolver() {
        if (solver == null)
            solver = getDefaultSolver();
        return solver;
    }

    /**
     * @param solver
     *            The solver to set.
     */
    public void setSolver(Solver<S, A> solver) {
        this.solver = solver;
        reporter.debug(3, "Solver set to " + solver);
    }

    /**
     * The class that extends MDP must define the default solver to use.
     * 
     * @return the solver to use for this problem.
     */
    protected abstract Solver<S, A> getDefaultSolver();

    /**
     * Returns the optimal ValueFunction. This causes the problem to be solved
     * if it has not been solved.
     * 
     * @return Returns the valueFunction.
     * @throws SolverException
     *             This exception is thrown if a solution cannot be found.
     */
    public ValueFunction<S> getOptimalValueFunction() throws SolverException {
        return getSolver().getOptimalValueFunction();
    }

    /**
     * The Operator between present and future costs. By default is sum, but can
     * be changed by the user, by overriding this method.
     * 
     * @param present
     *            Cost of current transition
     * @param future
     *            Cost of future transitions.
     * @return By default it returns present + future.
     */
    public double operation(double present, double future) {
        return (present + future);
    }

    /**
     * @return Returns the reporter.
     */
    public DebugReporter getReporter() {
        return reporter;
    }

    /**
     * @param reporter
     *            The reporter to set.
     */
    public void setReporter(DebugReporter reporter) {
        this.reporter = reporter;
    }

    /**
     * Prints a message in the reporter.
     * 
     * @param level
     *            maximum debug level at which to show message
     * @param message
     *            message
     * @see jmarkov.DebugReporter
     */
    public void debug(int level, String message) {
        reporter.debug(level, message);
    }

    /**
     * Prints debug information in the reporter.
     * 
     * @param level
     *            the level for the info
     * @param s
     *            Message
     * @param newline
     *            true if a new line is to be inserted
     * @param indent
     *            true if the info is indented according to level
     * @see jmarkov.DebugReporter
     * @see jmarkov.DebugReporter#debug(int, java.lang.String, boolean, boolean)
     */
    public void debug(int level, String s, boolean newline, boolean indent) {
        reporter.debug(level, s, newline, indent);
    }

    /**
     * Prints debug information in the reporter.
     * 
     * @param level
     *            the level for the info
     * @param s
     *            Message
     * @param newline
     *            true if a new line is to be inserted
     * @see jmarkov.DebugReporter
     * @see jmarkov.DebugReporter#debug(int, java.lang.String, boolean)
     */
    public void debug(int level, String s, boolean newline) {
        reporter.debug(level, s, newline);
    }

    /**
     * Gets the current debug level.
     * 
     * @return The current debug level
     * @see jmarkov.DebugReporter
     * @see jmarkov.DebugReporter#getDebugLevel()
     */
    public int getDebugLevel() {
        return reporter.getDebugLevel();
    }

    /**
     * Sets teh current level
     * 
     * @param level
     *            The new level to level
     * @see jmarkov.DebugReporter
     * @see jmarkov.DebugReporter#setDebugLevel(int)
     */
    public void setDebugLevel(int level) {
        reporter.setDebugLevel(level);
    }

    /**
     * Prints the solution to Standard output.
     */
    public void printSolution() {
        PrintWriter pw = new PrintWriter(System.out, true);
        printSolution(pw);
    }

    /**
     * Prints the solution to the given PrintWriter
     * 
     * @param pw
     *            The PrintWriter where the solution will be printed. It must
     *            have been initialized.
     */

    public void printSolution(PrintWriter pw) {
        getSolver().printSolution(pw);
    }
}
