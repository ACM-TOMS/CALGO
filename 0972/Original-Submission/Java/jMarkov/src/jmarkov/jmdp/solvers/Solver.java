package jmarkov.jmdp.solvers;

import java.io.PrintWriter;

import jmarkov.basic.Action;
import jmarkov.basic.JMarkovElement;
import jmarkov.basic.Policy;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.MDP;

/*
 * Created on 15/09/2004
 */

/**
 * Structural class for every solver. Any solver that a user implements must
 * extend this class.
 * 
 * @author Andres Sarmiento, Germán Riaño - Universidad de Los Andes
 * @param <S>
 *            States class.
 * @param <A>
 *            Actions class.
 */
public abstract class Solver<S extends State, A extends Action> implements
        JMarkovElement {

    /** The policy Object. This should be written by the solver. */
    protected Policy<S, A> policy = null;
    /** The value function, to be written by the solver */
    protected ValueFunction<S> valueFunction = new ValueFunction<S>();
    /** Marker to indicate that the problem has been solved */
    protected boolean solved = false;
    /** The problem to be solved */
    protected MDP<S, A> problem = null;

    /** True if the process time is to be reported */
    protected boolean printProcessTime = false;
    /** true if the value function is to be reported */
    protected boolean printValueFunction = false;

    /**
     * Default constructor. Receives the problem to solve. All sub classes MUST
     * call this constructor.
     * 
     * @param problem
     *            to be solved.
     */
    protected Solver(MDP<S, A> problem) {
        this.problem = problem;
        problem.setSolver(this);
    }

    /**
     * Returns the problem associated wit this solver.
     * 
     * @return the problem associated with this solver.
     */
    public MDP<S, A> getProblem() {
        return problem;
    }

    /**
     * Called to solve the problem. This method MUST write the local variable
     * policy and valueFunction.
     * 
     * @return The solution Object taht contains the plicy and value fuenction.
     * @throws SolverException
     *             This exception is thrown if the solver cannot find a solution
     *             for some reason.
     */
    public abstract Solution<S, A> solve() throws SolverException;

    /**
     * Gets the optimal policy. It solves the problem if it has not been solved.
     * 
     * @return the optimal Policy.
     * @throws SolverException
     * @see Policy
     */
    public final Policy<S, A> getOptimalPolicy() throws SolverException {
        if (!solved) {
            Solution<S, A> sol = solve();
            policy = sol.getPolicy();
        }
        return policy;
    }

    /**
     * If the problem is solved, it will return the optimal value function.
     * Otherwise it returns the current valueFunction
     * 
     * @return the value function in the solver.
     */
    public final ValueFunction<S> getValueFunction() {
        return valueFunction;
    }

    /**
     * Gets the optimal ValueFunction.
     * 
     * @return the optimal ValueFunction.
     * @throws SolverException
     * @see ValueFunction
     */
    public final ValueFunction<S> getOptimalValueFunction()
            throws SolverException {
        if (!solved) {
            Solution<S, A> sol = solve();
            valueFunction = sol.getValueFunction();
        }
        return valueFunction;
    }

    /**
     * Tells whether the problem has been solved.
     * 
     * @return true if the problem has been solved
     */
    public final boolean isSolved() {
        return solved;
    }

    /**
     * The sub classes must return the Solver name.
     * 
     * @see #toString()
     */
    public abstract String label();

    public String description() {
        return "Solver: " + label() + "\nClass: " + this.getClass().getName();
    }

    /**
     * This calls label().
     * 
     * @see #toString()
     */
    @Override
    public final String toString() {
        return label();
    }

    /**
     * @return Returns the processTime of the last solve. Use
     *         <code>System.currentTimeMillis()</code> to get the current
     *         time.
     */
    public abstract long getProcessTime();

    /**
     * Option to print the time spent solving the problem. It is set to false by
     * default.
     * 
     * @param val
     *            True if the Process tiem is to be reported, false otherwise.
     */
    public void setPrintProcessTime(boolean val) {
        this.printProcessTime = val;
    }

    /**
     * Option to print the final value function for each state. It is set to
     * false by default.
     * 
     * @param val
     *            True if the value function is to be reported.
     */
    public void setPrintValueFunction(boolean val) {
        this.printValueFunction = val;
    }

    /**
     * Prints the solution on a given PrintWriter.
     * 
     * @param pw
     * @see java.io.PrintWriter
     */
    public void printSolution(PrintWriter pw) {
        pw.println(this);
        try {
            getOptimalPolicy().print(pw);

            if (printValueFunction)
                valueFunction.print(pw);
            if (printProcessTime) {
                pw.println("Process time = " + getProcessTime()
                        + " milliseconds");
            }
        } catch (SolverException e) {
            pw.print(" Error solving the problem :" +e);
        }
    }

    /**
     * Prints the solution in the default PrintWriter (System.out)
     * 
     * @throws Exception
     */
    public void printSolution() throws Exception {
        printSolution(new PrintWriter(System.out));
    }

}
