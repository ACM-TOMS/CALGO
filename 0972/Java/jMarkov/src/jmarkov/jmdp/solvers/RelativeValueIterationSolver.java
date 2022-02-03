package jmarkov.jmdp.solvers;

import java.io.PrintWriter;

import jmarkov.basic.Action;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.jmdp.CTMDP;
import jmarkov.jmdp.DTMDP;

/**
 * This class solves the average cost criteria for infinite horizon problems
 * 
 * @author Andres Sarmiento, Germán Riaño and Daniel F. Silva
 * 
 * @param <S>
 *            States class.
 * @param <A>
 *            Actions class.
 */
public class RelativeValueIterationSolver<S extends State, A extends Action>
        extends AbstractAverageSolver<S, A> {

    private ValueIterationSolver<S, A> valueSolver;

    //
    // Constructors
    //

    /**
     * The constructor method exclusively receives a discrte time infinite
     * horizon problem of the type DTMDP.
     * 
     * @param problem
     *            the structure of the problem of type InfiniteMDP
     */

    public RelativeValueIterationSolver(DTMDP<S, A> problem) {
        super(problem);
        valueSolver = new ValueIterationSolver<S, A>(problem, false);
        // valueSolver.isAverage = true;
    }

    /**
     * Creates a new solver for the given discrete time, infinite horizon
     * problem. It uses the modified relative value iteration method. For
     * details, consult the User's Manual. The factor helps avoiding periodicity
     * in the chain.
     * 
     * @param problem
     *            problem
     * @param factor
     *            factor
     */
    public RelativeValueIterationSolver(DTMDP<S, A> problem, double factor) {
        super(problem);
        valueSolver = new ValueIterationSolver<S, A>(problem, true);
        // valueSolver.isAverage = true;
        // valueSolver.modifiedAverage = true;
        setFactor(factor);
    }

    /**
     * Creates a new solver for a continuous time, infinite horizon problem.
     * 
     * @param problem
     *            continuous time, infinite horizon problem
     */
    public RelativeValueIterationSolver(CTMDP<S, A> problem) {
        super(problem);
        valueSolver = new ValueIterationSolver<S, A>(problem, false);
        // valueSolver.isAverage = true;
    }

    /**
     * Creates a new solver for a continuous time, infinite horizon problem to
     * be solved with the modified relative value iteration method.
     * 
     * @param problem
     *            continuous time, infinite horizon problem
     * @param factor
     */
    public RelativeValueIterationSolver(CTMDP<S, A> problem, double factor) {
        super(problem);
        valueSolver = new ValueIterationSolver<S, A>(problem, true);
        // valueSolver.isAverage = true;
        // valueSolver.modifiedAverage = true;
        setFactor(factor);
    }

    //
    // Methods
    //

    @Override
    public void setPrintValueFunction(boolean val) {
        this.printValueFunction = val;
        valueSolver.setPrintValueFunction(val);
    }

    /**
     * Sets the factor for the modified relative value iteration method.
     * 
     * @param factor
     *            A number between 0 and 1.
     */
    public void setFactor(double factor) {

        if (factor < 1 && factor > 0) {
            valueSolver.discountFactor = factor;
            // valueSolver.modifiedAverage = true;
        } else if (!valueSolver.usesModifiedAverage()) {
            throw new IllegalArgumentException(
                    "Trying to set convergence factor when Mofified Average is not enabled!");
        } else

            throw new IllegalArgumentException(
                    "Factor set outside the valid range (0,1).");
    }

    //
    // Overriden methods
    //

    @Override
    public String label() {
        return "Relative Value Iteration Solver";
    }

    @Override
    public Solution<S, A> solve() {
        Solution<S, A> sol = valueSolver.solve();
        solved = true;
        valueFunction = sol.getValueFunction();
        policy = sol.getPolicy();
        return sol;
    }

    /**
     * @return Returns the processTime.
     */
    @Override
    public final long getProcessTime() {
        return valueSolver.getProcessTime();
    }

    /**
     * @return Returns the iterations.
     */
    @Override
    public final long getIterations() {
        return valueSolver.getIterations();
    }

    public void setPrintBias(boolean val) {
        valueSolver.setPrintBias(val);
    }
    
    public void setPrintGain(boolean val) {
        valueSolver.setPrintGain(val);
    }

    public double getGain() {
        return valueSolver.getGain();
    }

    
    @Override
    public void printSolution(PrintWriter pw) {
    	valueSolver.printSolution(pw);
    }

    /**
     * Prints the solution in the default PrintWriter (System.out)
     * 
     * @throws Exception
     */
    @Override
    public void printSolution() throws Exception {
        printSolution(new PrintWriter(System.out));
    }
    
}
