/*
 * Created on 22/11/2004
 *
 */
package jmarkov.jmdp.solvers;

import java.util.ArrayList;
import java.util.List;

import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.Solution;
import jmarkov.basic.StateC;
import jmarkov.basic.States;
import jmarkov.basic.exceptions.StructureException;
import jmarkov.jmdp.StochasticShortestPath;

/**
 * This solver gives a solution for the minimization of the total cost criterion
 * for an infinite horizon MDP.
 * 
 * @author Juan F. Redondo, Andrés Sarmiento , German Riaño - Universidad de los
 *         Andes
 * @param <S>
 *            States class
 * @param <A>
 *            Action class.
 */
public class StochasticShortestPathSolver<S extends StateC, A extends Action>
        extends AbstractTotalSolver<S, A> {

    // this is a partner solver.
    private ValueIterationSolver<S, A> valueSolver;

    /**
     * Default contructor.
     * 
     * @param problem
     *            the structure of the problem of type StochasticShortestPath
     */
    public StochasticShortestPathSolver(StochasticShortestPath<S, A> problem) {
        super(problem);
        valueSolver = new ValueIterationSolver<S, A>(problem, 1);
    }

    /**
     * 
     * @see jmarkov.jmdp.solvers.Solver#getProblem()
     */
    @Override
    public StochasticShortestPath<S, A> getProblem() {
        return (StochasticShortestPath<S, A>) problem;
    }

    /**
     * This method calculates the expected value of valueFunction for the
     * current state and a specified action.
     * 
     * @param i
     *            Current State
     * @param a
     *            Action taken
     * @return Future value from this state.
     * @throws StructureException
     */
    public final double future(S i, A a) throws StructureException {
        double sum = 0;
        States<S> reachableStates = getProblem().reachable(i, a);
        for (S j : reachableStates) {
            sum += getProblem().modifiedProb(i, j, a) * valueFunction.get(j);
        }
        return sum;
    }

    /*
     * This is the method which activated the solver, a well-specified problem
     * must be created, then StochasticShortestPathSolver can solve the problem
     * 
     * @throws StructureException
     */
    @Override
    public Solution<S, A> solve() throws StructureException {
        valueSolver.init();
        long i = 0;
        double beforeDifference = Double.MAX_VALUE;
        double actualDifference = Double.MAX_VALUE;
        List<Double> l = new ArrayList<Double>();
        long initialTime = 0;
        int maxIteartions = 100 * getProblem().getAllStates().size();
        if (printProcessTime) {
            initialTime = System.currentTimeMillis();
        }

        while (actualDifference > valueSolver.getEpsilon()) {
            if (beforeDifference + valueSolver.getEpsilon() < actualDifference) {
                String msg = "\t Error: The problem solution won't be reached in a finite number of iterations";
                msg += "\t The optimal policy must be acyclic for finite termination";
                msg += "\t Possible solution: Verify the problem structure";
                throw new StructureException(msg);
            }
            if (i > maxIteartions) {
                StructureException ex = new StructureException(
                        "The value iteration method won't yield the optimal cost vector in a finite number of iterations");

                throw ex;

            }
            beforeDifference = actualDifference;
            if (valueSolver.usesErrorBounds())
                actualDifference = valueSolver.computeWithErrorBounds();
            else
                actualDifference = valueSolver.computeNoErrorBounds();
            l.add(new Double(actualDifference));
            i++;
        }
        if (printProcessTime) {
            valueSolver.processTime = System.currentTimeMillis() - initialTime;
        }
        return new Solution<S, A>(valueFunction, policy);
    }

    /**
     * Sets the best action to take in state i, in the variable bestAction. Note
     * that in this case StochasticShortestPathProblem Bertsekas expose a
     * transformation for the graph which modify the immediate reward function
     * and the transition probability, only to make a graph without
     * self-transition states. This will increase the finite termination
     * probability for the algorithm.
     * 
     * @param i
     *            state for which the best action is being determined
     * @return the new ValueFunction for this state.
     * @throws StructureException
     */

    protected double bestAction(S i) throws StructureException {
        double immediateRewardT;

        Actions<A> act = getProblem().feasibleActions(i);
        double val = 0;
        double maxSoFar = -Double.MAX_VALUE;
        for (A a : act) {
            if (i.isTerminal()) {
                immediateRewardT = getProblem().immediateCost(i, a);
            } else {
                immediateRewardT = getProblem().immediateCost(i, a)
                        + ((getProblem().immediateCost(i, a))
                                * getProblem().prob(i, i, a) / (1 - getProblem()
                                .prob(i, i, a)));
            }
            val = getProblem().operation(immediateRewardT, future(i, a));
            if (val > maxSoFar) {
                maxSoFar = val;
                valueSolver.bestAction = a;
            }
        }
        return maxSoFar;
    }

    /**
     * @see java.lang.Object#toString()
     */
    @Override
    public String description() {
        StringBuffer buf = new StringBuffer();
        buf
                .append("\n\t ***Stochastic Shortest Path Solver***\n"
                        + "____________________________________________________________");
       // buf.append("\nDiscount Factor = 1");
        if (valueSolver.usesGaussSeidel())
            buf.append("using Gauss-Seidel modification\n");
        if (valueSolver.usesErrorBounds())
            buf.append("using Error Bounds convergence\n");
        return buf.toString();
    }

    /**
     * @see Solver#label()
     */
    @Override
    public String label() {
        StringBuffer buf = new StringBuffer();
        buf.append("Stochastic Shortest Path Solver");
        // buf.append(", Discount Factor = 1");
        if (valueSolver.usesGaussSeidel())
            buf.append(", using Gauss-Seidel modification");
        if (valueSolver.usesErrorBounds())
            buf.append(", using Error Bounds convergence");
        buf.append(".");
        return buf.toString();
    }

    //
    // Overriden methods
    //

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
}
