package jmarkov.jmdp.solvers;

import java.util.Iterator;
import java.util.Map;

import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.DecisionRule;
import jmarkov.basic.Policy;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.NonStochasticException;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDP;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.sparse.BiCGstab;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;
import no.uib.cipr.matrix.sparse.IterativeSolver;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
import no.uib.cipr.matrix.sparse.SparseVector;

/*
 * Created on 18/09/2004
 */

/**
 * This class solves infinite horizon discounted problems using the
 * policy iteration algorithm. It extends Solver and should only be
 * used on infinite horizon problems. The objective function the
 * solver uses is the discounted cost. The result is a deterministic
 * optimal policy for the given structure. Policy Iteration is a
 * solver method that is always convergent in a finite number of
 * iterations. The algorithm has to solve a linear system of equations
 * as big as the amount of states. When there are too many states, it
 * is recommended to use other solvers, or using the modified policy
 * iteration (by using the second constructor). The advantage of using
 * Policy Iteration is that the result is the true optimal solution
 * and not an approximation as in other common methods. The method
 * starts with a policy. It solves the system of linear equations for
 * the value functions for that policy. With this values it looks for
 * a better policy. It then solves the value functions again and looks
 * for a better policy. If this policy is equal to the last policy
 * tried, it stops, in any other case it keeps improving the policy
 * and updating the value functions.
 * @author Andres Sarmiento, Germán Riaño, Daniel F. Silva
 * @param <S> States class.
 * @param <A> Actions class.
 */

public class PolicyIterationSolver<S extends State, A extends Action> extends
        AbstractDiscountedSolver<S, A> {

    // DenseVector vecValueFunction = null;

    /** Used to store local states */
    // protected List<S> localStates = new ArrayList<S>();
    private DenseVector costs;
    private boolean isOptimal = false;
    private boolean modifiedPolicy = false;
    private long initialIterations = 20;
    private double increasingFactor = 1.1;
    private long maxIterations;
    /** Used to store the number of iterations */
    protected long iterations;
    /** Used to store process time */
    protected long processTime = 0;

    private double epsilon = 0.0001;
    private boolean gaussSeidel = true;
    private boolean errorBounds = false;

    private DenseVector vecValueFunction = null;
    private FlexCompRowMatrix matrix = null;
    private DecisionRule<S, A> currentDecisionRule = null;

    /**
     * The constructor method exclusively receives a problem of the
     * type InfiniteMDP because this solver is only designed to work
     * on infinite horizon problems. This solver solves the discounted
     * objective function problem.
     * @param problem the structure of the problem of type InfiniteMDP
     * @param discountFactor represents how much less is the reward
     *        received in the next period instead of receiving it in
     *        the present period.
     */
    public PolicyIterationSolver(DTMDP<S, A> problem, double discountFactor) {
        this(problem, discountFactor, false);
    }

    /**
     * The constructor method exclusively receives a problem of the
     * type InfiniteMDP because this solver is only designed to work
     * on infinite horizon problems. This solver solves the discounted
     * objective function problem.
     * @param problem the structure of the problem of type InfiniteMDP
     * @param discountFactor represents how much less is the reward
     *        received in the next period instead of receiving it in
     *        the present period.
     * @param setModifiedPolicy
     */
    public PolicyIterationSolver(DTMDP<S, A> problem, double discountFactor,
            boolean setModifiedPolicy) {
        super(problem, discountFactor);
        modifiedPolicy = setModifiedPolicy;

        // SequentialBLAS theSeq = new SequentialBLAS();
        // matrix.scale(-discountFactor);
        // // matrix.addDiagonal(1.0);
        // int n = getProblem().getNumStates();
        // for (int i = 0; i < n; i++) {
        // matrix.add(i, i, 1.0);
        // }
    }

    /**
     * @return increasing factor of the maximum iterations.
     */
    public double getIncreasingFactor() {
        return increasingFactor;
    }

    /**
     * Sets the increasing factor of the maximum iterations of the
     * Modified policy iteration method. The first iterations are a
     * vague aproximation to the real value functions and need not be
     * exhaustive. But the last iterations must refine the value
     * functions in order to get better precision. The increasing
     * factor determines how many iteratinos are to be done in each
     * iteration. Faster growth will be more precise but
     * computationaly more expensive.
     * @param increasingFactor greater that 1. Determines max
     *        iterations growth.
     */
    public void setIncreasingFactor(double increasingFactor) {
        this.increasingFactor = increasingFactor;
    }

    /**
     * @return initial maximum iterations of the modified policy
     *         iteration algorithm.
     */
    public double getInitialIterations() {
        return initialIterations;
    }

    /**
     * Sets maximum iterations for the first run of the modified
     * policy iteration.
     * @param initialIterations
     */
    public void setInitialIterations(int initialIterations) {
        this.initialIterations = initialIterations;
    }

    @Override
    public Solution<S, A> solve() throws SolverException {
        long initialTime = System.currentTimeMillis();
        currentDecisionRule = initialDecisionRuleMyopic();
        // currentDecisionRule = initialDecisionRuleFirst();
        policy = new Policy<S, A>(currentDecisionRule);

        vecValueFunction = new DenseVector(getDiscreteProblem().getNumStates());
        matrix = buildMatrix(currentDecisionRule);

        iterations = 0;
        while (!isOptimal) {
            problem.debug(2, "Iteration " + iterations);
            getProblem().debug(3, "Current Rule = " + currentDecisionRule);
            valueFunction = policyEvaluation();
            getProblem().debug(3, "Current Value function = " + valueFunction);
            currentDecisionRule = policyImprovement();
            // policy.setDecisionRule(currentDecisionRule);
            iterations++;
            // assert (iterations < 20);
        }
        // step for precision
        if (modifiedPolicy) {
            matrix = buildMatrix(currentDecisionRule);
            valueFunction = solveMatrix();
            // valueFunction =
            // solveMatrixModified(currentDecisionRule);
        }
        // valueFunction = buildValueFunction(vecValueFunction);
        policy = new Policy<S, A>(currentDecisionRule);
        solved = true;
        processTime = System.currentTimeMillis() - initialTime;
        return new Solution<S, A>(valueFunction, policy);
    }

    /**
     * Builds the initial policy, setting the first available action
     * to each state.
     */
    private DecisionRule<S, A> initialDecisionRuleFirst() {
        valueFunction = new ValueFunction<S>();
        DecisionRule<S, A> localDecisionRule = new DecisionRule<S, A>();
        States<S> states = getProblem().getAllStates();
        for (S i : states) {
            Actions<A> availableActions = getProblem().feasibleActions(i);
            // build an initial policy
            localDecisionRule.set(i, availableActions.iterator().next());
            valueFunction.set(i, 0.0);
        }
        return localDecisionRule;
    }

    /**
     * Builds the initial policy, setting the action with lower
     * immediate cost.
     */
    private DecisionRule<S, A> initialDecisionRuleMyopic() {
        valueFunction = new ValueFunction<S>();
        DecisionRule<S, A> localDecisionRule = new DecisionRule<S, A>();
        States<S> states = getProblem().getAllStates();
        A bestAction = null;
        double bestVal = Double.MAX_VALUE;
        for (S state : states) {
            Actions<A> availableActions = getProblem().feasibleActions(state);
            // build an initial policy
            for (A action : availableActions) {
                double val = getDiscreteProblem().immediateCost(state, action);
                if (val < bestVal) {
                    bestAction = action;
                    bestVal = val;
                }
            }
            localDecisionRule.set(state, bestAction);
            valueFunction.set(state, bestVal);
        }
        return localDecisionRule;
    }

    private ValueFunction<S> policyEvaluation() throws SolverException {
        if (modifiedPolicy)
            valueFunction = solveMatrixModified(currentDecisionRule);
        else
            valueFunction = solveMatrix();
        return valueFunction;
    }

    private DecisionRule<S, A> policyImprovement() throws SolverException {

        States<S> sts = getProblem().getAllStates();
        DecisionRule<S, A> newDecisionRule;
        newDecisionRule = new DecisionRule<S, A>(currentDecisionRule);
        double val;
        Iterator<Map.Entry<S, A>> itCurDR = currentDecisionRule.iterator();
        Iterator<Map.Entry<S, A>> itNewDR = newDecisionRule.iterator();
        for (S i : sts) {
            Actions<A> actions = getProblem().feasibleActions(i);
            A bestAction = null;
            double bestValue = Double.MAX_VALUE;
            for (A a : actions) {
                val = getProblem().operation(
                        getDiscreteProblem().immediateCost(i, a),
                        future(i, a, discountFactor));
                if (val < bestValue) {// minimization
                    bestValue = val;
                    bestAction = a;
                }
            }
            if (!modifiedPolicy) {
                Map.Entry<S, A> curDRentry = itCurDR.next();
                A curAction = curDRentry.getValue();
                if (!bestAction.equals(curAction)) {
                    matrix.setRow(i.getIndex(), buildRowVector(i, bestAction));
                    costs.set(i.getIndex(), getDiscreteProblem().immediateCost(
                            i, bestAction));
                }
            }
            Map.Entry<S, A> newDRentry = itNewDR.next();
            newDRentry.setValue(bestAction);
        }
        isOptimal = (currentDecisionRule.equals(newDecisionRule));
        return newDecisionRule;
    }

    /**
     * Builds the i-th row vector of the matrix (I-beta P)
     * corresponding to the given action.
     * @param i the state
     * @param a the action
     * @return A SparseVector.
     */
    private SparseVector buildRowVector(S i, A a) {
        int n = getDiscreteProblem().getNumStates();
        States<S> reachableStates = getDiscreteProblem().reachable(i, a);
        SparseVector vec = new SparseVector(n, reachableStates.size());
        double sum = 0.0;
        for (S j : reachableStates) {
            double probability = getDiscreteProblem().prob(i, j, a);
            sum += probability;
            assert (probability >= 0);
            // make sure is the same.
            j = getDiscreteProblem().getAllStates().get(j);
            if (probability > 0) {
                vec.set(j.getIndex(), probability);
            }

        }
        vec.scale(-discountFactor);
        vec.add(i.getIndex(), 1.0);
        if (Math.abs(sum - 1.0) > 1e-5) {
            throw new NonStochasticException(
                    "Probabilities do not add up to 1 for state " + i
                            + ", and action " + a + ", sum = " + sum);
        }
        return vec;
    }

    /**
     * This method builds the Probability Transision matrix for a
     * specified policy. The solver then transforms this matrix and
     * uses it to solve the value functions for each state.
     * @param currentDecisionRule the policy under which the
     *        probability matrix is to be built.
     * @return the probability matrix.
     */

    private FlexCompRowMatrix buildMatrix(DecisionRule<S, A> currentDecisionRule) {
        // DenseRowMatrix M = new DenseRowMatrix(localStates.size(),
        // localStates
        // .size());
        StatesSet<S> stts = getDiscreteProblem().getAllStates();
        int n = stts.size();
        // int[] nonZeros = new int[n];
        // for (S i : stts) {
        // A a = localPolicy.getAction(i);
        // States<S> reached = getDiscreteProblem().reachable(i, a);
        // nonZeros[i.getIndex()] = reached.size();
        // }
        // Matrix M = new CompRowMatrix(n, n, nonZeros);
        FlexCompRowMatrix matrix = new FlexCompRowMatrix(n, n);
        costs = new DenseVector(n);
        for (S i : stts) {
            A a = currentDecisionRule.getAction(i);
            matrix.setRow(i.getIndex(), buildRowVector(i, a));
            costs.set(i.getIndex(), getDiscreteProblem().immediateCost(i, a));
        }
        return matrix;
    }

    /**
     * This method is used by the PolicyIterationSolver to solve the
     * linear system of equations to determine the value functions of
     * each state for a given policy.
     * @return a DenseVector (type defined in the JMP package
     *         documentation) with the value functions for each state.
     *         The index for each state are the same ones determined
     *         in the localStates ArrayList
     * @throws SolverException
     */
    protected ValueFunction<S> solveMatrix() throws SolverException {
        getProblem().debug(4, "Matrix to solve:\n" + matrix);
        try {
            IterativeSolver solver = new BiCGstab(vecValueFunction);
            solver.solve(matrix, costs, vecValueFunction);
        } catch (IterativeSolverNotConvergedException e) {
            throw new SolverException(
                    "Policy iteration Solver: error solving linear system.", e);
        }
        return buildValueFunction(vecValueFunction);
    }

    /**
     * This method is used by the PolicyIterationSolver to solve the
     * linear system of equations to determine the value functions of
     * each state for a given policy.
     * @return a DenseVector (type defined in the JMP package
     *         documentation) with the value functions for each state.
     *         The index for each state are the same ones determined
     *         in the localStates ArrayList declared as static.
     */
    protected ValueFunction<S> solveMatrixModified(
            DecisionRule<S, A> localDecisionRule) {
        States<S> st = getProblem().getAllStates();
        ValueFunction<S> vf = new ValueFunction<S>(valueFunction);
        ValueFunction<S> vf2 = new ValueFunction<S>(valueFunction);
        double maxDifference = 0;
        int localIterations = 0;
        boolean toContinue = true;

        while (toContinue) {
            Iterator<Map.Entry<S, Double>> it = vf.iterator();
            Iterator<Map.Entry<S, Double>> it2 = vf2.iterator();
            maxDifference = 0;
            double bound = (1 - discountFactor) * epsilon
                    / (2 * discountFactor);
            for (S i : st) {
                A a = localDecisionRule.getAction(i);
                double val = getProblem().operation(
                        getDiscreteProblem().immediateCost(i, a),
                        future(i, a, discountFactor, vf));
                Map.Entry<S, Double> entry = it.next();
                Map.Entry<S, Double> entry2 = it2.next();
                double diff = Math.abs(val - entry.getValue());
                if (maxDifference < diff)
                    maxDifference = diff;
                if (gaussSeidel)
                    entry.setValue(val);
                else
                    entry2.setValue(val);
            }
            localIterations++;
            if (!gaussSeidel)
                vf = vf2;
            toContinue = (localIterations < initialIterations);
            if (isOptimal)
                toContinue = (maxDifference > bound);
        }
        initialIterations = (long) Math.ceil(increasingFactor
                * initialIterations);
        // isOptimal = (maxDifference <=
        // epsilon*(1-discountFactor)/(2*discountFactor));
        return vf;
    }

    /**
     * Stores the vector information in the valueFunction object
     */
    private ValueFunction<S> buildValueFunction(DenseVector vec) {
        ValueFunction<S> vf = new ValueFunction<S>();
        States<S> stts = getDiscreteProblem().getAllStates();
        int i = 0;
        for (S s : stts) {
            vf.set(s, vec.get(i));
            i++;
        }
        return vf;
    }

    // /**
    // * @see jmdp.solvers.Solver#getProcessTime()
    // */
    // @Override
    // public long getProcessTime() {
    // return processTime;
    // }
    //
    //	
    // /**
    // * @return Returns the iterations.
    // */
    // @Override
    // public final int getIterations() {
    // return iterations;
    // }

    /**
     * Activates the modified policy iteration algorithm.
     * @param val True if the modified policy iteration is to be used.
     */
    public void setModifiedPolicy(boolean val) {
        modifiedPolicy = val;
    }

    @Override
    public String description() {
        return (modifiedPolicy) ? "Modified " : ""
                + "Policy Iteration Solver\n" + "Discount Factor = "
                + discountFactor;
    }

    @Override
    public String label() {
        return "Policy Iter. Solver(disc)";
    }

    //
    // Overriden methods
    //

    /**
     * @return Returns the processTime.
     */
    @Override
    public final long getProcessTime() {
        return processTime;
    }

    /**
     * @return Returns the iterations.
     */
    @Override
    public final long getIterations() {
        return iterations;
    }

}// class end
