package jmarkov.jmdp.solvers;

import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.DecisionRule;
import jmarkov.basic.Policy;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.FiniteMDP;

/*
 * Created on 14-sep-2004
 */

/**
 * This class belongs to the set of default solvers included in the
 * jmdp package. It extends Solver and should only be used on FINITE
 * horizon problems. The objective function of the solver is to
 * minimize the the total cost. The result is a deterministic optimal
 * policy for the given structure.
 * @author Germán Riaño and Andres Sarmiento - Universidad de Los
 *         Andes
 * @param <S> States class.
 * @param <A> Actions class.
 */
public class FiniteSolver<S extends State, A extends Action> extends
        AbstractFiniteSolver<S, A> {

    private long processTime = 0;

    /**
     * Initialized the solver with the given getProblem().
     * @param problem The problem to be solved.
     */
    public FiniteSolver(FiniteMDP<S, A> problem) {
        super(problem);
        policy = new Policy<S, A>(getProblem().getHorizon());
    }

    @Override
    public final Solution<S, A> solve() {
        init();
        for (int t = getProblem().getHorizon() - 1; t >= 0; t--) {
            compute(t);
        }
        solved = true;
        return new Solution<S, A>(valueFunction, policy);
    }

    /**
     * This method calculates the expected value of valueFunction for
     * the current state i and a specified action a at the given stage
     * t.
     * @param i Current State
     * @param a Action taken
     * @param t Time stage
     * @return The value.
     * @throws NullPointerException
     */
    protected final double future(S i, A a, int t) throws NullPointerException {
        double sum = 0;
        States<S> reachableStates = getProblem().reachable(i, a, t);
        try {
            for (S j : reachableStates) {
                sum += getProblem().prob(i, j, a, t)
                        * getValueFunction().get(j);
            }
        } catch (NullPointerException e) {
            throw e;
        }
        return sum;
    }

    private void init() {
        States<S> st = getProblem().getStates(getProblem().getHorizon());
        policy = new Policy<S, A>(getProblem().getHorizon());
        double cost;
        for (S i : st) {
            cost = getProblem().finalCost(i);
            valueFunction.set(i, cost);
            getProblem().debug(4,
                    "Final cost for state " + i + " set at " + cost);
        }
    }

    private void compute(int t) {

        States<S> st = getProblem().getStates(t);
        DecisionRule<S, A> decisionRuleCompute = new DecisionRule<S, A>();
        ValueFunction<S> vF = new ValueFunction<S>();
        double val = 0;

        for (S i : st) {
            Actions<A> act = getProblem().feasibleActions(i, t);
            getProblem().debug(4, "Actions found for state " + i +":" + act);
            A Best_a = null;
            double minSoFar = Double.MAX_VALUE;
            for (A a : act) {
                try {
                    val = getProblem().operation(
                            getProblem().immediateCost(i, a, t),
                            future(i, a, t));
                } catch (NullPointerException e) {
                    continue;
                }
                if (val < minSoFar) {// minimization
                    minSoFar = val;
                    Best_a = a;
                }
            }
            vF.set(i, minSoFar);
            decisionRuleCompute.set(i, Best_a);
        }
        policy.setDecisionRule(decisionRuleCompute, t);
        valueFunction = vF;
    }

    /**
     * Prints out the policy
     * @param initial
     * @return a string with the optimal policy
     * @throws SolverException
     */
    public String bestPolicy(S initial) throws SolverException {
        StringBuffer buf = new StringBuffer(40);
        buf.append("*** Best Policy (starting in " + initial + " )***\n");
        buf.append(print(initial, 0));
        return buf.toString();
    }

    private StringBuffer print(S i, int t) throws SolverException {
        StringBuffer str = new StringBuffer(100);
        Policy<S, A> pol = getProblem().getOptimalPolicy();
        if (t < getProblem().getHorizon()) {
            // A a = getProblem().policy.getAction(i, t);
            A a = pol.getAction(i, t);
            str.append("Stage: " + t);
            str.append(" State: " + i);
            str.append(" Take Action: " + a + "\n");

            States<S> sttes = getProblem().reachable(i, a, t);
            for (S j : sttes) {
                str.append(print(j, t + 1));
            }
        }
        return str;
    }

    /**
     * @return Returns the processTime.
     */
    @Override
    public final long getProcessTime() {
        return processTime;
    }

    /**
     * @see jmarkov.jmdp.solvers.Solver#toString()
     */
    @Override
    public String label() {
        return "Default Finite horizon solver";
    }

}// class end
