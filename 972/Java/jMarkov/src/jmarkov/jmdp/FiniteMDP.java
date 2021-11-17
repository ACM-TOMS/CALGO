package jmarkov.jmdp;

import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.jmdp.solvers.FiniteSolver;
import jmarkov.jmdp.solvers.Solver;

/**
 * This class should ONLY be used in FINITE horizon problems. It must
 * be extended in order to represent the appropriate structure for
 * each FINITE Dynamic Programming problem. The user must implement at
 * least the functions that have been declared abstract. It´s also
 * necessary to create one of the extensions of the class Solver. By
 * default, the program includes the FiniteSolver class to solve
 * finite horizon problems. PolicyIterationSolver and
 * ValueIterationSolver are only for infinite horizon problems. To
 * solve the problem follow the instructions in each of the solvers´
 * instructions.
 * @author Andres Sarmiento, Germán Riaño - Universidad de Los Andes
 * @param <S> States class
 * @param <A> Actions class
 * @see jmarkov.jmdp.solvers.FiniteSolver
 */

public abstract class FiniteMDP<S extends State, A extends Action> extends
        MDP<S, A> {

    /**
     * Time horizon. The last stage in the problem.
     */
    protected int horizon;
    private States<S>[] states;

    //
    // constructors
    //

    /**
     * Creates a new FINITE horizon (MDP) Problem.
     * @param initial set of initial states
     * @param horizon last stage at which actions can be taken
     */

    @SuppressWarnings("unchecked")
    public FiniteMDP(States<S> initial, int horizon) {
        this.initial = initial;
        this.horizon = horizon;
        this.states = new States[horizon + 1];
        this.finite = true;
    }

    /**
     * Creates a finite horizon MDP.
     * @param initial a initial state
     * @param horizon horizon.
     */
    public FiniteMDP(S initial, int horizon) {
        this(new StatesSet<S>(initial), horizon);
    }

    //
    // Abstract Methods
    //
    @Override
    protected Solver<S, A> getDefaultSolver() {
        return new FiniteSolver<S, A>(this);
    }

    /**
     * This function must return the Immediate cost incurred when
     * taking action a from state i
     * @param i Current state
     * @param a Action
     * @param t Current time stage
     * @return Cost value
     */
    public abstract double immediateCost(S i, A a, int t);

    /**
     * This is the probability of going from state i to state j by
     * taking the action a at stage t.
     * @param i Current state
     * @param j Destination state
     * @param a Action taken
     * @param t Current time stage
     * @return Probability
     */

    public abstract double prob(S i, S j, A a, int t);

    /**
     * Set of States that can be reached from this state i, at this
     * stage t, after taking the acton a. The user must implement this
     * method.
     * @param i Current state
     * @param a Action taken
     * @param t Time stage
     * @return Set of reachable states.
     */
    public abstract States<S> reachable(S i, A a, int t);

    /**
     * All the states that are available at stage t.
     * @param t time stage
     * @return States available at stage t.
     */
    public States<S> getStates(int t) {
        if (states[t] == null) {
            int n = t;
            while ((states[n] == null) && (n > 1))
                n--;
            if (n > 1) {
                states[n + 1] = oneStageReachable(states[n], n);
            } else if (n == 1) {
                states[1] = oneStageReachable(initial, n);
            } else
                // n=0
                states[0] = initial;
            for (; n < t; n++) {
                states[n + 1] = oneStageReachable(states[n], n);
            }
        }
        return states[t];
    }

    private States<S> oneStageReachable(States<S> initSet, int t) {
        StatesSet<S> stts = new StatesSet<S>();
        for (S s : initSet) {
            Actions<A> act = feasibleActions(s, t);
            for (A a : act) {
                States<S> reached = reachable(s, a, t);
                for (S s1 : reached)
                    stts.add(s1);
            }
        }
        return stts;
    }


    /**
     * This method returns the cost incurred if the last stage ends
     * with the system at state i. The user must extend this method.
     * @param i Ending state
     * @return Cost.
     */
    public abstract double finalCost(S i);

    /**
     * This method returns the cost incurred if the last stage ends
     * with the system at state i choosing the best immediate cost.
     * @param i Ending State
     * @return Cost.
     */
    public final double defaultFinalCost(S i) {
        double val = 0, maxSoFar;
        Actions<A> availableActions = feasibleActions(i, horizon);
        maxSoFar = -Double.MAX_VALUE;
        for (A a : availableActions) {
            val = immediateCost(i, a, horizon);
            if (val > maxSoFar) {
                maxSoFar = val;
            }
        }
        return val;
    }

    /**
     * Returns the actions available at this state i and at this stage
     * t . The user must implement this method.
     * @param i Current State
     * @param t Time stage
     * @return Set of feasible actions.
     */
    public abstract Actions<A> feasibleActions(S i, int t);

    //
    // Other Methods
    //

    /**
     * Sets the time lastStage at which decisions can be taken
     */
    protected void setHorizon(int T) {
        horizon = T;
    }

    /**
     * Returns the time lastStage
     * @return Time horizon
     */

    public int getHorizon() {
        return (horizon);
    }

}// class end

