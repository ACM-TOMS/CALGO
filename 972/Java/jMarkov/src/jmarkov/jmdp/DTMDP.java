package jmarkov.jmdp;

import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.DecisionRule;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.NonStochasticException;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.solvers.AbstractDiscountedSolver;
import jmarkov.jmdp.solvers.ProbabilitySolver;

/**
 * This class represents a discrete time infnite horizon MDP problem.
 * It must be extended in order to represent the appropriate structure
 * for each problem. The user must implement at least the functions
 * that have been declared abstract.
 * @author Andres Sarmiento, Germán Riaño - Universidad de Los Andes
 * @param <S> State class
 * @param <A> Action class
 * @see jmarkov.jmdp.solvers.PolicyIterationSolver
 * @see jmarkov.jmdp.solvers.ValueIterationSolver
 */

public abstract class DTMDP<S extends State, A extends Action> extends
        InfiniteMDP<S, A> {


    //
    // constructors
    //

    /**
     * Creates a new infinite horizon discrete time (MDP) Problem.
     * @param initial set of initial states for the exploration
     *        algorithm
     */

    public DTMDP(States<S> initial) {
        super(initial);
    }

    /**
     * Creates a new infinite horizon discrete time (MDP) Problem.
     * @param initial An initial state for the exploration algorithm
     */

    public DTMDP(S initial) {
        this(new StatesSet<S>(initial));
    }

    //
    // Abstract methods
    //

    /**
     * Cost incurred when taking action a from state i
     * @param i Current State
     * @param a Current Action
     * @return The cost incurred per transition
     */
    public abstract double immediateCost(S i, A a);

    /**
     * Probability of going from state i to state j by taking the
     * action a
     * @param i Current state.
     * @param j Destination State
     * @param a Action
     * @return The probability.
     */

    public abstract double prob(S i, S j, A a);

    /**
     * Set of states that can be reached from this state i, after
     * taking the action a.
     * @param i Current State
     * @param a Action taken
     * @return The reachable states.
     */
    abstract public States<S> reachable(S i, A a);

    //
    // Implemented methods
    //

    
    // TODO: Check absorbent states stuff
    /**
     * TFinds the states reached in one step.
     * @param initSet
     * @return Set of states reached in one step.
     */
    protected StatesSet<S> oneStageReachable(States<S> initSet) {
        StatesSet<S> stts = new StatesSet<S>();
        for (S i : initSet) {
            Actions<A> act = feasibleActions(i);
            for (A a : act) {
                States<S> reached = reachable(i, a);
                double sum = 0.0;
                for (S j : reached) {
                    stts.add(j);
                    // hasAbsorbingState = ((reached.size() == 1 &&
                    // j.equals(i))
                    // || hasAbsorbingState);
                    // absorbingStates.add(j);
                    sum += prob(i, j, a);
                }
                if (Math.abs(sum - 1.0) > 1.0E-5) {
                    throw (new NonStochasticException(
                            "Non stochastic row of matrix for state " + i
                                    + " and action " + a + ". Sum = " + sum));
                }
            }
        }
        return stts;
    }

    // TODO: This code could be neater.
    @Override
    protected StatesSet<S> generate() {
        StatesSet<S> unexplored = new StatesSet<S>(initial);
        StatesSet<S> explored = new StatesSet<S>();
        // absorbingStates = new StatesSet<S>();
        long initialTime = System.currentTimeMillis();
        while (unexplored.size() > 0) {
            StatesSet<S> stts = oneStageReachable(unexplored);
            for (S s : unexplored) {
                explored.add(s);
            }
            for (S s : explored) {
                stts.remove(s);
            }
            unexplored = stts;
        }
        explorationTime = System.currentTimeMillis() - initialTime;
        debug(1, explored.size() + " states found.\n");
        numStates = explored.numerateStates();
        return explored;
    }

   

    /**
     * @return a map with the probability for each state.
     * @throws SolverException
     */
    public ValueFunction<S> getSteadyStateProbabilities()
            throws SolverException {
        if (probabilitySolver == null) {
            DecisionRule<S, A> dr = getOptimalPolicy().getDecisionRule();
            probabilitySolver = new ProbabilitySolver<S, A>(this, dr);
        }
        if (!probabilitySolver.isSolved()) {
            probabilitySolver.solve();
            probability = probabilitySolver.getProbability();
        }
        return probability;
    }

    /**
     * @param solv Sets the solver that solves the steady state
     *        probabilities.
     */
    public void setProbabilitySolver(ProbabilitySolver<S, A> solv) {
        probabilitySolver = solv;
    }

    /**
     * Solves the problem with the given interest rate
     * @param interestRate the interest rate parameter to solve the
     *        problem.
     * @return The soultion to the problem.
     * @throws SolverException
     */
    @SuppressWarnings("unchecked")
    public final Solution<S, A> solve(double interestRate)
            throws SolverException {
        if (!(getSolver() instanceof AbstractDiscountedSolver))
            setSolver(getDefaultDiscountedSolver(interestRate));
        ((AbstractDiscountedSolver) getSolver()).setInterestRate(interestRate);
        // this.interestRate = interestRate;
        return getSolver().solve();
    }

}// class end

