package jmarkov.jmdp;

import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.DecisionRule;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.solvers.AbstractDiscountedSolver;
import jmarkov.jmdp.solvers.ProbabilitySolver;

/*
 * Created 07/08/05
 */

/**
 * This class represents a continuous time MDP. It should ONLY be used
 * in INFINITE horizon Problems. It must be extended in order to
 * represent the appropriate structure for each INFINITE horizon MDP
 * problem. The user must implement at least the functions that have
 * been declared abstract. It is also necessary to create one of the
 * extensions of the class Solver. By default, the program includes
 * PolicyIterationSolver and ValueIterationSolver classes to solve
 * infinite horizon problems. The FiniteSolver class is only for
 * finite horizon problems. To solve the problem follow the
 * instructions in each of the solvers's instructions.
 * @author Andres Sarmiento, Germán Riaño - Universidad de Los Andes
 * @param <S> The state class
 * @param <A> The Action class
 * @see jmarkov.jmdp.solvers.PolicyIterationSolver
 * @see jmarkov.jmdp.solvers.ValueIterationSolver
 */

public abstract class CTMDP<S extends State, A extends Action> extends
        InfiniteMDP<S, A> {
    /** Tha maxRate used for uniformization */
    protected double maxRate = -1;
    // TODO: These two fields can probably be eliminated. 
    //protected TreeMap<A,Double> exitRates = null;
    protected S activeState = null;
    
    /** The converter used to map the problem to a DTMDP. */
    protected CT2DTConverter<S, A> converter;

    //
    // constructors
    //

    /**
     * Creates a new continuous time infinite horizon MDP Problem.
     * @param initial set of initial states for the exploration
     *        algorithm
     */

    public CTMDP(States<S> initial) {
        super(initial);
        // TODO: what the heck is this??
        if (initial != null) {
            for (S s : initial) {
                activeState = s;
                break;
            }
        }

    }

    //
    // implemented methods
    //
    /**
     * Complete set of states explored
     * @return set of states explored
     */
    @Override
    public StatesSet<S> getAllStates() {
        if (states == null) {
            this.states = generate();
            converter.states = this.states;
        }
        return states;
    }

    /**
     * Sets the class in charge of making a DTMDP equivalent to the
     * CTMDP
     * @param converter class that makes a DTMDP equivalent to the
     *        CTMDP
     */
    public void setConverter(CT2DTConverter<S, A> converter) {
        this.converter = converter;
    }

    /**
     * @return maximum exit rate for all states and all actions
     */
    public double getMaxRate() {
        if (maxRate < 0)
            states = getAllStates();
        return maxRate;
    }

    /**
     * Finds the states reachable in one step.
     * @param initSet
     * @return States reachable from this set
     */
    protected StatesSet<S> oneStageReachable(States<S> initSet) {
        // TODO: change back to type States, check if it works with absorbing states
        StatesSet<S> stts = new StatesSet<S>();
        for (S i : initSet) {
            Actions<A> act = feasibleActions(i);
            for (A a : act) {
                States<S> reached = reachable(i, a);
                double tempRate = 0;
                for (S j : reached) {
                    stts.add(j);
                    tempRate += rate(i, j, a);
                    if (tempRate > maxRate)
                        maxRate = tempRate;
                }
                if (tempRate == 0) {
                    // hasAbsorbingState = true;
                    // absorbingStates.add(i);
                }
            }
        }
        return stts;
    }

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
        maxRate *= 1.1;
        debug(1, explored.size() + " states found.\n");
        numStates = explored.numerateStates();
        return explored;
    }

    /**
     * @return The steady state probability for each state
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
     * Solves the problem with the given interest rate
     * @param interestRate the interest rate parameter to solve the
     *        problem.
     * @return The Optimal solution.
     * @throws SolverException
     */

    @SuppressWarnings("unchecked")
    public Solution<S, A> solve(double interestRate) throws SolverException {
        if (!(getSolver() instanceof AbstractDiscountedSolver))
            setSolver(getDefaultDiscountedSolver(interestRate));
        ((AbstractDiscountedSolver) getSolver()).setInterestRate(interestRate);

        return getSolver().solve();
    }

    //
    // abstract methods
    //

    /**
     * Cost incurred instantaneously in the moment when action a is
     * taken from state i.
     * @param i State
     * @param a Action
     * @return Lump cost received.
     */
    public abstract double lumpCost(S i, A a);

    /**
     * Cost incurred continuously in time until the next transition
     * from state i given that action a is taken.
     * @param i State
     * @param a Action
     * @return Rate at which cost is incurred when action a is taken.
     */
    public abstract double continuousCost(S i, A a);

    /**
     * Set of States that can be reached from this state i, after
     * taking the action a.
     * @param i current State
     * @param a action taken
     * @return the reachable states.
     */
    abstract public States<S> reachable(S i, A a);

    /**
     * Rate of going from state i to state j by taking the action a
     * @param i current state
     * @param j Destination state.
     * @param a Action taken
     * @return The rate
     */
    public abstract double rate(S i, S j, A a);

}// class end

