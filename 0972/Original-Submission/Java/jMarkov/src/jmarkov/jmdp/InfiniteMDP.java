package jmarkov.jmdp;

import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.ValueFunction;
import jmarkov.jmdp.solvers.AbstractAverageSolver;
import jmarkov.jmdp.solvers.AbstractDiscountedSolver;
import jmarkov.jmdp.solvers.AbstractInfiniteSolver;
import jmarkov.jmdp.solvers.ProbabilitySolver;
import jmarkov.jmdp.solvers.RelativeValueIterationSolver;
import jmarkov.jmdp.solvers.ValueIterationSolver;

/**
 * This class is a structural class. It represents a general
 * Infinite horizon MDP problem. It is extended for discrete and
 * continuous problems.
 * @author Daniel F. Silva, Andres Sarmiento, German Riano
 * @param <S> States class.
 * @param <A> Actions class.
 * @see jmarkov.jmdp.solvers.PolicyIterationSolver
 * @see jmarkov.jmdp.solvers.ValueIterationSolver
 */

public abstract class InfiniteMDP<S extends State, A extends Action> extends
        MDP<S, A> {
    /** set of states */
    protected StatesSet<S> states;
    /** The value function */
    protected ValueFunction<S> probability;
    // TODO: IMPORTANT: Remove the probability solver, can be handled from jMarkov
    /**
     * This solver is not really needed since JMarkov can handle the job
     */
    protected ProbabilitySolver<S, A> probabilitySolver = null;
    /** Whether an absorving state was found */
    protected boolean hasAbsorbingState = false;
    /** Set of absorving states. */
    protected States<S> absorbingStates;
    /** Time used to explore the system. */
    protected long explorationTime = 0;
    /** Number of states. Set when calling generate. */
    protected int numStates = -1;
    /** Number of actions. */
    protected int numActions = -1;

    //
    // constructors
    //

    /**
     * Creates a new INFINITE Dynamic Programming (DP) Problem.
     * @param initial set of initial states for the exploration
     *        algorithm
     */

    public InfiniteMDP(States<S> initial) {
        this.finite = false;
        this.initial = initial;
    }

    //
    // Abstract methods
    //

    /**
     * Returns the set of actions available at this state.
     * @param i Current State
     * @return set of Actions that can be taken at this state.
     */
    public abstract Actions<A> feasibleActions(S i);

    /**
     * Returns the number of states in the model. It causes the model
     * to be generated.
     * @return The number of ststes in the system.
     */
    public final int getNumStates() {
        if (numStates == -1)
            generate();
        return numStates;
    }

    /**
     * @return The set of states found.
     */
    abstract protected StatesSet<S> generate();

    //
    // Implemented Functions
    //

    // TODO: Check these warnings
    /**
     * Sets the interest rate to be used in the problem solving if the
     * objective is to minimze the discounted cost.
     * @param interestRate effective interest rate
     */
    @SuppressWarnings("unchecked")
    protected void setInterestRate(double interestRate) {
        if (!(getSolver() instanceof AbstractDiscountedSolver))
            setSolver(getDefaultDiscountedSolver(interestRate));
        ((AbstractDiscountedSolver) getSolver()).setInterestRate(interestRate);
        if (!(this instanceof CTMDP)) {
            // double maxRate = ((CTMDP) this).getMaxRate();
            ((AbstractDiscountedSolver) getSolver())
                    .setInterestRate(interestRate);
        } else if (!(this instanceof DTMDP)) {
            ((AbstractDiscountedSolver) getSolver())
                    .setInterestRate(interestRate);
        }
    }

    /**
     * @see jmarkov.jmdp.MDP#getDefaultSolver()
     */
    protected AbstractDiscountedSolver<S, A> getDefaultDiscountedSolver(
            double interestRate) {
        if (this instanceof CTMDP)
            return new ValueIterationSolver<S, A>((CTMDP<S, A>) this,
                    interestRate);
        else if (this instanceof DTMDP)
            return new ValueIterationSolver<S, A>((DTMDP<S, A>) this,
                    interestRate);
        return null;
    }

    /**
     * @see jmarkov.jmdp.MDP#getDefaultSolver()
     */
    protected AbstractAverageSolver<S, A> getDefaultAverageSolver() {
        if (this instanceof CTMDP)
            return new RelativeValueIterationSolver<S, A>((CTMDP<S, A>) this);
        else if (this instanceof DTMDP)
            return new RelativeValueIterationSolver<S, A>((DTMDP<S, A>) this);
        return null;
    }

    @Override
    protected AbstractInfiniteSolver<S, A> getDefaultSolver() {
        return getDefaultAverageSolver();
    }

    /**
     * Complete set of states explored
     * @return set of states explored
     */
    public StatesSet<S> getAllStates() {
        if (states == null)
            states = generate();
        if (!states.isClosed())
            states.numerateStates();
        return states;
    }

    /**
     * Complete set of actions available in any state
     * @return set of actions available in any state
     */
    public ActionsSet<A> getAllActions() {
	    StatesSet<S> stts = getAllStates();
	    ActionsSet<A> acts = new ActionsSet<A>();
	    for(S s: stts){
	    	for(A a: feasibleActions(s)){
	    		acts.add(a);
	    	}
	    }
	    numActions=acts.size();
	    return acts;
    }

    
    /**
     * Returns the number of actions in the model. It causes the model
     * to be generated.
     * @return The number of actions in all states combined (union, not sum)
     */
    public final int getNumActions() {
    	if (numStates == -1)
    		generate();
    	if (numActions == -1)
            getAllActions();
        return numActions;
    }
    
    /**
     * @return solver of the current problem
     */
    @Override
    public AbstractInfiniteSolver<S, A> getSolver() {
        return (AbstractInfiniteSolver<S, A>) super.getSolver();
    }

}// class end

