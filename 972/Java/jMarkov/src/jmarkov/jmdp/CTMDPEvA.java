package jmarkov.jmdp;

import java.util.TreeMap;

import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.Event;
import jmarkov.basic.Events;
import jmarkov.basic.State;
import jmarkov.basic.StateEvent;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;

/*
 * Created 30/07/05
 */

/**
 * This class represents an Infinite horizon, continuous time Markov Decision
 * Process with events where actions depend on events. It allows the definition
 * of events that can occur in a given state and this makes the cost and
 * probability definition easier than in the cases where no events are defined.
 * 
 * @author Andres Sarmiento and Germán Riaño. Universidad de los Andes.
 * @param <S>
 *            The staes class
 * @param <A>
 *            Actions class
 * @param <E>
 *            Events class
 * 
 */
public abstract class CTMDPEvA<S extends State, A extends Action, E extends Event>
        extends CTMDP<StateEvent<S, E>, A> {


    /** Initial set of States. */
    protected States<S> initSet;

    //
    // constructors
    //
    /**
     * Creates a new continuous time infinite horizon MDP Problem with events
     * 
     * @param initial
     *            set of initial states for the exploration algorithm
     */
    public CTMDPEvA(States<S> initial) {
        super(null);
        initSet = initial;
    }

    //
    // implemented methods
    //

    @Override
    protected StatesSet<StateEvent<S, E>> generate() {
        // TODO: This method could be neater. Absorbing states not active
        States<StateEvent<S, E>> unexplored = null;
        StatesSet<StateEvent<S, E>> explored = new StatesSet<StateEvent<S, E>>();
        // absorbingStates = new StatesSet<StateEvent<S, E>>();
        if (initial == null) {
            StatesSet<StateEvent<S, E>> statesSet = new StatesSet<StateEvent<S, E>>();
            for (S s : this.initSet) {
                for (E e2 : activeEvents(s))
                    statesSet.add(new StateEvent<S, E>(s, e2));
            }
            this.initial = statesSet;
            unexplored = new StatesSet<StateEvent<S, E>>(this.initial);
            for (StateEvent<S, E> s : this.initial) {
                activeState = s;
                break;
            }

            long initialTime = System.currentTimeMillis();
            while (unexplored.size() > 0) {
                StatesSet<StateEvent<S, E>> stts = oneStageReachable(unexplored);
                for (StateEvent<S, E> s : unexplored) {
                    explored.add(s);
                }
                for (StateEvent<S, E> s : explored) {
                    stts.remove(s);
                }
                unexplored = stts;
            }
            explorationTime = System.currentTimeMillis() - initialTime;
            maxRate *= 1.1;
            debug(1, explored.size() + " states found.\n");
        }
        return explored;
    }

    // p(e|i,a) = sum(k in reachable(i,a), rate(i,k,a,e)) / exitRate(i,a)
    /**
     * Conditional probability. Probability that event e occurs given that the
     * current state is i.
     * 
     * @param i
     *            current state
     * @param e
     *            event that occurs
     * @return Conditional probability
     */
    private double prob(StateEvent<S, E> i, A a) {
        States<S> st = reached(i.getState(), a, i.getEvent());
        double sum = 0;
        for (S s : st) {
            sum += rate(i.getState(), s, a, i.getEvent());
        }
        return sum / converter.exitRate(i, a);
    }

    @Override
    public double lumpCost(StateEvent<S, E> i, A a) {
        return prob(i, a) * lumpCost(i.getState(), a, i.getEvent());
    }

    @Override
    public double continuousCost(StateEvent<S, E> i, A a) {
        return prob(i, a) * continuousCost(i.getState(), a, i.getEvent());
    }

    @Override
    public double rate(StateEvent<S, E> i, StateEvent<S, E> j, A a) {
        return rate(i.getState(), j.getState(), a, i.getEvent()) * prob(i, a);
    }

    @Override
    public States<StateEvent<S, E>> reachable(StateEvent<S, E> i, A a) {
        StatesSet<StateEvent<S, E>> statesSet = new StatesSet<StateEvent<S, E>>();
        States<S> reached = reached(i.getState(), a, i.getEvent());
        for (S j : reached) {
            for (E e2 : activeEvents(j))
                statesSet.add(new StateEvent<S, E>(j, e2));
        }
        // implies change in i state activeState = i;
        converter.setExitRates(new TreeMap<A, Double>());
        return statesSet;
    }

    @Override
    public Actions<A> feasibleActions(StateEvent<S, E> s) {
        return feasibleAct(s.getState());
    }

    //
    // abstract methods
    //

    /**
     * Rate. Rate of going of reaching state j given that the current state is
     * i, the action taken is a and the event that occurs is e.
     * 
     * @param i
     *            current state
     * @param j
     *            state to reach
     * @param a
     *            action taken (given)
     * @param e
     *            event that occurs (given)
     * @return Rate
     */
    abstract public double rate(S i, S j, A a, E e);

    /**
     * Set of reachable states from state i given that action a is taken and
     * event e occurs.
     * 
     * @param i
     *            current state
     * @param a
     *            action taken
     * @param e
     *            event that occurs
     * @return set of reachable states.
     */
    abstract public States<S> reached(S i, A a, E e);

    /**
     * Set of events that are active from state i given that action a is taken.
     * 
     * @param i
     *            current state
     * @return set of events that can occur
     */
    abstract public Events<E> activeEvents(S i);

    /**
     * Reward instantaneously gained in the moment when action a is taken from
     * state i.
     * 
     * @param i
     *            current state
     * @param a
     *            action taken
     * @param e
     *            event that occurs
     * @return instanteneous reward.
     */
    public abstract double lumpCost(S i, A a, E e);

    /**
     * Reward obtained continuously in time during the sojourn time in state i
     * until an action a is taken and a transition is triggered.
     * 
     * @param i
     *            current state
     * @param a
     *            action taken
     * @param e
     *            event that occurs
     * @return instanteneous reward.
     */
    public abstract double continuousCost(S i, A a, E e);

    /**
     * Returns the set of actions available at this state. The user must
     * implement this method.
     * 
     * @param i
     *            current state
     * @return set of feasible actions
     */

    public abstract Actions<A> feasibleAct(S i);

}
