package jmarkov.jmdp;

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
 * This class represents an infinite horizon, discrete time, Markov Decision
 * Process with events, where actions depend on events. It allows the definition
 * of events that can occur in a given state and this makes the cost and
 * probability definition easier to define than in the cases where no events are
 * defined.
 * 
 * @author Andres Sarmiento and Germán Riaño. Universidad de los Andes.
 * @param <S>
 *            States class
 * @param <A>
 *            Action Class
 * @param <E>
 *            Events class
 * 
 */
public abstract class DTMDPEvA<S extends State, A extends Action, E extends Event>
        extends DTMDP<StateEvent<S, E>, A> {

    //
    // constructors
    //
    /**
     * Creates a new infinite horizon discrete time (MDP) Problem with events
     * 
     * @param initial
     *            set of initial states for the exploration algorithm
     */
    @SuppressWarnings("unchecked")
    public DTMDPEvA(States<S> initial) {
        super((States<StateEvent<S, E>>) new StatesSet<S>());
        StatesSet<StateEvent<S, E>> statesSet = new StatesSet<StateEvent<S, E>>();
        for (S s : initial) {
            for (A a : feasibleAct(s)) {
                for (E e2 : activeEvents(s, a))
                    statesSet.add(new StateEvent<S, E>(s, e2));
            }
        }
        this.initial = statesSet;
    }

    //
    // implemented methods
    //

    @Override
    public final States<StateEvent<S, E>> reachable(StateEvent<S, E> i, A a) {
        // TODO: usar iteradores
        StatesSet<StateEvent<S, E>> statesSet = new StatesSet<StateEvent<S, E>>();
        States<S> reached = reachable(i.getState(), a, i.getEvent());
        for (S j : reached) {
            for (E e2 : activeEvents(j, a))
                statesSet.add(new StateEvent<S, E>(j, e2));
        }
        return statesSet;
    }

    @Override
    public final double immediateCost(StateEvent<S, E> i, A a) {
        return prob(i.getState(), i.getEvent())
                * immediateCost(i.getState(), a, i.getEvent());
    }

    @Override
    public final double prob(StateEvent<S, E> i, StateEvent<S, E> j, A a) {
        // p(i', j', a) = p(i, j, a, i.e) * p(j.e)
        return prob(i.getState(), j.getState(), a, i.getEvent())
                * prob(i.getState(), i.getEvent());
    }

    @Override
    public final Actions<A> feasibleActions(StateEvent<S, E> i) {
        return feasibleAct(i.getState());
    }

    //
    // abstract methods
    //

    /**
     * Reward received when the current state is i, the action taken is a and
     * event e occurs.
     * 
     * @param i
     *            current state
     * @param a
     *            action taken
     * @param e
     *            event that occurs
     * @return reward
     */
    abstract public double immediateCost(S i, A a, E e);

    /**
     * Conditional destination probability. Probability of reaching state j
     * given that the current state is i, the action taken is a and the event
     * that occurs is e.
     * 
     * @param i
     *            current state
     * @param j
     *            state to reach
     * @param a
     *            action taken (given)
     * @param e
     *            event that occurs (given)
     * @return Conditional probability
     */
    abstract public double prob(S i, S j, A a, E e);

    /**
     * Conditional event probability. Probability that event e occurs given that
     * the current state is i.
     * 
     * @param i
     *            current state
     * @param e
     *            event that occurs
     * @return Conditional probability
     */
    abstract public double prob(S i, E e);

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
    abstract public States<S> reachable(S i, A a, E e);

    /**
     * Set of events that are active from state i given that action a is taken.
     * 
     * @param i
     *            current state
     * @param a
     *            action taken
     * @return set of events that can occur
     */
    abstract public Events<E> activeEvents(S i, A a);

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