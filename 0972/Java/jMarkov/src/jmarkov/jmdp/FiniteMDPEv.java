package jmarkov.jmdp;

import jmarkov.basic.Action;
import jmarkov.basic.Event;
import jmarkov.basic.Events;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;

/*
 *
 * Created 30/07/05
 *
 */
/**
 * This class represents a finite horizon discrete time MDP with events.
 * 
 * @author Andres Sarmiento and Germán Riaño - Universidad de Los Andes
 * @param <S>
 *            States class
 * @param <A>
 *            Actions class
 * @param <E>
 *            Events class
 */
public abstract class FiniteMDPEv<S extends State, A extends Action, E extends Event>
        extends FiniteMDP<S, A> {

    /**
     * @param initial
     * @param horizon
     */
    public FiniteMDPEv(States<S> initial, int horizon) {
        super(initial, horizon);
    }

    @Override
    public double immediateCost(S i, A a, int t) {
        double sum = 0;
        Events<E> eventSet = activeEvents(i, a, t);
        for (E e : eventSet) {
            sum += prob(i, e, t) * immediateCost(i, a, e, t);
        }
        return sum;
    }

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
     * @param t
     *            current stage
     * @return reward
     */
    abstract public double immediateCost(S i, A a, E e, int t);

    @Override
    public double prob(S i, S j, A a, int t) {
        double sum = 0;
        Events<E> eventSet = activeEvents(i, a, t);
        for (E e : eventSet) {
            sum += prob(i, j, a, e, t) * prob(i, e, t);
        }
        return sum;
    }

    /**
     * Conditional probability. Probability of reaching state j given that the
     * current state is i, the action taken is a and the event that occurs is e.
     * 
     * @param i
     *            current state
     * @param j
     *            state to reach
     * @param a
     *            action taken (given)
     * @param e
     *            event that occurs (given)
     * @param t
     *            current stage
     * @return conditional probability
     */
    abstract public double prob(S i, S j, A a, E e, int t);

    /**
     * Conditional probability. Probability that event e occurs given that the
     * current state is i.
     * 
     * @param i
     *            current state
     * @param e
     *            event that occurs
     * @param t
     *            current stage
     * @return Conditional probability
     */
    abstract public double prob(S i, E e, int t);

    @Override
    public States<S> reachable(S i, A a, int t) {
        // TODO: This can be done with java iterators
        Events<E> eventSet = activeEvents(i, a, t);
        StatesSet<S> statesSet = new StatesSet<S>();
        for (E e : eventSet) {
            States<S> reached = reachable(i, a, e, t);
            for (S s : reached)
                statesSet.add(s);
        }
        return statesSet;
    }

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
     * @param t
     *            current stage
     * @return set of reachable states.
     */
    abstract public States<S> reachable(S i, A a, E e, int t);

    /**
     * Set of events that are active from state i given that action a is taken.
     * 
     * @param i
     *            current state
     * @param a
     *            action taken
     * @param t
     *            current stage
     * @return set of events that can occur
     */
    abstract public Events<E> activeEvents(S i, A a, int t);

}
