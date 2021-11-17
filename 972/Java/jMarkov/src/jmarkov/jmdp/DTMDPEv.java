package jmarkov.jmdp;

import jmarkov.basic.Action;
import jmarkov.basic.Event;
import jmarkov.basic.Events;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;

/*
 * Created 30/07/05
 */

/**
 * This class represents an infinite horizon, discrete time, Markov Decision
 * Process with events. It allows the definition of events that can occur in a
 * given state and this makes the cost and probability definition easier to
 * define than in the cases where no events are defined.
 * 
 * @author Andres Sarmiento and German Riaño. Universidad de los Andes.
 * @param <S>
 *            The States class
 * @param <A>
 *            Actions class
 * @param <E>
 *            Events class
 * 
 */
public abstract class DTMDPEv<S extends State, A extends Action, E extends Event>
        extends DTMDP<S, A> {

    //
    // constructors
    //
    /**
     * Creates a new infinite horizon discrete time (MDP) Problem with events
     * 
     * @param initial
     *            set of initial states for the exploration algorithm
     */
    public DTMDPEv(States<S> initial) {
        super(initial);
    }

    //
    // implemented methods
    //

    @Override
    public final double immediateCost(S i, A a) {
        double sum = 0;
        Events<E> eventSet = activeEvents(i, a);
        for (E e : eventSet) {
            sum += prob(i, e) * immediateCost(i, a, e);
        }
        return sum;
    }

    @Override
    public final double prob(S i, S j, A a) {
        double sum = 0;
        Events<E> eventSet = activeEvents(i, a);
        for (E e : eventSet) {
            sum += prob(i, j, a, e) * prob(i, e);
        }
        return sum;
    }

    @Override
    public final States<S> reachable(S i, A a) {
        Events<E> eventSet = activeEvents(i, a);
        StatesSet<S> statesSet = new StatesSet<S>();
        for (E e : eventSet) {
            States<S> reached = reachable(i, a, e);
            for (S s : reached)
                statesSet.add(s);
        }
        return statesSet;
    }

    //
    // abstract methods
    //

    /**
     * Cost incurred received when the current state is i, the action taken is a
     * and event e occurs.
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
     * Conditional destination probability. Probability of reaching state j given that the
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
     * @return Conditional probability
     */
    abstract public double prob(S i, S j, A a, E e);

    /**
     * Conditional Event probability. Probability that event e occurs given that the
     * current state is i.
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

}
