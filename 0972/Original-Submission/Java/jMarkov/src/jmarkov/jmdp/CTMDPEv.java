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
 * This class represents an Infinite horizon, continuous time Markov Decision
 * Process with events. It allows the definition of events that can occur in a
 * given state and this makes the reward and probability definition easier than
 * in the cases where no events are defined.
 * 
 * @author Andres Sarmiento and Germ�n Ria�o. Universidad de los Andes.
 * @param <S>
 *            The States class
 * @param <A>
 *            Tha Action class
 * @param <E>
 *            the Events class
 * 
 */
public abstract class CTMDPEv<S extends State, A extends Action, E extends Event>
        extends CTMDP<S, A> {

    //
    // constructors
    //

    /**
     * This constructor builds a continuous time MDP with events.
     * 
     * @param initial
     *            set of initial states for the exploration algorithm
     */
    public CTMDPEv(States<S> initial) {
        super(initial);
    }

    //
    // implemented methods
    //


    // p(e|i,a) = sum(k in reachable(i,a)) rate(i,k,a,e)) / exitRate(i,a)
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
    private double prob(S i, A a, E e) {
        States<S> st = reachable(i, a);
        double sum = 0;
        for (S s : st) {
            sum += rate(i, s, a, e);
        }
        return sum / converter.exitRate(i, a);
    }

    @Override
    public final double lumpCost(S i, A a) {
        double sum = 0;
        Events<E> eventSet = activeEvents(i, a);
        for (E e : eventSet) {
            sum += prob(i, a, e) * lumpCost(i, a, e);
        }
        return sum;
    }

    @Override
    public final double continuousCost(S i, A a) {
        double sum = 0;
        Events<E> eventSet = activeEvents(i, a);
        for (E e : eventSet) {
            sum += prob(i, a, e) * continuousCost(i, a, e);
        }
        return sum;
    }

    @Override
    public final double rate(S i, S j, A a) {
        double sum = 0;
        Events<E> eventSet = activeEvents(i, a);
        for (E e : eventSet) {
            sum += rate(i, j, a, e);
        }
        return sum;
    }

    @Override
    public final States<S> reachable(S i, A a) {
        Events<E> eventSet = activeEvents(i, a);
        // TODO:manejar esto con iteradores
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

}
