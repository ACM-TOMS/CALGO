/**
 * TransitionsSet.java
 * Created: Feb 26, 2006
 */
package jmarkov.basic;

import java.util.Iterator;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 * @param <S>
 *            States class
 * 
 */
public final class TransitionsSet<S extends State> implements JMarkovElement,
        Transitions<S> {

    private SortedMap<S, Double> theSet = new TreeMap<S, Double>();

    /**
     * Default Constructor.
     */
    public TransitionsSet() {
    }

    /**
     * @see jmarkov.basic.Transitions#add(jmarkov.basic.Transition)
     */
    public boolean add(Transition<S> t) {
        return (theSet.put(t.getState(), t.getRate()) == null);
    }

    /**
     * Adds all the given Transtions to the current set.
     * 
     * @param trans
     *            A collection of Transitions
     * @return true if none of the elements was in the set.
     * @see java.util.Set#addAll(java.util.Collection)
     */
    public boolean add(Transitions<S> trans) {
        boolean result = true;
        for (Transition<S> t : trans) {
            result &= add(t);
        }
        return result;
    }

    /**
     * Adds a transition with the given state and rate.
     * @see jmarkov.basic.Transitions#add(State, double)
     */
    public boolean add(S state, double rate) {
        return add(new Transition<S>(state, rate));
    }

    /**
     * Adds the given rate to the transition to this state.
     * @param state 
     * @param rate 
     * @return The old value associated with this state.
     */
    public double addRate(S state, double rate) {
        double oldVal = getRate(state);
        add(new Transition<S>(state, oldVal +rate));
        return oldVal;
    }

    /**
     * Gets the rate for this state. It returns 0.0 if this state is not in the
     * Transitions.
     * 
     * @param state
     * @return The rate for this state
     */
    public double getRate(S state) {
        Double val = theSet.get(state);
        return (val == null) ? 0.0 : val;
    }

    /*
     * @see Transitions#size()
     */
    public int size() {
        return theSet.size();
    }

    /*
     * @see jmarkov.basic.JMarkovElement#label()
     */
    public String label() {
        String result = "";
        if (size() > 15)
            result = "Transtions with " + theSet.size() + " states.";
        else
            result = description();
        return result;
    }

    /*
     * @see jmarkov.basic.JMarkovElement#description()
     */
    public String description() {
        return theSet.toString();
    }

    @Override
    public final String toString() {
        return label();
    }

    /**
     * Returns an iterator used to walk through the Transitions.
     * @see java.lang.Iterable#iterator()
     */
    public Iterator<Transition<S>> iterator() {
        return new Iterator<Transition<S>>() {
            Iterator<Map.Entry<S, Double>> it = theSet.entrySet().iterator();

            public Transition<S> next() {
                Map.Entry<S, Double> me = it.next();
                return new Transition<S>(me.getKey(), me.getValue());
            }

            public void remove() {
                it.remove();
            }

            public boolean hasNext() {
                return it.hasNext();
            }
        };
        // return theSet.entrySet().iterator();
    }

}
