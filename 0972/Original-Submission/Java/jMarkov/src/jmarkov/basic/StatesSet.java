/**
 * StatesSet.java
 * Created: Jul 11, 2005
 * This file replaces the one in jmdp.basic
 */
package jmarkov.basic;

import java.util.Iterator;
import java.util.SortedMap;
import java.util.TreeMap;

import jmarkov.jmdp.CTMDP;
import jmarkov.jmdp.DTMDP;

/**
 * This class represent a set of States. It is used as a convinence to
 * build multiple destinations in the method <code>reachable</code>.
 * @see DTMDP#reachable(State, Action)
 * @see CTMDP#reachable(State, Action)
 * @author Germán Riaño, Andrés Sarmiento. Universidad de los Andes.
 * @param <S> The States class for the states in this set.
 */
public class StatesSet<S extends State> implements States<S> {

    private SortedMap<S, S> theSet = new TreeMap<S, S>();
    private boolean closed = false;

    /**
     * Creates an empty set of States;
     */
    public StatesSet() {
    }

    /**
     * Creates set of States with only this State;
     * @param s The state to include in the set.
     */
    public StatesSet(S s) {
        this();
        add(s);
    }

    /**
     * Creates a set of objects S from the given States<S>.
     * @param states a set of States of type States.
     */
    public StatesSet(Iterable<S> states) {
        this();
        for (S s : states) {
            add(s);
        }
    }

    /**
     * Creates a set of objects S from a given set of States<S>.
     * @param states a set of States of type States.
     */
    public StatesSet(S[] states) {
        this();
        for (S s : states) {
            add(s);
        }
    }

    /**
     * Creates a set of objects S from the given States <S>.
     * @param states a set of States of type States.
     */

    public StatesSet(States<S> states) {
        this();
        for (S s : states) {
            add(s);
        }
    }

    /**
     * This method returns a safe way to walk through the states in a
     * particular set. Collections and their implementations (Set,
     * List, and Map) have iterators defined by default.
     * @return iterator over the states.
     */
    public final Iterator<S> iterator() {
        return theSet.values().iterator();
    }

    /**
     * Adds the State s to the set.
     * @param s The State to be added.
     * @return True if the set did not already contained this event.
     */
    public boolean add(S s) {
        if (closed)
            throw new RuntimeException(
                    "This set was closed and no new elemns can be added");
        return (theSet.put(s, s) == null);
    }

    /**
     * Adds the States in the iterator to the set.
     * @param states a set of States of type States.
     * @return True if the set did not contain ANY of the elements.
     */
    public boolean add(Iterable<S> states) {
        boolean result = true;
        for (S s : states) {
            result &= add(s);
        }
        return result;
    }

    /**
     * Adds the States in the iterator to the set.
     * @param states a set of States of type States.
     * @return True if the set did not contain ANY of the elements.
     */

    public boolean add(States<S> states) {
        boolean result = true;
        for (S s : states) {
            result &= add(s);
        }
        return result;
    }

    /**
     * Returns the element that is equal (according to equals() ) to
     * the given element.
     * @param state The given state
     * @return The state in the set, or null if it was not defined in the set.
     */
    public S get(S state) {
        return theSet.get(state);
    }

    /**
     * Removes an object from the set.
     * @param s The element to remove.
     * @return If the remove was successful (i.e. the element was in
     *         the set).
     */
    public boolean remove(S s) {
        return (theSet.remove(s) == null);
    }

    /**
     * Returns true if the set contains this State.
     * @param s A State
     * @return true if the state is contained in the set.
     */
    public boolean contains(S s) {
        return theSet.containsKey(s);
    }

    /**
     * This method numerates all states and returns the number of
     * states found. Afther this method is called it is illegal to add
     * more states to the set.
     * @return The number of states.
     */
    public int numerateStates() {
        int numStates = 0;
        for (S s : theSet.keySet())
            s.setIndex(numStates++);
        closed = true;
        return numStates;
    }

    /**
     * Returns an array with the States in the set.
     * @return An array representation of the states.
     */

    @SuppressWarnings("unchecked")
    public S[] toStateArray() {
        int size = theSet.size();
        S[] result = null;
        if (size != 0) {
            S elem = theSet.firstKey();
            result = (S[]) java.lang.reflect.Array.newInstance(elem.getClass(),
                    size);

            Iterator<S> it = theSet.values().iterator();
            for (int i = 0; i < size; i++)
                result[i] = it.next();
        }
        return result;
    }

    /**
     * @return the amount of states in the set.
     */

    public int size() {
        return theSet.size();
    }

    /*
     * @see java.lang.Object#toString()
     */
    @Override
    public String toString() {
        return theSet.toString();
    }

    /**
     * @return Returns the closed.
     */
    public boolean isClosed() {
        return closed;
    }

}
