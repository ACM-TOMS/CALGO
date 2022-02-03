/**
 * NewSimpleMarkovProcess.java
 * Created: Feb 26, 2006
 */
package jmarkov;

import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.Transitions;
import jmarkov.basic.TransitionsSet;

/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 * @param <S>
 *            States Class.
 * @param <E>
 *            Events Class.
 * 
 */
public abstract class SimpleMarkovProcess<S extends State, E extends Event>
        extends MarkovProcess<S, E> {

    /**
     * @param i0
     * @param eSet
     * @param name
     */
    public SimpleMarkovProcess(S i0, EventsSet<E> eSet, String name) {
        super(i0, eSet, name);
    }

    /**
     * @param i0
     * @param eSet
     */
    public SimpleMarkovProcess(S i0, EventsSet<E> eSet) {
        super(i0, eSet);
    }

    /**
     * 
     */
    public SimpleMarkovProcess() {
        super();
    }

    /**
     * Determines if event e is active when the system is in state i. It has to
     * be implemented by a subclass.
     * 
     * @param i
     *            The current State
     * @param e
     *            The current Event.
     * @return True if the Event is Active.
     */
    public abstract boolean active(S i, E e);

    /**
     * Determines the destination set of States when events e occurs. It has to
     * be implemented by the subclass.
     * 
     * @param i
     *            current State.
     * @param e
     *            The Event that ocurred.
     * @return The destination States
     */
    public abstract States<S> dests(S i, E e);

    /**
     * Returns the rate to go from State i to j when Event e occurs. This is
     * called only if Event <code>e</code> is indeed active and j is a valid
     * destination.
     * 
     * @param i
     *            Current state
     * @param j
     *            Destination State
     * @param e
     *            The occuring event
     * @return The rate at which the system goes from i to j when e occurs.
     * 
     * @see #dests(State, Event)
     */
    public abstract double rate(S i, S j, E e);

    /**
     * This method calls active, dests and rate to create the set of
     * transitions. The user cannot override this method and would rarely call
     * it.
     * 
     * @see jmarkov.MarkovProcess#activeTransitions(State, Event)
     */
    @Override
    public final Transitions<S> activeTransitions(S i, E e) {
        Transitions<S> trans = new TransitionsSet<S>();
        boolean isActive = active(i, e);
        if (isActive) {
            debug(3, "Event " + e + ((isActive) ? " is " : " is not")
                    + " active");
            States<S> dests = dests(i, e);
            if (dests == null || dests.size() == 0)
                debug(0, "WARNNING: Event [" + e + "] in state " + i
                        + "has no destinations!");
            else {
                for (S j : dests) {
                    Thread.yield();
                    assert (j.isConsistent());
                    boolean result = trans.add(j, rate(i, j, e));
                    assert (result);
                } // end for
            } // end if

        }
        return trans;
    }
    //
    // /**
    // * Updates the current rate thru R(i,j) := R(i,j) + rate(e,i)
    // *
    // * @param i
    // * current State.
    // * @param j
    // * destination State.
    // * @param e
    // * event.
    // */
    // private void oldUpdateRates(Stte i, Stte j, Evt e) {
    // double curVal, newVal;
    // curVal = i.getRateToState(j);
    // newVal = curVal + rate(i, j, e);
    // debug(5, "Rate(" + i + "," + j + ") = " + newVal + " (It was " + curVal
    // + ")");
    // i.setRateToState(j, newVal);
    // }
    //
    // /**
    // * Generates a single state step. (i.e. completely explores an state). It
    // * assumes that there are indeed states needing to be explored
    // */
    // private void OldGenerateStep() {
    // Stte i;
    // Stte curj;
    // States<Stte> dests;
    // i = uncheckedStates.first();
    // if (reporter.getDebugLevel() == 1) {
    // if (cnt % 100 == 0) {
    // debug(1, "", true); // newline
    // debug(1, "", false, true); // indent
    // }
    // debug(1, ".", false);
    // }
    // cnt++;
    // debug(2, "Building State " + i);
    // debug(5, "S = " + statesSet);
    // debug(5, "U = " + uncheckedStates);
    // uncheckedStates.remove(i);
    // statesSet.add(i);
    // foundStates.put(i, i);
    // i.computeMOPs();
    // for (Evt e : theEvents) {
    // boolean isActive = active(i, e);
    // if (isActive) {
    // debug(3, "Event " + e + ((isActive) ? " is " : " is not")
    // + " active");
    // dests = dests(i, e);
    // if (dests == null)
    // debug(0, "WARNNING: Event [" + e + "] in state " + i
    // + "has no destinations!");
    // else {
    // for (Stte j : dests) {
    // Thread.yield();
    // debug(4, "Event [" + e + "]: Adding arc i = " + i
    // + " j= " + j);
    // curj = foundStates.get(j);
    // if (curj == null) { // not in the set
    // debug(5, "New state found j = " + j);
    // uncheckedStates.add(j);
    // foundStates.put(j, j);
    // // j.computeMOPs();
    // assert (j.isConsistent());
    // } else { // it is already in the set
    // j = curj; // we ensure it is the same
    // }
    // oldUpdateRates(i, j, e);
    // } // end for
    // } // end if
    // } // fi
    // } // end for
    // Thread.yield();
    // }

}
