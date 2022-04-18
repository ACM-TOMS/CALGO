/*
 * Created on 19/06/2003
 */
package jmarkov.basic;

/**
 * The class Event allows the user to define the implementation of the Events
 * that can alter the States of the Markov Chain.
 * 
 * @author Germán Riaño. Universidad de los Andes.
 * 
 */
public abstract class Event implements Comparable<Event>, JMarkovElement {

    private int index = -1; // EVENT NUMBER
    // the set where it belongs:
    private EventsSet eSet = null;

    /**
     * Sets a link to conatainer set
     * 
     * @param eSet
     *            The set where this events belong.
     */
    void setSet(EventsSet<? extends Event> eSet) {
        this.eSet = eSet;
    }

    /**
     * Returns the set of Events to which this event belongs.
     * 
     * @return the set to which this event belongs.
     */
    public EventsSet getSet() {
        return eSet;
    }

    /**
     * Gives the position of the Event in the Events set. Returns -1 if this
     * events has not yet been added to the set.
     * 
     * @return The position of the Event in the Events set. Returns -1 if this
     *         events has not yet been added to the set.
     */
    public int getIndex() {
        return index;
    }

    /**
     * sets the number when added to the set. Should be changed only by add in
     * class EventsSet.
     * 
     * @param num
     *            The number of this event on its set.
     */
    void setIndex(int num) {
        this.index = num;
    }

    /**
     * Returns positive if this Event has a higher number then the given event.
     * 
     * @see java.lang.Comparable#compareTo(Object)
     */
	public int compareTo(Event ev) {
        return this.index - ev.getIndex();
    }

    /**
     * If this function is not overriden by the user it returns the Event
     * number. The user should override to give a short label description of the
     * Event. It is highly recommended that the user overrides it to give a more
     * descriptive label to be used when reporting the occurrance rates of the
     * events.
     * 
     * @return A short string description of the Event.
     * @see #description()
     */

	public String label() {
        return "Event " + getIndex();
    }

    /**
     * It is highly recommended that the user overrides it to give a description
     * to be used when reporting the occurrance rates of the events, and GUI.
     * 
     * @return a String description
     * @see #label()
     */
	public String description() {
        return label();
    }

    @Override
    public final String toString() {
        return label();
    }

    /**
     * This method calls compareTo to check if the Action are equal.
     * 
     * @see java.lang.Object#equals(java.lang.Object)
     */
    @Override
    public final boolean equals(Object o) {
        if (!(o instanceof Event))
            return false;
        return (compareTo((Event) o) == 0);
    }
}
