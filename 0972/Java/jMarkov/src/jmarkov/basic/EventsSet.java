/*
 * Created on 20/06/2003
 *
 */
package jmarkov.basic;

import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * This class represent a set of Events. The set should be build at
 * the beggining and should not be changed in any way afterwards.
 * @author Germán Riaño. Universidad de los Andes.
 * @param <E> The type of Events in this set.
 */
public class EventsSet<E extends Event> implements Events<E> {

    private SortedSet<E> theSet = new TreeSet<E>();

    /**
     * Creates an empty set of Events;
     */
    public EventsSet() {
    }

    /**
     * Creates an empty set of Events;
     * @param eventArray an array representation of the set.
     */
    public EventsSet(E[] eventArray) {
        for (E e : eventArray) {
            add(e);
        }
    }

    /**
     * This method returns a safe way to walk through the events in a
     * particular set. Collections and their implementations (Set,
     * List, and Map) have iterators defined by default.
     * @return iterator over the events.
     */
	public final Iterator<E> iterator() {
        return theSet.iterator();
    }

//    /**
//     * Creates a set with numE Events.
//     * @param numE The number of events.
//     * @return An event set with events clled
//     */
//    public static EventsSet<Event> createEventSet(int numE) {
//        EventsSet<Event> es = new EventsSet<Event>();
//        for (int i = 0; i < numE; i++) {
//            es.add(new Event());
//        }
//        return es;
//    }

    /**
     * Adds the Event e to the set.
     * @param e The event to be added.
     * @return True if the set did not already contained this event.
     */
	public boolean add(E e) {
        e.setSet(this);
        e.setIndex(theSet.size());
        return theSet.add(e);
    }

    /**
     * Returns true if the set contains this Event.
     * @param e The event.
     * @return whether the set contains this event.
     */
    public boolean contains(Event e) {
        return theSet.contains(e);
    }

    /**
     * Returns an array with the Events in the set.
     * @return array representation of the set.
     */
    @SuppressWarnings("unchecked")
    public E[] toEventArray() {
        int size = theSet.size();
        E[] result = null;
        if (size != 0) {
            E elem = theSet.first();
            result = (E[]) java.lang.reflect.Array.newInstance(elem.getClass(),
                    size);

            Iterator<E> it = theSet.iterator();
            for (int i = 0; i < size; i++)
                result[i] = it.next();
        }
        return result;
    }

    /*
     * @see java.lang.Object#toString()
     */
    @Override
    public String toString() {
        return theSet.toString();
    }

    /**
     * @return the amount of events in the set.
     */

	public int size() {
        return theSet.size();
    }

}
