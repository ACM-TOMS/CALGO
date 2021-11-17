package jmarkov.basic;

import java.util.Iterator;

/**
 * 
 * This class represents a set of objects Event. The user must choose his own
 * data structure and define the constuctors. For an easy way to declare and 
 * use a set of events see <\c>EventsCollection<\c>, which is an extension of Events.
 *  
 * 
 * @author Andres Sarmiento and Germán Riaño - Universidad de Los Andes
 * @param <E> The event Element.
 */
abstract public interface Events<E extends Event> extends Iterable<E> {

    /**
     * This function must be implemented. Must return an iterator over the
     * events.
     */
    abstract public Iterator<E> iterator();
    
    /**
     * This method adds an object to the set of events.
     * @param s object to be added.
     * @return True if the set did not contained this element.
     */
    abstract public boolean add(E s);
    /**
     * Returns the number of elements.
     * @return the number of Event elements.
     */
    public int size();
}
