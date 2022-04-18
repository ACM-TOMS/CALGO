package jmarkov.basic;

import java.util.Iterator;

/**
 * 
 * 
 * This interface represents a set of objects <c>Action</c>. The user must choose his own
 * data structure and define the constuctors. It´s recommended to use the Set
 * structure to avoid repeated actions. The ActionsSet class extends this
 * class and explotes the goodnesses of Collections. It is recommended to use
 * that class instead of this one for beginner users.
 * 
 * @author Andres Sarmiento, Germán Riaño. - Universidad de Los Andes
 * @param <A> The Action class.
 * 
 * @see java.lang.Iterable
 * @see jmarkov.basic.ActionsSet
 */

abstract public interface Actions<A extends Action> extends Iterable<A> {

    /**
     * This function must be implemented. Must return an iterator over the
     * Actions.
     */
    abstract public Iterator<A> iterator();
    /**
     * Returns the number of elements.
     * @return the number of elements.
     */
    public int size();
    
    /**
     * 
     * @return A string representation of the Actions in this set.
     */
    public String toString();    
}