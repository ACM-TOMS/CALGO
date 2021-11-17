package jmarkov.basic;

import java.util.Iterator;

/**
 * 
 * This interface represents a set of objects State. The user must choose his
 * own data structure and define the constuctors, or provide a mechanism to
 * generate the sates on-the-fly. A convinience class, StatesSet, is provided if
 * the states are to be stored.
 * 
 * @author Andres Sarmiento and Germán Riaño - Universidad de Los Andes
 * @param <S>
 *            States set
 * @see StatesSet
 */
abstract public interface States<S extends State> extends Iterable<S> {

    /**
     * This function must be implemented. Must return an iterator over the
     * states.
     */
    abstract public Iterator<S> iterator();

    /**
     * This method adds an object to the set of states.
     * 
     * @param s
     *            object to be added.
     * @return True if the element was already in the set.
     */
   // abstract public boolean add(S s);

    /**
     * Removes an object from the set.
     * 
     * @param s
     *            State to be removed
     * @return True if the set contined the element.
     */
   // abstract public boolean remove(S s);

    /**
     * Returns the number of elements.
     * 
     * @return the number of State elements.
     */
    public int size();
    
    /**
     * This method numerates all states and returns the number of states found.
     * Afther this method is called it is illegal to add more states to the set.
     * 
     * @return The number of states.
     */
    public int numerateStates() ;
    
    /**
     * The set is closed if all elements have been added. 
     * @return true if the set is closed.
     */
    public boolean isClosed();

    
}
