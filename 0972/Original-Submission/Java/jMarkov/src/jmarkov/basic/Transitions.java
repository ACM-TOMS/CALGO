/**
 * Transitions.java
 * Created: Feb 26, 2006
 */
package jmarkov.basic;


/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 * 
 * @param <S>
 */
public interface Transitions<S extends State> extends JMarkovElement,
        Iterable<Transition<S>> {

    /**
     * @param t
     *            The trenasition to add.
     * @return true of the element was already in the set.
     */
    public abstract boolean add(Transition<S> t);

    /**
     * Adds a ne transition to the given state
     * 
     * @param state
     *            State the transition goes to
     * @param rate
     *            The rate at which this transition occurs.
     * @return true if the state was already on the set.
     */
    public abstract boolean add(S state, double rate);

    /**
     * Adds the given rate to the transition to this state.
     * @param state 
     * @param rate 
     * @return The old value associated with this state.
     */
    public double addRate(S state, double rate) ;
    
    
    /**
     * Gets the rate for this state. It returns 0.0 if this state is not in the
     * Transitions.
     * 
     * @param state
     * @return The rate for this state
     */
    public double getRate(S state) ;
    
    /**
     * Returns the number of Transtions represented by this object.
     * @return The number of Transitions.
     * @see java.util.Set#size()
     */
    public abstract int size();

}