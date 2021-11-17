/**
 * Transition.java
 * Created: Feb 26, 2006
 */
package jmarkov.basic;

/**
 * This class represent a transition to a given state. It has an associated rate
 * and state.
 * 
 * @author German Riano. Universidad de los Andes. (C) 2006
 * @param <S>
 *            State Class
 * 
 */
public final class Transition<S extends State> implements JMarkovElement {

    private S state = null;
    private double rate;

    /**
     * Basic constructor.
     * 
     * @param state
     * @param rate
     */
    public Transition(S state, double rate) {
        super();
        this.rate = rate;
        this.state = state;
    }

    /*
     * @see jmarkov.basic.JMarkovElement#label()
     */
    public String label() {
        return state + " -> " + rate;
    }

    /*
     * @see jmarkov.basic.JMarkovElement#description()
     */
    public String description() {
        return "Transition to state: " + state + "at rate " + rate;
    }

    @Override
    public final String toString() {
        return label();
    }

    /**
     * Returns the rate.
     * 
     * @return Returns the rate.
     */
    public final double getRate() {
        return rate;
    }

    /**
     * Returns the state.
     * 
     * @return Returns the state.
     */
    public final S getState() {
        return state;
    }

    // /**
    // * @param state The state to set.
    // */
    // public final void setState(S state) {
    // this.state = state;
    // }

}
