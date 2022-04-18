/*
 * Created on 14/11/2004
 *
 */
package jmarkov.basic;

/**
 * State to model shortest path problems.
 * 
 * @author Juan F Redondo, German Riaño - Universidad de los Andes
 * 
 */
public abstract class StateC extends State {

    private boolean terminal = false;

    /**
     * @return Returns true if this a terminal state.
     */
    public final boolean isTerminal() {
        return terminal;
    }

    /**
     * Default constructor
     * 
     */
    public StateC() {
        this(false);
    }

    /**
     * General constructor. 
     * @param t
     *            Whether it is a terminal state or not.
     */

    public StateC(boolean t) {
        terminal = t;
    }

}
