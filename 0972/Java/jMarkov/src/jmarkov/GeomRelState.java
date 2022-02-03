/*
 * Created on 11-jul-2005
 */
package jmarkov;

import jmarkov.basic.State;

/**
 * This class is used to build destinations which are relative to a given
 * GeomState. The user shoulld not extend this class, but rather use it when
 * building destination states.
 * 
 * @author Julio Góez, Germán Riaño. Universidad de los Andes. (C) 2005
 * @param <Sub>
 *            Is the State that represents the Sub-states
 */
public final class GeomRelState<Sub extends State> extends State {

    /** Relitive Level */
    protected int rLevel;
    /** Whether it is boundary */
    protected boolean boundary = false;
    /**
     * subState represnts the background states in every rLevel.
     */
    protected Sub subState;

    /**
     * Creates a Non boundary GeomState with the given relative level rLevel,
     * and subState.
     * 
     * @param rLevel
     * @param subState
     */
    public GeomRelState(Sub subState, int rLevel) {
        super();
        this.rLevel = rLevel;
        this.boundary = false;
        this.subState = subState;
    }

    /**
     * Creates a boundary GeomState with the given relative level rLevel, and
     * subState.
     * 
     * @param subState
     */
    public GeomRelState(Sub subState) {
        super();
        this.rLevel = Integer.MIN_VALUE;
        this.boundary = true;
        this.subState = subState;
    }

    /**
     * @return Returns the rLevel.
     */
    public int getRelLevel() {
        return rLevel;
    }

    /**
     * This method determines fi the State is a boundary state.
     * 
     * @return Whether it is Boundary
     */
    public boolean isBoundary() {
        return boundary;
    }

    /**
     * @return Returns the subState.
     */
    public Sub getSubState() {
        return subState;
    }

    /**
     * Compares GeomStates according to rLevel first and then according to the
     * subStates comparator.
     * 
     * @param s
     *            state to compare to.
     */
    @Override
    @SuppressWarnings("unchecked")
    public int compareTo(State s) {
        GeomRelState gs = (GeomRelState) s;
        // First Compare rLevel
        if (rLevel < gs.getRelLevel()) {
            return -1;
        } else if (rLevel > gs.getRelLevel()) {
            return +1;
        }
        // compare the rest
        else {
            return getSubState().compareTo(gs.getSubState());
        }
    }

    /**
     * @see jmarkov.basic.State#computeMOPs(MarkovProcess)
     */
    @Override
    public void computeMOPs(MarkovProcess mp) {
        subState.computeMOPs(mp);
    }

    /**
     * @see jmarkov.basic.State#label()
     */
    @Override
    public String label() {
        return subState.label();
    }

    /**
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        return -1 <= rLevel && rLevel <= +1 && subState.isConsistent();
    }

}
