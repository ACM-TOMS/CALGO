package jmarkov.basic;

import jmarkov.MarkovProcess;

/**
 * 
 * This class represents a state compounded of a state and an event. It is used
 * for state expansion for the problems where actions can depend on the event
 * that happens in a transition. Only future events that can occur from the
 * state state should be allowed as events event.
 * 
 * @author Andres Sarmiento, Germán Riaño - Universidad de Los Andes
 * @param <S>
 *            States set
 * @param <E>
 *            Event set
 * 
 */

public class StateEvent<S extends State, E extends Event> extends State {

    private S state;
    private E event;

    /**
     * Builds a new state with the event information
     * 
     * @param state
     *            state
     * @param event
     *            event
     */
    public StateEvent(S state, E event) {
        this.state = state;
        this.event = event;
    }

    /**
     * Gets the state.
     * @return the original state from the state
     */
    public S getState() {
        return state;
    }

    /**
     * Gets the event.
     * @return the original event from the state
     */
    public E getEvent() {
        return event;
    }

    @Override
    public String label() {
        return "( " + state.label() + " , " + event.label() + ")";
    }

    @Override
    public int compareTo(State i) {
        if (!(i instanceof StateEvent))
            throw new IllegalArgumentException("Comparing wrong type of Objects");
        StateEvent s1 = (StateEvent) i;
        if (state.compareTo(s1.getState()) == 0)
            return event.compareTo(s1.getEvent());
        return state.compareTo(s1.getState());
    }

 

    /**
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        return state.isConsistent();
    }

    /**
     * @see jmarkov.basic.State#computeMOPs(jmarkov.MarkovProcess)
     */
    @Override
    public void computeMOPs(MarkovProcess<?, ?> model) {
        state.computeMOPs(model);
        
    }

}
