package jmarkov.basic;

/**
 * This class represents a single Action in Markov Decision Process (MDP). It
 * implements Comparable in order to be easily organized and searched.
 * 
 * @author Germán Riaño and Andres Sarmiento - Universidad de Los Andes
 * @see java.lang.Comparable
 * @see jmarkov.jmdp
 */

public abstract class Action implements Comparable<Action>, JMarkovElement {

    /**
     * The user MUST override this method to give a (hopefully short) label for
     * the Action.
     * 
     * @return short description of the Action.
     */
    public abstract String label();

    /**
     * The user SHOULD override this method to give a complete description for
     * the action.
     * 
     * @return short description of the state.
     */
    public String description() {
        return "Action = " + label();
    }

    /**
     * This calls label().
     */
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
        if (!(o instanceof Action))
            return false;
        Action a1 = (Action) o;
        return (compareTo(a1) == 0);
    }

}
