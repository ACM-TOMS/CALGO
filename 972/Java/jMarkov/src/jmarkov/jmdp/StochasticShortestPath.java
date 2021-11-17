/*
 * Created on 24/11/2004
 *
 */
package jmarkov.jmdp;

import jmarkov.basic.Action;
import jmarkov.basic.StateC;
import jmarkov.basic.States;
import jmarkov.basic.exceptions.StructureException;

/**
 * This class represents an infinite horizon shortest path problem. 
 * @author Juan F. Redondo, Andres Sarmiento, German Riaño -- - Universidad de los Andes
 * @param <S> States class
 * @param <A> Actions class
 */

public abstract class StochasticShortestPath<S extends StateC, A extends Action>
        extends DTMDP<S, A> {
    // States<StateC> states;

    /**
     * @param states Constructor
     */
    public StochasticShortestPath(States<S> states) {
        super(states);
        this.finite = false;
    }

    /**
     * This method was specially created to eliminate in a existent graph the
     * self-transition probabilities.
     * 
     * @param i
     * @param j
     * @param a
     * @return the modified probability
     * @throws StructureException
     */
    public double modifiedProb(S i, S j, A a) throws StructureException {
        if (i.equals(j)) {
            if (i.isTerminal() == false) {
                return 0;
            }
            return prob(i, j, a);
        }
        if (prob(i, i, a) == 1) {
            StructureException ex = new StructureException(
                    "Assumptions violation: This algorithm does not accept self-transition with probability equals 1");
            throw ex;
        }
        return prob(i, j, a) / (1 - prob(i, i, a));
    }

}
