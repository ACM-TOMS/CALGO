package jmarkov.basic.exceptions;

/**
 * This Exception indicates that the transition probability matrix is not
 * stochastic for the state and action computed. The matrix is not stochastic
 * when the sum of the transition probabilities (row) is not 1.0.
 * 
 * @author Andrés sarmiento, Germán Riaño. Universidad de los Andes.
 * 
 */
public class NonStochasticException extends RuntimeException {

    private static final long serialVersionUID = 98657353068057784L;

    /**
     * Default constructor.
     * @param message
     */
    public NonStochasticException(String message) {
        super(message);
    }

}
