/*
 * Created on 08-jul-2005
 */
package jmarkov.basic.exceptions;


/**
 * This Exception should be thrown by the SteadyStateSolver if it detects that
 * there is not a unique solution to the stationary probabilities. This occurs
 * when there are multiple closed communicating classes in the system, and
 * therefore the corresponding linear system has multiple solutions.
 * 
 * @author Germán Riaño. Universidad de los Andes.
 * @see jmarkov.solvers.SteadyStateSolver
 */
public class NotUnichainException extends SolverException {

    /**
     * 
     */
    private static final long serialVersionUID = 3470470902826744678L;

    /**
     * Default constructor.
     * @param message
     */
    public NotUnichainException(String message) {
        super(message);
    }

    /**
     * Constructor with cause.
     * @param message
     * @param cause
     */
    public NotUnichainException(String message, Throwable cause) {
        super(message, cause);
    }

}
