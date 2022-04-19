/**
 * SolverException.java
 * Created: Mar 10, 2006
 */
package jmarkov.basic.exceptions;

/**
 * This exception is thrown by solve methods.
 * @author German Riano. Universidad de los Andes. (C) 2006
 *
 */
public class SolverException extends Exception {

    /**        */
    private static final long serialVersionUID = -4635090230228775900L;


    /**
     * @param message
     */
    public SolverException(String message) {
        super(message);
    }

    /**
     * @param message
     * @param cause
     */
    public SolverException(String message, Throwable cause) {
        super(message, cause);
    }


}
