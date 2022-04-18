/*
 * Created on 7/12/2004
 *
 */
package jmarkov.basic.exceptions;

/**
 * This exception is produced in shortest path problems if the conditions for
 * convergence are not met.
 * 
 * @author Juan F Redondo
 */
public class StructureException extends SolverException {

    /**
     * 
     */
    private static final long serialVersionUID = 8715059036058763687L;

    /**
     * Default constructor.
     * @param message
     */
    public StructureException(String message) {
        super(message);
    }

}
