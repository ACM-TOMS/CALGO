package jmarkov.jmdp.solvers;

import jmarkov.basic.Action;
import jmarkov.basic.State;
import jmarkov.jmdp.CTMDP;
import jmarkov.jmdp.DTMDP;
/**
 * Structural class for average cost solvers to extend. 
 * @author Andres Sarmiento and Germán Riaño - Universidad de Los Andes
 *
 * @param <S> state
 * @param <A> action
 */
public abstract class AbstractAverageSolver<S extends State, A extends Action>
        extends AbstractInfiniteSolver<S,A> {

    /**
     * Creates a solver for an infinite horizon discrete time MDP 
     * @param problem discrete time problem
     */
    protected AbstractAverageSolver(DTMDP<S,A> problem){
        super(problem);
    }

    /**
     * Creates a solver for an infinite horizon continuous time MDP
     * @param problem continuous time problem
     */
    protected AbstractAverageSolver(CTMDP<S,A> problem){
        super(problem);
    }
}
