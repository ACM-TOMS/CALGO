package jmarkov.jmdp.solvers;
/*
 * Created Nov 27 2005
 */
import jmarkov.basic.Action;
import jmarkov.basic.State;
import jmarkov.jmdp.DTMDP;
/**
 * Structural class to be extended by solvers in order to solve the total cost
 * criteria for an infinite horizon problem
 * 
 * @author Andres Sarmiento and Germán Riaño - Universidad de Los Andes
 *
 * @param <S> state
 * @param <A> action
 */
public abstract class AbstractTotalSolver<S extends State, A extends Action>
        extends AbstractInfiniteSolver<S,A> {

    /**
     * Creates a solver for a discrete time problem
     * @param problem discrete time problem
     */
    public AbstractTotalSolver(DTMDP<S,A> problem){
        super(problem);
    }
    
    

}
