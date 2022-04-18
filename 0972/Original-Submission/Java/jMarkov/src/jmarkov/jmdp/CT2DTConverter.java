package jmarkov.jmdp;

import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.jmdp.solvers.AbstractDiscountedSolver;

/**
 * This class formulates a DTMDP equivalent to a CTMDP.
 * 
 * @author Andres Sarmiento and Germán Riaño - Universidad de Los Andes
 * 
 * @param <S>
 *            state
 * @param <A>
 *            action
 */
public class CT2DTConverter<S extends State, A extends Action> extends
        DTMDP<S, A> {

    private double interestRate = -1;
    /**
     * This  map stores the outbound rate of the activeState for all actions. This reduces solver's 
     * computation time per iteration. Each time activeState changes the TreeMap must be initialized 
     */
    private Map<A, Double> exitRates = null;
    /**
     * This  is the state currently being analyzed.  
     */
    private S activeState;
    /**
     * This  is the original problem  
     */
    private CTMDP<S, A> problem;

    
    /**
     * Constructor is not public because it should only be invoked by CTMDP in
     * this same package.
     * 
     * @param problem
     *            an infinite horizon continuous time problem
     */
    public CT2DTConverter(CTMDP<S, A> problem) {
        super(problem.initial);
        this.problem = problem;
        problem.setConverter(this);
        if (initial != null) {
            for (S s : initial) {
                activeState = s;
                break;
            }
        }
    }

    //
    // Own methods
    //

    /**
     * This method calculates the exit rate for a given state and action. It
     * sums all rates for all reachable states under that action.
     * 
     * @param i
     *            current state
     * @param a
     *            current action
     * @return The soujurn rate for a given state and action
     */
    public double exitRate(S i, A a) {
        if (i.equals(activeState) && exitRates.containsKey(a)) {
            double rate = exitRates.get(a);
            if (rate > 0)
                return rate;
        }
        activeState = i;
        double exitRate = 0;
        States<S> reached = problem.reachable(i, a);
        for (S s : reached)
            exitRate += problem.rate(i, s, a);
        exitRates.put(a, exitRate);
        return exitRate;
    }

    //
    // Implemented methods
    //

    @Override
    public final double immediateCost(S i, A a) {
        if (interestRate < 0.0) {
            if (!(getSolver() instanceof AbstractDiscountedSolver))
                interestRate = ((AbstractDiscountedSolver<S, A>) problem
                        .getSolver()).getInterestRate();
            else
                interestRate = 0.0; // criteria is total or average
        }
        // c+gamma/(beta+lambda_i(a))
        return problem.lumpCost(i, a) + problem.continuousCost(i, a)
                / (interestRate + exitRate(i, a));
    }

    @Override
    public final States<S> reachable(final S i, A a) {
        final States<S> reached = problem.reachable(i, a);
        exitRates = new TreeMap<A, Double>();// implies change in i state
        return new States<S>() {

            public Iterator<S> iterator() {
                return new Iterator<S>() {
                    private boolean wasAdded = false;
                    private boolean justAdded = false;
                    private S theNext = null;

                    public boolean hasNext() {
                        return (!wasAdded) || reached.iterator().hasNext();
                    }

                    public S next() {
                        // TODO needs testing!!
                        if (!justAdded)
                            theNext = reached.iterator().next();
                        if ((!wasAdded) && theNext.compareTo(i) == -1) {
                            wasAdded = true;
                            justAdded = true;
                            return i;
                        } else {
                            justAdded = false;
                            return theNext;
                        }
                    }

                    public void remove() {
                        throw new UnsupportedOperationException(
                                "Method not supported");
                    }
                };
            }

            public int size() {
                throw new UnsupportedOperationException("Method not supported");
            }

            public int numerateStates() {
                throw new UnsupportedOperationException("Method not supported");
            }

            public boolean isClosed() {
                return true;
            }

        };
    }

    @Override
    public final double prob(S i, S j, A a) {
        if (i != j) {
            return problem.rate(i, j, a) / problem.getMaxRate();
        }
        return 1.0 - exitRate(i, a) / problem.getMaxRate();
    }

    @Override
    public final Actions<A> feasibleActions(S i) {
        return problem.feasibleActions(i);
    }

    /**
     * @return Returns the exitRates.
     */
    private Map<A, Double> getExitRates() {
        // TODO: Shouldnt this be private (or not existent?)
        return exitRates;
    }

    /**
     * @param exitRates
     *            The exitRates to set.
     */
    void setExitRates(TreeMap<A, Double> exitRates) {
        // TODO: Shouldnt this be private (or not existent?)
        this.exitRates = exitRates;
    }

}