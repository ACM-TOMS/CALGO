package jmarkov.jmdp.solvers;

/*
 * Nov 27 - 2005
 */
import jmarkov.basic.Action;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.ValueFunction;
import jmarkov.jmdp.CTMDP;
import jmarkov.jmdp.DTMDP;

/**
 * This is a structural class that must be extended by classes solving the
 * dicounted cost minimization problem.
 * 
 * @author Andres Sarmiento and Germán Riaño - Universidad de Los Andes
 * 
 * @param <S>
 *            state
 * @param <A>
 *            action
 */
public abstract class AbstractDiscountedSolver<S extends State, A extends Action>
        extends AbstractInfiniteSolver<S, A> {
    /** The discount factor used to bring to present value from next period. */
    protected double discountFactor = -1;

    /**
     * @param problem
     *            The problem associated with this solver.
     * @param interestRate
     *            The interest rate per period
     */
    protected AbstractDiscountedSolver(DTMDP<S, A> problem, double interestRate) {
        super(problem);
        setInterestRate(interestRate);
    }

    /**
     * @param problem
     *            The problem associated with this solver.
     * @param interestRate
     *            The interest rate (nominal compunded continuously). For
     *            example, if you measure time in months, and the APR is A, then
     *            this rate satisfies exp(i/12) = 1 + A. Therefore i=12*ln(1+A).
     */
    protected AbstractDiscountedSolver(CTMDP<S, A> problem, double interestRate) {
        super(problem);
        this.discountFactor = problem.getMaxRate()
                / (interestRate + problem.getMaxRate());
    }

    /**
     * Sets a new Interest Rate
     * 
     * @param interestRate
     *            set.
     */
    public final void setInterestRate(double interestRate) {
        if (interestRate < 0)
            throw new IllegalArgumentException(
                    "The interest rate must be positive.");
        double newDiscountFactor = 1 / (1 + interestRate);
        if (newDiscountFactor != discountFactor) {
            policy = null; // invalidate current solution
            valueFunction = new ValueFunction<S>();
        }
        this.discountFactor = newDiscountFactor;
    }

    /**
     * Returns the current value of the discount factor.
     * 
     * @return The current value of the discount factor.
     */
    public double getInterestRate() {
        return (1.0 / discountFactor) - 1.0;
    }

    /**
     * @param discountFactor
     *            The discountFactor to set.
     */
    protected final void setDiscountFactor(double discountFactor) {
        if (this.discountFactor != discountFactor) {
            policy = null; // invalidate current solution
        }
        this.discountFactor = discountFactor;
    }
    
    /**
     * Expected value of valueFunction for the current state and a
     * specified action.
     * @param discountF is the rate for discounting from one period to
     *        another. It means how much less it would represent to
     *        receive one unit of the reward in the next period
     *        instead of receiving it in the present period.
     */
    protected final double future(S i, A a, double discountF,
            ValueFunction<S> vf) {
        double sum = 0.0;
        States<S> reachableStates = getDiscreteProblem().reachable(i, a);
        for (S j : reachableStates) {
            sum += getDiscreteProblem().prob(i, j, a) * vf.get(j);
        }
        return discountF * sum;
    }

    /**
     * Expected value of valueFunction for the current state and a
     * specified action.
     * @param i Ths State
     * @param a Action taken
     * @param discountF is the rate for discounting from one period to
     *        another. It means how much less it would represent to
     *        receive one unit of the reward in the next period
     *        instead of receiving it in the present period.
     * @return Expected value of valueFunction.
     */
    protected final double future(S i, A a, double discountF) {
        return future(i, a, discountF, getValueFunction());
    }


}
