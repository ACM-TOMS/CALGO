/**
 * ContContPhaseVar.java
 * Created: Feb 18, 2006
 */
package jphase;

/**
 * This interface defines the methods and attributes that any Continuous
 * Phase-Type distribution should have
 * @author Juan F. Perez
 * @author German Riano. Universidad de los Andes. (C) 2006
 */
public interface ContPhaseVar extends PhaseVar {


    /**
     * Returns the sum of a Phase-type-distributed number of Continuous Phase-type distributions
     * 
     * @param B
     *            Discrete-Phase Type Distribution that determines the number of
     *            Continuous Phase-Type Distributions to sum
     * @param res
     *            Continuous Phase Variable to store the resulting distribution
     * @return Sum of a Phase number of Continuous Phase-type distributions
     */
    public abstract ContPhaseVar sumPH(DiscPhaseVar B, ContPhaseVar res);

    /**
     * Returns the sum of a Phase-type-distributed number of Continuous Phase-type distributions
     * 
     * @param B
     *            Discrete-Phase Type Distribution that determines the number of
     *            Continuous Phase-Type Distributions to sum
     * @return Sum of a Phase number of Continuous Phase-type distributions
     */
    public abstract ContPhaseVar sumPH(DiscPhaseVar B);
    
    
    
	/**
     * Evaluates the probability density function at x
     * @param x Evaluation point
     * @return Probability density function at x
     */
    public abstract double pdf(double x);

    
    /**
     * Evaluates the probability density function at n values of x, 
     * starting with x=0, step delta  
     * @param n number of evaluation points
     * @param delta distance between evaluation points
     * @return Evaluation of the probability density function at x = 0,d,2d,..,(n-1)d
     */
    public abstract double[] pdf(int n, double delta);
    
    /**
     * Computes the sum of this variable and B
     * 
     * @param B Variable to sum to the original
     * @param res Variable to store the result
     * @return Sum of this variable and B
     */
    public ContPhaseVar sum(ContPhaseVar B, ContPhaseVar res);

    /**
     * Computes the sum of this variable and B
     * 
     * @param B Variable to sum to the original
     * @return Sum of this variable and B
     */
    public ContPhaseVar sum(ContPhaseVar B);
    
    /**
     * Returns the sum of a geometric number of independent copies of this
     * variable
     * 
     * @param p Parameter of the geometric variable
     * @return Sum of a geometric number of independent copies of this variable
     */
    public ContPhaseVar sumGeom(double p);
    
    /**
     * Computes the distribution of the mix: res = A*p + B*(1-p)
     * 
     * @param B Variable to mix with the original
     * @param p Portion of this variable in the mix (0<=p<=1)
     * @param res Variable to store the resulting distribution with the same
     *            number of phases of the original distribution
     * @return Distribution of the mix: res = A*p + B*(1-p)
     */
    public ContPhaseVar mix(double p, ContPhaseVar B, ContPhaseVar res);

    /**
     * Computes the distribution of the mix: res = A*p + B*(1-p)
     * 
     * @param B Variable to mix with the original
     * @param p Portion of this variable in the mix (0<=p<=1)
     * @return Distribution of the mix: res = A*p + B*(1-p)
     */
    public ContPhaseVar mix(double p, ContPhaseVar B);
    
    /**
     * Returns the minimum between the variable B and the original: res =
     * min(A,B)
     * 
     * @param B Variable to compare with the original
     * @param res Variable to store the resulting distribution
     * @return res = min(A,B)
     */
    public ContPhaseVar min(ContPhaseVar B, ContPhaseVar res);

    /**
     * Returns the minimum between the variable B and the original: res =
     * min(A,B)
     * 
     * @param B Variable to compare with the original
     * @return res = min(A,B)
     */
    public ContPhaseVar min(ContPhaseVar B);
       
    /**
     * Returns the maximum between the variable B and the original: res =
     * max(A,B)
     * 
     * @param B Variable to compare with the original
     * @param res Variable to store the resulting distribution
     * @return res = max(A,B)
     */
    public ContPhaseVar max(ContPhaseVar B, ContPhaseVar res);

    /**
     * Returns the maximum between the variable B and the original: res =
     * max(A,B)
     * 
     * @param B Variable to compare with the original
     * @return res = max(A,B)
     */
    public ContPhaseVar max(ContPhaseVar B);
    
    /**
     * Returns a Phase continuous variable that is the original one times c
     * @param c Scale factor to be applied to the original Phase continuous
     *            distribution
     * @return Phase continuous variable that is the original one times c
     */
    public ContPhaseVar times(double c);

    /**
     * Computes the Residual Time Distribution
     * @param x evaluation point
     * @return Distribution of P(X - tau <= x | X > tau
     */
    public abstract ContPhaseVar residualTime(double x);

    /**
     * Computes the Equilibrium Residual Distribution
     * @return Fo(x) = integ(0,t,(1 - F(t))) / E(X)
     */
    public abstract ContPhaseVar eqResidualTime();

    /**
     * Computes the distribution of the waiting time in queue
     * @param rho Server utilization
     * @return Phase Variable that describes the waiting time in Queue
     */
    public abstract ContPhaseVar waitingQ(double rho);
    
    /**
     * Computes the variable (X-a)+, i.e. the distribution
     * takes the value of the original distribution if it
     * is greater or equal to a. Otherwise, it is equal 
     * to null.
     * @param a Parameter for determining loss variable
     * @return Phase Variable that describes (X-a)+
     */
    public abstract ContPhaseVar residualVar(double a);

    /**
     * Creates a deep copy of the original Phase-Type Variable
     * @return A deep copy of the original Phase-Type Variable
     */
    public ContPhaseVar copy();
    
    /**
     * Creates a new variable of the same class of the original 
     * Continuous Phase-Type Variable
     * @param n number of Phases of the new Variable
     * @return A new variable of the same class of the original 
     * Continuous Phase-Type Variable
     */
    public ContPhaseVar newVar(int n);

    //@Override
    public abstract String toString();

}