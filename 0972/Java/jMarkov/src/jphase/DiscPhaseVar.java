package jphase;


/**
 * This interface defines the behaviour that any Discrete
 * Phase-Type distribution should have
 * @author Juan F. Perez
 * @author German Riano. Universidad de los Andes. (C) 2006
 */
public interface DiscPhaseVar extends PhaseVar{

    /**
     * Returns the sum of a Phase-type-distributed  number of Discrete Phase-type
     * distributions
     * @param B Discrete-Phase Type Distribution that determines
     * the number of Discrete Phase-Type Distributions to sum
     * @param res Discrete Phase Variable to store the resulting
     * distribution 
     * @return Sum of a Phase number of Discrete Phase-type
     * distributions
     */
    public abstract DiscPhaseVar sumPH(DiscPhaseVar B,
            DiscPhaseVar res);

    /**
     * Returns the sum of a Phase-type-distributed number of Discrete Phase-type
     * distributions
     * @param B Discrete-Phase Type Distribution that determines
     * the number of Discrete Phase-Type Distributions to sum
     * @return Sum of a Phase number of Discrete Phase-type
     * distributions
     */
    public abstract DiscPhaseVar sumPH(DiscPhaseVar B);
    
    
    
    /**
     * Evaluates the probability mass function at k
     * @param k Evaluation point
     * @return Evaluation of the probability
     * mass function at k
     */
    public abstract double pmf(int k);

    /**
     * Evaluates the probability mass function at 
     * n values of x, from zero to n times delta 
     * @param n number of evaluation points
     * @param delta distance between evaluation points
     * @return Evaluation of the survival Function at x = 0,d,2d,..,(n-1)d
     */
    public abstract double[] pmf(int n, int delta);

    
    /**
     * Computes the sum of variables: res = A +B
     * @param B Variable to sum to the original 
     * @param res Variable to store the result
     * @return Sum of Variables: res = A +B
     */
    public DiscPhaseVar sum(DiscPhaseVar B, DiscPhaseVar res);
    
    /**
     * Computes the sum of variables: res = A +B
     * @param B Variable to sum to the original 
     * @return Sum of Variables: res = A +B
     */
    public DiscPhaseVar sum(DiscPhaseVar B);
    
    /**
     * Returns the sum of a geometric number of 
     * independent copies of this variable
     * @param p Parameter of the geometric variable 
     * @return Sum of a geometric number of 
     * independent copies of this variable
     */
    public DiscPhaseVar sumGeom(double p);
   
    /**
     * Computes the distribution of the mix: 
     * res = A*p + B*(1-p)
     * @param B Variable to mix with the original
     * @param p Portion of this variable in the mix (0<=p<=1)
     * @param res Variable to store the resulting distribution
     * with the same number of phases of the original 
     * distribution
     * @return Distribution of the mix: res = A*p + B*(1-p)
     */
    public DiscPhaseVar mix(double p, DiscPhaseVar B, DiscPhaseVar res);

    /**
     * Computes the distribution of the mix: 
     * res = A*p + B*(1-p)
     * @param B Variable to mix with the original
     * @param p Portion of this variable in the mix (0<=p<=1)
     * with the same number of phases of the original 
     * distribution
     * @return Distribution of the mix: res = A*p + B*(1-p)
     */
    public DiscPhaseVar mix(double p, DiscPhaseVar B);
    
    /**
     * Returns the minimum between the variable B and
     * the original: res = min(A,B)
     * @param B Variable to compare with the original 
     * @param res Variable to store the resulting distribution
     * @return res = min(A,B)
     */
    public DiscPhaseVar min(DiscPhaseVar B, DiscPhaseVar res);
    
    /**
     * Returns the minimum between the variable B and
     * the original: res = min(A,B)
     * @param B Variable to compare with the original 
     * @return res = min(A,B)
     */
    public DiscPhaseVar min(DiscPhaseVar B);
    
    /**
     * Returns the maximum between the variable B and
     * the original: res = max(A,B)
     * @param B Variable to compare with the original 
     * @param res Variable to store the resulting distribution
     * @return res = max(A,B)
     */
    public DiscPhaseVar max(DiscPhaseVar B, DiscPhaseVar res);
    
    /**
     * Returns the maximum between the variable B and
     * the original: res = max(A,B)
     * @param B Variable to compare with the original 
     * @return res = max(A,B)
     */
    public DiscPhaseVar max(DiscPhaseVar B);
    
    /**
     * Creates a deep copy of the original Phase-Type Variable
     * @return A deep copy of the original Phase-Type Variable
     */
    public DiscPhaseVar copy();
    
    /**
     * Creates a new variable of the same class of the original 
     * Discrete Phase-Type Variable
     * @param n number of Phases of the new Variable
     * @return A new variable of the same class of the original 
     * Discrete Phase-Type Variable
     */
    public DiscPhaseVar newVar(int n);
    
    //@Override
    public abstract String toString();

}