package jphase;

import jmarkov.basic.JMarkovElement;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

/**
 *  This interface defines the methods and attributes that any Phase-Type
 *  distribution should have 
 *  @author  German Riaño
 *  @author  Juan F. Pérez
 *  @version  1.0
 */
public interface PhaseVar extends JMarkovElement{

	/**
	 * Returns the transition matrix of the Phase-Type Distribution
	 * @return Transition matrix for transient states of the 
	 * Phase-Type Distribution
	 */
	public Matrix getMatrix();

	/**
	 * Sets the transition matrix of the Phase-type distribution to be A
	 * @param A Transition matrix for the transient states of the 
	 * Phase-Type Distribution
	 */
	public void setMatrix(Matrix A);

	/**
	 *	Returns the initial probability mass vector
	 *  @return Initial probability mass vector
	 */
	public Vector getVector();


	/**
	 * Sets the initial probability vector to be alpha
	 * @param alpha Initial probability mass vector
	 */
	public void setVector(Vector alpha);

	/**
	 * Returns the number of phases of the Phase-type distribution
	 * @return Number of phases of the Phase-type distribution
	 */
	public int getNumPhases();
	
	
    /**
     * Returns the probability mass at zero (alpha_0)
     * @return Probability mass at zero (alpha_0)
     */
    public double getVec0();
    
    /**
     * Returns the exit vector from the transient states into absorption
     * @return Exit vector from the transient states into absorption
     */
    public Vector getMat0();
	
	/**
	 * Returns the transition matrix in double format
	 * @return Transition matrix for transient states of the 
	 * Phase-Type Distribution in double[][] format
	 */
	public double[][] getMatrixArray();
	
	/**
	 * Returns the initial probability 
	 * mass vector in double[] format
	 * @return Initial probability mass vector in double[] format
	 */
	public double[] getVectorArray();
	
	/**
	 * Returns the exit vector in double[] format
	 * @return Exit vector from the transient states to absorption
	 * in double[] format 
	 */
	public double[] getMat0Array();	
	
	/**
	 * Creates a deep copy of the original Phase-type variable
	 * @return A deep copy of the original Phase-type variable
	 */
	public PhaseVar copy();
    
	/**
	 * Computes the Expected Value of the Phase-type variable
	 * @return Expected Value of the Phase-type variable
	 */
	public double expectedValue();

	/**
	 * Computes the variance of the Phase-type variable
	 * @return Variance of the Phase-type variable
	 */
	public double variance();
	
	/**
	 * Computes the Standard deviation of the Phase-type variable
	 * @return Standard deviation of the Phase-type variable
	 */
	public double stdDeviation();
	
	/**
	 * Computes the Coefficient of Variation 
	 * of the Phase-type variable
	 * @return Coefficient of Variation of the Phase-type variable
	 */
	public double CV();
	
	/**
	 * Compuetes the k-th Moment of the Phase-type variable
	 * @param k Moment to compute
	 * @return k-th Moment of the Phase-type variable
	 */
	public double moment(int k);


	/**
     * Evaluates the cumulative distribution function at x
     * @param x Evaluation point
     * @return Cumulative density function at x
     */
    public double cdf(double x);
    
    /**
     * Evaluates the cumulative distribution function at n values of x, 
     * starting with x=0, step delta  
     * @param n number of evaluation points
     * @param delta distance between evaluation points
     * @return Evaluation of the survival Function at x = 0,d,2d,..,(n-1)d
     */
    public double[] cdf(int n, double delta);
    
    /**
     * Computes the probability that this variable 
     * takes a value between a and b
     * @param a inferior limit
     * @param b superior limit
     * @return Probability that this variable 
     * takes a value between a and b
     */
    public double prob(double a, double b);
    
    /**
     * Evaluates the survival function at x
     * @param x Evaluation point
     * @return Evaluation of the survival 
     * Function at x = 1-F(x)=P(X>x)
     */
    public double survival(double x);

    /**
     * Evaluates the survival function at n values of x, 
     * starting with x=0, step delta  
     * @param n number of evaluation points
     * @param delta distance between evaluation points
     * @return Evaluation of the survival Function at x = 0,d,2d,..,(n-1)d
     */
    public double[] survival(int n, double delta);
    
    /**
     * Evaluates the loss function of order 1
     * at x
     * @param x Evaluation point
     * @return Evaluation of the loss 
     * function of order 1
     */
    public double lossFunction1(double x);
    
    /**
     * Evaluates the loss function of order 2
     * at x
     * @param x Evaluation point
     * @return Evaluation of the loss 
     * function of order 2
     */
    public double lossFunction2(double x);
    
    /**
     * Computes the quantile q of the distribution, 
     * such that F(q) = p 
     * @param p probability such that F(q) = p
     * @return The quantile q of the distribution, 
     * such that F(q) = p
     */
    public double quantil(double p);

    /**
     * Computes the median of the distribution
     * @return The median of the distribution
     */
    public double median();



}
