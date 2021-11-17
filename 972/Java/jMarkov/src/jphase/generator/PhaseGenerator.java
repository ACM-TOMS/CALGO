package jphase.generator;

import jphase.PhaseVar;


/**
 * This abstract class defines the behavior of classes 
 * that implement Phase-Type random number generators
 * @author Juan F. Perez
 * @version 1.0 
 */
public abstract class PhaseGenerator {
	
	/**
	 * Phase variable from which the random numbers must be generated
	 */
	protected PhaseVar var;
	
	/**
	 * Constructs a new PhaseGenerator, initializing it  
	 * @param var variable from which the random numbers must be generated
	 */
	public PhaseGenerator(PhaseVar var){
		this.var=var;
		this.initialize();
	}
	
	
	/**
	 * @return A random number that has a probability
	 * distribution of Phase-Type
	 */
	public abstract double getRandom();
	
	/**
	 * @param num Number of variates to be generated
	 * @return A vector of random numbers that have a probability
	 * distribution of Phase-Type
	 */	
	public abstract double[] getRandom(int num);
	
	/**
	 * Initialize the cutoff values and the aliases for the initial
	 * probability distribution and the transition probability matrix
	 */
	protected abstract void initialize();
	
	/**
	 * @return Phase variable that is being used to generate the random
	 * numbers
	 */
	public PhaseVar getVar(){
		return this.var;		
	}
}
