package jphase.fit;

/**
 * This class defines the behavior for a class that
 * implements a moment-matching algorithm for 
 * Discrete Phase-Type distributions
 * @author Juan F. Pérez
 * @version 1.0 
 */
public abstract class MomentsDiscPhaseFitter extends DiscPhaseFitter{
	
	/**
	 * First factorial moment
	 */
	protected double m1;
	
	/**
	 * Second factorial moment
	 */
	protected double m2;
	
	/**
	 * Third factorial moment
	 */
	protected double m3;
	
	/**
	 * @see jphase.fit.DiscPhaseFitter#DiscPhaseFitter(int[])
	 */
	public MomentsDiscPhaseFitter(int[] data){
		super(data);
		this.m1 = FitterUtils.factMomentK(data, 1);
		this.m2 = FitterUtils.factMomentK(data, 2);
		this.m3 = FitterUtils.factMomentK(data, 3);
	}
	
	/**
	 * 
	 * @param m1 first factorial moment
	 * @param m2 second factorial moment
	 * @param m3 third factorial moment
	 */
	public MomentsDiscPhaseFitter(double m1, double m2, double m3){
		super(null);
		this.m1 = m1;
		this.m2 = m2;
		this.m3 = m3;
	}
}
