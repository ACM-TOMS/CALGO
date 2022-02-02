package jphase.fit;

/**
 * This class defines the behavior for a class that
 * implements a moment-matching algorithm for  
 * Continuous Phase-Type distributions
 * @author Juan F. Pérez
 * @version 1.0 
 */
public abstract class MomentsContPhaseFitter extends ContPhaseFitter{
	
	/**
	 * First power moment to match
	 */
	protected double m1;
	
	/**
	 * Second power moment to match
	 */
	protected double m2;
	
	/**
	 * Third power moment to match
	 */
	protected double m3;
	
	/**
	 * @see jphase.fit.ContPhaseFitter#ContPhaseFitter(double[])
	 */
	public MomentsContPhaseFitter(double[] data){
		super(data);
		this.m1 = FitterUtils.powerMomentK(data, 1);
		this.m2 = FitterUtils.powerMomentK(data, 2);
		this.m3 = FitterUtils.powerMomentK(data, 3);
	}
	
	/**
	 * 
	 * @param m1 first power moment to match
	 * @param m2 second power moment to match
	 * @param m3 third power moment to match
	 */
	public MomentsContPhaseFitter(double m1, double m2, double m3){
		super(null);
		this.m1 = m1;
		this.m2 = m2;
		this.m3 = m3;
	}
}
