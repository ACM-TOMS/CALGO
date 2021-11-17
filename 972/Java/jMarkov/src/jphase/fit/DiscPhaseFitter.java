package jphase.fit;

import jphase.DiscPhaseVar;

/**
 * This class defines the methods and attributes for any class implementing 
 * fitting methods for Discrete Phase-Type distributions 
 * @author Juan F. Pérez
 * @version 1.0
 *
 */
public abstract class DiscPhaseFitter implements PhaseFitter{
	
	/**
	 * Fitted Discrete Phase-Type variable 
	 */
	protected DiscPhaseVar var;
	
	/**
	 * Non-negative data trace from independent experiments
	 */
	protected int data[];
	
	/**
	 * Constructor of the fitter from a data array 
	 * @param data array with data to perform the fitting procedures
	 */
	public DiscPhaseFitter(int[] data){
		if(data==null)System.out.println("No data loaded");
		else{ 
			for(int i = 0; i < data.length; i++){
				if(data[i] < 0){
					System.out.println("Negative data. No data loaded");
					break;
				}
			}
			this.data = data;
		}
	}
	
	/**
	 * 
	 * @return -1 if there is no data associated to the algorithm,
	 * 0 if there has not been found a ContPhaseVar yet, or the 
	 * likelihood.
	 */
	public double getLogLikelihood(){
		if(data == null)return -1;
		else if(var == null){
			System.out.println("Not variable has been fitted yet");
			return 0;
		}else{
			double logLH = 0;
			for(int i = 0; i < data.length; i++)
				logLH+=Math.log(this.var.pmf(data[i]));
			return logLH;
		}
	}
	
	/**
     * @see jphase.fit.PhaseFitter#fit()
     */
	public abstract DiscPhaseVar fit();
}
