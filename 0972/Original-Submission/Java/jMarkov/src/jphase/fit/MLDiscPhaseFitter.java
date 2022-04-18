package jphase.fit;

/**
 * This class defines the behavior for a class that
 * implements a maximum likelihood algorithm for 
 * fitting data to a Discrete Phase-Type distribution
 * @author Juan F. Pérez
 * @version 1.0 
 */
public abstract class MLDiscPhaseFitter extends DiscPhaseFitter{
	
	/**
	 * Log-likelihood of the variable found
	 */
	double logLH = 0;
	
	/**
	 * @see jphase.fit.DiscPhaseFitter#DiscPhaseFitter(int[])
	 */
	public MLDiscPhaseFitter(int[] data){
		super(data);
	}
	
	@Override
	public double getLogLikelihood(){
		if(this.logLH != 0)return this.logLH;
		else{
			if(data == null)return -1;
			else if(var == null){
				System.out.println("No variable has been fitted yet");
				return 0;
			}else{
				double logLH = 0;
				for(int i = 0; i < data.length; i++)
					logLH+=Math.log(this.var.pmf(data[i]));
				return logLH;
			}
		}
	}
}