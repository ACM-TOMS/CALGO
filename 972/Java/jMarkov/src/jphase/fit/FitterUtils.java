package jphase.fit;

/**
 * This class contains a set of methods for common
 * computations in the PhaseFitter classes
 * @author Juan F. Pérez
 * @version 1.0 
 */
public class FitterUtils {
	
	
	/**
	 * Computes the k-th power moment of the data trace
	 * @param data data trace
	 * @param k power moment to be calculated (>= 1)
	 * @return Data k-th power Moment
	 */
	public static double powerMomentK(double[] data, int k){
		if(k>=1){
			double moment =0;
			for(int i=0;i<data.length;i++){
				if(k==1)moment+=data[i];
				else if(k==2)moment+=data[i]*data[i];
				else if(k==3)moment+=data[i]*data[i]*data[i];
				else moment+=Math.pow(data[i],k);
			}
			return moment/data.length;
		}else{
			System.out.println("The required moment is not well-defined. It must be a positive integer");
			return 0;
		}
	}
	
	
	/**
	 * Computes the k-th factorial moment of the data trace
	 * @param data data trace
	 * @param k factorial moment to be calculated (>= 1)
	 * @return Data k-th factorial Moment
	 */
	public static double factMomentK(int[] data, int k){
		if(k>=1){
			double moment =0;
			for(int i=0;i<data.length;i++){
				double temp = data[i];
				for(int j =1;j<k;j++)temp*=(data[i]-j);
				moment+=temp;
			}
			return moment/data.length;
		}else{
			System.out.println("The required moment is not well-defined. It must be a positive integer");
			return 0;
		}
	}
	
	/**
	 * Computes the floor of a double with a
	 * predefined precision
	 * @param x
	 * @param epsilon precision
	 * @return The floor of a double with the
	 * predefined precision
	 */
	public static double floor(double x, double epsilon){
		if(Math.abs( x - ((int)x + 1.0) )< epsilon){
			return ((int)x+1 );
		}else{
			return ((int)x);
		}
	}
	
	/**
	 * Computes the ceiling of a double with a
	 * predefined precision
	 * @param x
	 * @param epsilon precision
	 * @return The ceiling of a double with the
	 * predefined precision
	 */
	public static double ceil(double x, double epsilon){
		if(Math.abs( x - ((int)x ) )< epsilon){
			return ((int)x );
		}else{
			return ((int)x+1);
		}
	}
	
	/**
	 * Computes the square root of a double with a
	 * predefined precision
	 * @param x
	 * @param epsilon precision
	 * @return The square root of a double with the
	 * predefined precision
	 */
	public static double sqrt(double x, double epsilon){
		if(Math.abs( x )< epsilon){
			return 0;
		}else{
			return Math.sqrt(x);
		}
	}
	
}
