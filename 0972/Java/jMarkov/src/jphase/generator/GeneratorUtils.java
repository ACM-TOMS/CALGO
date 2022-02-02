package jphase.generator;

import no.uib.cipr.matrix.Vector;
import java.util.Random;
import umontreal.iro.lecuyer.rng.RandomStream;


import jphase.PhaseVar;

/**
 * This class contains a set of methods with common 
 * computations for the PhaseGenerator classes
 * @author Juan F. Pérez
 * @version 1.0 
 */
public class GeneratorUtils {
	
	/**
	 * precision for calculations
	 */
	static final double precision = 10E-5;
	
	
	
	/**
	 * Returns the index of the minimum value in the data
	 * @param data array of data
	 * @return The index of the minimum value in the data
	 */
	public static int argmin(double[] data){
		int index=0;
		for(int i=0;i<data.length;i++)if(data[i]<data[index])index=i;
		return index;	
	}

	
	/**
	 * Returns the index of the maximum value in the data
	 * @param data array of data
	 * @return The index of the maximum value in data
	 */
	public static int argmax(double[] data){
		int index=0;
		for(int i=0;i<data.length;i++)if(data[i]>data[index])index=i;
		return index;	
	}
	
	/**
	 * Returns the sum of the elements of the data array
	 * @param data array of data
	 * @return Sum of the elements of the data array
	 */
	public static double sum(double[] data){
		double sum=0;
		for(int i=0;i<data.length;i++)sum+=data[i];
		return sum;	
	}
	
	/**
	 * Returns the sum of the absolute values of the elements 
	 * of the data array
	 * @param data array of data
	 * @return Sum of the absolute values of the elements 
	 * of the data array
	 */
	public static double sumAbs(double[] data){
		double sum=0;
		for(int i=0;i<data.length;i++)sum+=Math.abs(data[i]);
		return sum;	
	}
	
	
	/**
	 * This method generates the aliases and cutoff values 
	 * according to the distribution specified. This distribution
	 * is an MTJ Vector and may sum less than one (as in the 
	 * initial probability vector of a 
	 * Phase-type distribution). The distribution is then adjusted 
	 * by adding a position in zero that completes the
	 * probability mass to sum to zero.
	 * @param dist Distribution from which the aliases and cutoff vectors
	 * must be generated. Represented by an MTJ Vector. 
	 * @param alias vector of the distribution
	 * @param cutoff values of the distribution
	 */
	public static void aliasCut(Vector dist, int[] alias, double[] cutoff){
		int n = dist.size();
		double[] dist2= new double[n+1];
		for(int i =0; i<n;i++)dist2[i+1]=dist.get(i);
		dist2[0]=GeneratorUtils.sum(dist2);
		aliasCut(dist2, alias, cutoff);	
	} 
	
	
	/**
	 * This method generates the aliases and cutoff values 
	 * according to the distribution specified. 
	 * @param dist Distribution from which the aliases and cutoff vectors
	 * must be generated. Represented by array of doubles.
	 * @param alias vector of the distribution
	 * @param cutoff values of the distribution
	 */
	public static void aliasCut(double[] dist, int[] alias, double[] cutoff){
		if(dist.length==alias.length && alias.length==cutoff.length){
			int n = dist.length;
			double[] b = new double[n];
			int i=0;
			for(i=0;i<n;i++){
				alias[i]=i;
				cutoff[i]=0;
				b[i]=dist[i]-1/(new Integer(n)).doubleValue();
			}
			double c, d = 0;
			int k, m = 0;
			boolean ready = false;
			for(i=0;i<n && !ready ;i++){
				k=GeneratorUtils.argmin(b);
				c=b[k];
				m=GeneratorUtils.argmax(b);
				d=b[m];
				alias[k]=m;
				cutoff[k]=1+c*n;
				b[k]=0;
				b[m]=c+d;
				if(GeneratorUtils.sumAbs(b)<precision)ready=true;
			}
		}
	}
	
	/**
	 * Returns a random number with discrete distribution dist in {0,...,n}
	 * @param dist Discrete distribution in {0,...,n}
	 * @param alias Aliases of the distributions points
	 * @param cutoff Cutoff values to generate the random numbers
	 * @param rand Random type object to use as generator of random numbers
	 * @return A random number with discrete distribution dist in {0,...,n}
	 */
	public static int getNumber(double[] dist, int[] alias, double[] cutoff, RandomStream rand) {
		double u1 = rand.nextDouble();
		double u2 = rand.nextDouble();
		int entero = new Double(Math.floor((dist.length)*u1)).intValue();
		return u2 < cutoff[entero]? entero : alias[entero];
	}
	
	/**
	 * Returns a random number with Erlang(lambda, r) distribution.
	 * @param lambda Erlang rate
	 * @param r number of phases in the Erlang Distribution
	 * @param rand Random number source
	 * @return random number with Erlang(lambda, r) distribution. 
	 */
	public static double erlang(double lambda, int r, RandomStream rand){
		if(r>0){
		double x=rand.nextDouble();
		for(int i = 1; i < r; i++)x*=rand.nextDouble(); 
		return -(1/lambda)*Math.log(x);
		}else return 0;
	}
	
	/**
	 * This method implements the KS algorithm proposed by
	 * González, Sahni and Franta in "An efficient algorithm for the
	 * Kolmogorov-Smirnov and Lilliefors Tests" in ACM Transactions
	 * on Mathematical Software, Vol 3, No. 1, March 1977, pages 60-64.
	 * @param data data trace to be tested 
	 * @param var theoretical phase variable to be compared against the trace
	 * @return absolute maximum deviation from the data to the phase
	 * variable
	 */
	public static double algorKS(double[] data, PhaseVar var){
		double Kmax=0;
		int n = data.length;
		int[] NUM = new int[n+1];
		double[] MAX = new double[n+1];
		double[] MIN = new double[n+1];
		//Step 1: Initialize
		for(int i = 0; i<n+1;i++){
			MIN[i]=1;
			MAX[i]=0;
			NUM[i]=0;
		}
		//Step 2 Input observations into bins
		for(int i = 0; i<n;i++){
			double f = var.cdf(data[i]);
			double j = Math.ceil(f*n);//bin for data i
			int jd = new Double(j).intValue();
			NUM[jd]++; 
			if(MAX[jd]<f)MAX[jd]=f;
			if(MIN[jd]>f)MIN[jd]=f;
		}		
		//Step 3 Maximum positive and negative deviates
		double j = 0;
		double DP = 0;//Positive deviation
		double DN = 0;//Negative deviation
		for(int i = 0; i<n+1;i++){
			if(NUM[i]>0){
				double z = MIN[i]-j/n;
				if(z>DN)DN=z;
				j=j+NUM[i];
				z=j/n - MAX[i];
				if(z>DP)DP=z;
			}
		}
		//Step 4 Calculate max deviation
		Kmax = Math.max(DP,DN);	
		
		return Kmax;
	}
}
