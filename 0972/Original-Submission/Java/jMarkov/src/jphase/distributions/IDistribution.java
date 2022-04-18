package jphase.distributions;

import java.util.ArrayList;

import org.jfree.data.xy.DefaultTableXYDataset;

import cern.colt.list.DoubleArrayList;
import umontreal.iro.lecuyer.probdist.ContinuousDistribution;
import umontreal.iro.lecuyer.randvar.RandomVariateGen;
import umontreal.iro.lecuyer.rng.RandomStream;

/**
 * This interface defines the methods necessary to define a distribution
 * @author Andrés Sarmiento. Universidad de los Andes. (C) 2013
 */
public interface IDistribution {
	
	/**
	 * Return the generator for the distribution
	 * @return generator
	 */
	public ArrayList<String> darGeneradores();
	
	/**
	 * Returns the generator for the distribution
	 * @param s name of the U generator
	 * @param nombre name of the distribution
	 * @param paramatros parameters oh the distribution
	 * @return generator
	 */
	public RandomVariateGen darGenerador(RandomStream s, String nombre, ArrayList<Double> paramatros);
	
	/**
	 * Resturns the parameters of the distribution
	 * @return distribution parameters
	 */
	public ArrayList<String> darParametros();
	
	/**
	 * Returns the distribution
	 * @return the distrition
	 */
	public ContinuousDistribution darDistribucion();

	/**
	 * Returns the generator for the distribution
	 * @param s U generator
	 * @param modGenT Name of the distribution
	 * @return The generator
	 */
	public RandomVariateGen darGenerador(RandomStream s, String modGenT);

	/**
	 * Returns the actual distribution
	 * @return the actual distribution
	 */
	public RandomVariateGen darGenerador();

	/**
	 * Fits the parameters of the distribution to data
	 * @param data data
	 */
	public void ajustarParametros(DoubleArrayList data);
	
	/**
	 * Gets the moments from the distribution
	 * @return an array of the moments
	 */
	public double[] getMoments();
	
	/**
	 * Returns and string with the description of the distribution
	 * @return the description
	 */
	public String aString();
	
	/**
	 * Returns and dataset with the PDF of the distribution in different points, 
	 * this method is necessary for painting the distribution
	 * @param a minimum point
	 * @param b maximum point
	 * @param m number of x domain points
	 * @return The table
	 */
	public DefaultTableXYDataset getPDF( double a, double b, int m );
	
	/**
	 * Returns the name of the distribution
	 * @return the name of the distribution
	 */
	public String getName();
}
