package jphase.distributions;

import java.util.ArrayList;
import org.jfree.data.xy.DefaultTableXYDataset;
import org.jfree.data.xy.XYSeries;

import umontreal.iro.lecuyer.probdist.ContinuousDistribution;
import umontreal.iro.lecuyer.randvar.RandomVariateGen;

/**
 * This is an abstract class that implements 
 * some of the methods that define a distribution
 * @author Andrés Sarmiento. Universidad de los Andes. (C) 2013
 */
public abstract class EDistribution implements IDistribution{
	
	/**
	 * The distribution object
	 */
	public ContinuousDistribution distribucion;
	
	/**
	 * List of generators
	 */
	public ArrayList<String> generadores;
	
	/**
	 * List of parameters
	 */
	public ArrayList<String> parametros;
	
	/**
	 * The actual generator
	 */
	public RandomVariateGen generador;
	
	@Override
    public String toString( ){
		return distribucion.toString();
	}

	/**
	 * Returns the name of the distribution parameters
	 */
	public ArrayList<String> darParametros() {
		return parametros;
	}

	/**
	 * Returns the LeCruyer distribuion object
	 */
	public ContinuousDistribution darDistribucion() {
		return distribucion;
	}

	/**
	 * Return the current random variate generator
	 */
	public RandomVariateGen darGenerador() {
		return generador;
	}
	
	/**
	 * Return the set of available generators
	 */
	public ArrayList<String> darGeneradores( ) {
		return generadores;
	}
	
	/**
	 * Returns the double valor rounding to cifSignifica number significative numbers
	 * @param valor value to be rounded
	 * @param cifSignifica number of significance
	 * @return rounded number
	 */
	public double nSig(double valor, int cifSignifica){
		double pow = Math.pow(10, cifSignifica);
		double ret =((int)(valor*pow))/pow ;
		return ret;
	}

	public DefaultTableXYDataset getPDF( double a, double b, int m ){
        XYSeries dataPDF = new XYSeries(getName() + " (PDF)", true, false);
		double dif = (b-a)/m;
		for(int i = 0 ; i != m; i++){
			double point = a + i *dif;
			double pdf = distribucion.density(point);
        	dataPDF.add(point, pdf);
		}
		
        DefaultTableXYDataset datasetPDF = new DefaultTableXYDataset(); 
        datasetPDF.addSeries(dataPDF);
        return datasetPDF;
        
//        return ChartFactory.createXYLineChart("PDF - Theorical Distribution",
//        		"x", "F(x)", datasetPDF, PlotOrientation.VERTICAL, 
//        		true, true, true);        
	}
}
