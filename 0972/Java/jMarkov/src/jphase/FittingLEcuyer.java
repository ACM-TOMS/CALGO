package jphase;

import java.io.File;
import java.util.ArrayList;


import jphase.GUI.PhaseVarInfo;
import jphase.distributions.DistributionManager;
import jphase.distributions.IDistribution;
import jphase.values.FitResult;
import jphase.values.ReaderWriter;

import cern.colt.list.DoubleArrayList;

import umontreal.iro.lecuyer.gof.GofStat;
import umontreal.iro.lecuyer.probdist.ChiSquareDist;
import umontreal.iro.lecuyer.probdist.ContinuousDistribution;
import umontreal.iro.lecuyer.probdist.KolmogorovSmirnovDist;

/**
 * This class perform Goodness of Fit Tests
 * @author Andres Sarmiento Romero
 */
public class FittingLEcuyer {
	
	/**
	 * The is in charge of manage the distributions 
	 */
	private DistributionManager manager;
	
	/**
	 * Constructor of the class
	 */
	public FittingLEcuyer() {
		manager = new DistributionManager();
	}

	/**
	 * This method sorts a data List
	 */
	public static DoubleArrayList ordenamiento (DoubleArrayList unift){
		
        for( int i = unift.size( ); i > 0; i-- )
        {
            for( int j = 0; j < i - 1; j++ )
            {
                double p1 = unift.get( j );
                double p2 = unift.get( j + 1 );

                if( p1 > p2 )
                {
                	unift.set( j, p2 );
                	unift.set( j + 1, p1 );
                }
            }
        }
        return unift;
	}

	/**
	 * This method performs a GOF test and displays the result
	 * @param ruta Location of data
	 * @param n Number of groups for the Chi-Square test
	 * @param distribucion Distribution used in the initial hypothesis
	 * @param p Percetage of data when the GOF is finished
	 */
	public FitResult fit(int n, String ruta, String distribucion, double p ) {
		double[] d = ReaderWriter.saveFileToVector(ruta);
		DoubleArrayList data = new DoubleArrayList(d);
		ContinuousDistribution cont = ajusteDis(distribucion, data);
		DoubleArrayList unift = GofStat.unifTransform(data, cont);
		DoubleArrayList sortedData = ordenamiento(unift);
		
		int gl = n;
		
		int tamano = unift.size()/gl;
		
		double estim = GofStat.chi2Equal(unift, tamano);
		
		double[] ksEstim = GofStat.kolmogorovSmirnov(sortedData);
		
		ChiSquareDist dist = new ChiSquareDist(gl);
		KolmogorovSmirnovDist dist2 = new KolmogorovSmirnovDist(sortedData.size());
		
		double pvalue = 1 - dist.cdf(estim);
		double pvalueKS = 1 - dist2.cdf(ksEstim[2]);
		File f= new File(ruta);
		IDistribution distribu = manager.darDistribucion(distribucion, darArray(cont.getParams()));
		return new FitResult(estim, pvalue, pvalueKS, n,f.getName() , distribu , data.elements(), p);
		
	}

	/**
	 * This method convers a data vector to a data list
	 * @param params data Array
	 * @return data List
	 */
	private ArrayList<Double> darArray(double[] params) {
		ArrayList<Double> retorna = new ArrayList<Double>();
		for(int i = 0; i != params.length; i++){
			retorna.add(params[i]);
		}
		return retorna;
	}


	/**
	 * This method performs a GOF test and displays the result
	 * @param n Number of groups for the Chi-Square test
	 * @param d Data array used in the GOF test
	 * @param distribucion Name of thed distribution used in the initial hypothesis
	 * @param p Percetage of data when the GOF is finished
	 */
	public FitResult fit(int n, double[] d, String distribucion, double p) {
		DoubleArrayList data = new DoubleArrayList(d);
		ContinuousDistribution cont = ajusteDis(distribucion, data);
		DoubleArrayList unift = GofStat.unifTransform(data, cont);
		DoubleArrayList sortedData = ordenamiento(unift);
		
		int gl = n;
		
		int tamano = unift.size()/gl;
		
		double estim = GofStat.chi2Equal(unift, tamano);
		
		double[] ksEstim = GofStat.kolmogorovSmirnov(sortedData);
		
		ChiSquareDist dist = new ChiSquareDist(gl);
		KolmogorovSmirnovDist dist2 = new KolmogorovSmirnovDist(sortedData.size());
		
		double pvalue = 1 - dist.cdf(estim);
		double pvalueKS = 1 - dist2.cdf(ksEstim[2]);
		IDistribution distribu = manager.darDistribucion(distribucion, darArray(cont.getParams()));
		return new FitResult(estim, pvalue, pvalueKS, n, "null" , distribu , data.elements(), p);		
	}

	/**
	 * This method performs a GOF test and displays the result
	 * @param var Distribution used in the initial hypothesis
	 * @param file Location of data
	 * @param number Number of groups for the Chi-Square test
	 * @param p Percetage of data when the GOF is finished
	 */
	public FitResult fitPhases(PhaseVarInfo var, String file, int number, double p){
		double[] dataX = ReaderWriter.saveFileToVector(file);
		DoubleArrayList data = new DoubleArrayList(dataX);
		
		double[] obs = new double[number];
		double[] est = new double[number];		

		double min = dataX[0];
		double max = dataX[0];
		for(int j = 0 ; j != dataX.length; j++){
			double dato = dataX[j];
			if(dato < min)
				min = dato;
			if (dato > max)
				max = dato;
		}
		
		double inter = ( max - min ) / number;
		PhaseVar var0 = var.var;
		for (int i = 0; i != number; i++){
			double nmax = min + inter*(i+1);
			double cdf = var0.prob(inter*i,nmax);
			est[i] = cdf*dataX.length;
		}		
		
		for(int j = 0 ; j != dataX.length; j++){
			double dato = dataX[j];
			for (int i = 0; i != number; i++){
				if( min  + inter * i <= dato && dato <= min + inter * (i +1))
					obs[i] ++;
			}		
		}
		
		double sum = 0;
		double sum1 = 0;
		for(int i = 0; i != number ; i++){
			sum += (Math.pow(obs[i] - est [i], 2))/est[i];
			sum1 += (Math.pow(obs[i] - est [i], 2));
		}
		sum1 = sum1/number;
		
		ChiSquareDist dist = new ChiSquareDist(number);
		double pvalue = 0;
		
		if(!"NaN".equals("" + sum)){
			pvalue = 1 - dist.cdf(sum);
		}
		
		return new FitResult(sum1, pvalue, number, var, data.elements(), p);	
	}

	/**
	 * This method fits the a distribution to a data vector and returns the distribution
	 * @param distribucion Name of the desired distribution
	 * @param data data array
	 * @return the fitted distribution
	 */
	private ContinuousDistribution ajusteDis(String distribucion, DoubleArrayList data) {
		return manager.ajustarParametros(distribucion, data);
	}

	/**
	 * Returns the name of the distribution
	 * @return list of distributions
	 */
	public ArrayList<String> darDistribuciones() {
		return manager.darDistribuciones();
	}
}
