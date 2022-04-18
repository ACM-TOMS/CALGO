package jphase.generator;
import java.util.Random;
import umontreal.iro.lecuyer.rng.RandomStream;
import umontreal.iro.lecuyer.rng.MRG32k3a;


import jphase.AbstractDiscPhaseVar;


/**
 * This class implements the algorithm proposed 
 * by Neuts and Pagano "generating Random Variates  
 * of Phase-Type", 1981. This is also based in the 
 * so called alias method to generate a variate from 
 * a discrete distribution. This class implements 
 * the algorithm for the discrete case.
 * @author Juan F. Pérez
 * @version 1.0 
 */

public class NeutsDiscPHGenerator extends PhaseGenerator {
	
	/**
	 * Complete initial probability vector 
	 */
	double[] alphaDist;
	
	/**
	 * Aliases for the initial probability vector 
	 */
	int[] alphaAlias;
	
	/**
	 * Cutoff values for the initial probability vector 
	 */
	double[] alphaCutoff;
	
	/**
	 * Complete transition matrix   
	 */
	double[][] matrixDist;
	
	/**
	 * Aliases for the probability transition matrix 
	 */
	int[][] matrixAlias;
	
	/**
	 * Cutoff values for the probability transition matrix 
	 */
	double[][] matrixCutoff;
	
	/**
	 * Random generator 
	 */
	RandomStream rand;
	
	/**
	 * 
	 * @see jphase.generator.PhaseGenerator#PhaseGenerator(jphase.PhaseVar)
	 */
	public NeutsDiscPHGenerator(AbstractDiscPhaseVar var){
		super(var);
	}
	
	@Override
	public double getRandom() {
		double x = 0;
		boolean absorbed = false;
		int estado = GeneratorUtils.getNumber(this.alphaDist,this.alphaAlias,this.alphaCutoff, rand);
		if(estado==0)absorbed=true;
		while(absorbed == false){
			x++;
			estado=GeneratorUtils.getNumber(this.matrixDist[estado-1],this.matrixAlias[estado-1],this.matrixCutoff[estado-1], rand);
			if(estado==0)absorbed=true;
		}
		return x;
	}

	@Override
	public double[] getRandom(int num) {
		double[] numbers = new double[num];
		for(int i=0; i < num; i++)numbers[i]=this.getRandom();
		return numbers;
	} 
	
	@Override
	protected void initialize(){
		rand = new MRG32k3a();
		int n = this.var.getNumPhases();
		this.alphaDist=new double[n+1];
		this.alphaAlias=new int[n+1];
		this.alphaCutoff=new double[n+1];
		this.matrixDist=new double[n][n+1];
		this.matrixAlias=new int[n][n+1];
		this.matrixCutoff=new double[n][n+1];
		for(int i = 0; i<n; i++)this.alphaDist[i+1]=this.var.getVector().get(i);
		this.alphaDist[0]=1-Math.min(GeneratorUtils.sum(this.alphaDist),1);
		GeneratorUtils.aliasCut(this.alphaDist,this.alphaAlias, this.alphaCutoff);
		
		for(int j = 0; j < n; j++ ){
			for(int i = 0; i<n; i++)this.matrixDist[j][i+1]=this.var.getMatrix().get(j,i);
			this.matrixDist[j][0]=1-Math.min(GeneratorUtils.sum(this.matrixDist[j]), 1);
			GeneratorUtils.aliasCut(this.matrixDist[j],this.matrixAlias[j], this.matrixCutoff[j]);
		}		
	}
}
