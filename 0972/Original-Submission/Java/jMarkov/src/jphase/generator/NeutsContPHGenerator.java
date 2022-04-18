package jphase.generator;

import java.util.Random;
import umontreal.iro.lecuyer.rng.RandomStream;
import umontreal.iro.lecuyer.rng.MRG32k3a;

import jphase.AbstractContPhaseVar;

/**
 * This class implements the algorithm proposed 
 * by Neuts and Pagano "Generating Random Variates  
 * of Phase-Type", 1981. This is also based in the 
 * so-called alias method to generate random variates from 
 * a discrete distribution. This class implements 
 * the algorithm for the continuous case.
 * @author Juan F. Pérez
 * @version 1.0 
 */
public class NeutsContPHGenerator extends PhaseGenerator {
	
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
		 * Outgoing rates of each state 
		 */
		double[] rates;
		
		/**
		 * Complete embedded transition matrix 
		 */
		double[][] matrixDist;
		
		
		/**
		 * Aliases for the embedded transition matrix 
		 */
		int[][] matrixAlias;
		
		
		/**
		 * Cutoff values for the embedded transition matrix 
		 */
		double[][] matrixCutoff;
		
		
		/**
		 * Random generator 
		 */
		//Random rand;
		RandomStream rand;
		
		/**
		 * @see jphase.generator.PhaseGenerator#PhaseGenerator(jphase.PhaseVar)
		 */
		public NeutsContPHGenerator(AbstractContPhaseVar var){
			super(var);
		}
		
		@Override
		public double getRandom() {
			boolean absorbed = false;
			double x = 0;
			int n = this.var.getNumPhases();
			int[] k = new int[n+1];//number of visits to each phase
			int estado = GeneratorUtils.getNumber(this.alphaDist,this.alphaAlias,this.alphaCutoff, this.rand);
			if(estado==0)absorbed=true;
			while(absorbed == false){
				k[estado]++;
				estado=GeneratorUtils.getNumber(this.matrixDist[estado-1],this.matrixAlias[estado-1],this.matrixCutoff[estado-1], this.rand);
				if(estado==0)absorbed=true;
			}
			for(int i=0;i < n; i++)if(k[i+1]>0)x+=GeneratorUtils.erlang(this.rates[i],k[i+1], this.rand);
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
			//rand = new Random();
			rand = new MRG32k3a();
			
			int n = this.var.getNumPhases();
			this.alphaDist=new double[n+1];
			this.alphaAlias=new int[n+1];
			this.alphaCutoff=new double[n+1];
			this.rates=new double[n];
			this.matrixDist=new double[n][n+1];
			this.matrixAlias=new int[n][n+1];
			this.matrixCutoff=new double[n][n+1];
			
			for(int i = 0; i<n; i++)this.alphaDist[i+1]=this.var.getVector().get(i);
			this.alphaDist[0]=1-Math.min(GeneratorUtils.sum(this.alphaDist),1);
			GeneratorUtils.aliasCut(this.alphaDist,this.alphaAlias, this.alphaCutoff);
			
			for(int j = 0; j < n; j++ ){
				this.rates[j]=-this.var.getMatrix().get(j,j);
				for(int i = 0; i<n; i++)this.matrixDist[j][i+1]=this.var.getMatrix().get(j,i)/this.rates[j];
				this.matrixDist[j][j+1]=0;
				this.matrixDist[j][0]=1-Math.min(GeneratorUtils.sum(this.matrixDist[j]), 1);
				GeneratorUtils.aliasCut(this.matrixDist[j],this.matrixAlias[j], this.matrixCutoff[j]);
			}
				
		}
}
