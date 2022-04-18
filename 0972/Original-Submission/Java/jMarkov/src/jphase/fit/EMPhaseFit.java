package jphase.fit;

import java.util.Arrays;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;

import jphase.DenseContPhaseVar;
import jphase.MatrixUtils;

/**
 * This class implements the Maximum Likelihood method 
 * proposed by Asmussen, Nerman and Olsson in "Fitting
 * Phase-type Distributions via the EM algorithm", 1996.
 * The method returns a general continuous 
 * Phase-Type distributions.
 * @author Juan F. Pérez
 * @version 1.0 
 */
public class EMPhaseFit extends MLContPhaseFitter  {
	
	/**
	 * @see jphase.fit.MLContPhaseFitter#MLContPhaseFitter(double[])
	 */
	public EMPhaseFit(double[] data){
		super(data);
	}
	/**
	 * Precision for the convergence criterion
	 * in the algorithm
	 */
	public static double precision = 10E-6;
	
	/**
	 * 
	 */
	public static double logPrecision = 10E-5;
	
	/**
	 * Precision for the convergence criterion in
	 * the coefficient of variance
	 */
	public static double precisionParam = 10E-6;

	/**
	 * Constant to multiply the size of the data trace to 
	 * obtain the number of evaluation points.
	 * The number of evaluation points in the Runge-Kutta
	 * algorithm is the size of the data trace times evalPoints 
	 */
	public static int evalPoints = 10;
	
	/**
	 * Generator matrix of the Phase-type distribution found 
	 */
	private DenseMatrix A;
	
	/**
	 * Initial probability vector of the Phase-type distribution found 
	 */
	private DenseVector alpha;
	
	private DenseVector t;
	
	private int n;
	
	/**
	 * 
	 * @return Phase-Type variable with the best fit
	 */
	@Override
	public DenseContPhaseVar fit() {
		int N = 5;
		int best = 0;
		double[] LH = new double[N]; 
		DenseContPhaseVar[] vars = new DenseContPhaseVar[N];
		
		for(int k = 0; k < N; k++){
			System.out.println("\n\nITERATION: "+(k+1));
			System.out.println("N: "+(k+1));
			this.n=k+1;
			vars[k] = new DenseContPhaseVar(k+1);
			LH[k] = doFitN(data);
			vars[k] = new DenseContPhaseVar(this.alpha, this.A);
			System.out.println("Variable: "+vars[k].toString());
			System.out.println("Likelihood:\t"+LH[k]);
			System.out.println("Mean:\t"+vars[k].expectedValue());
			System.out.println("Variance:\t"+vars[k].variance());
			System.out.println("CV:\t"+vars[k].CV());
			if(LH[k]>LH[best])best=k;
		}
		System.out.println("The best variable found is: "
				+vars[best].toString());
		System.out.println("Likelihood: "+LH[best]);
		return vars[best];
	}
	
	
	
	/**
	 * Estimation method with a given number of phases
	 * @param data
	 * @return The log-likelihood
	 */
	public double doFitN(double[] data){
		Arrays.sort(data);

		//Initialization 
		this.A = new DenseMatrix(this.n, this.n);
		
		for(int i = 0; i < this.n; i++){
			double temp = 0;
			for(int j = 0; j < this.n ; j++){
				if(j != i)this.A.set(i, j, Math.pow(10,-j));
				temp+=Math.pow(10,-j);
			} 
			this.A.set(i, i, -temp);
		}
		
		this.alpha = new DenseVector(this.n);
		
		for(int i = 0; i < this.n; i++){
			this.alpha.set(i, (i+1.0)/(this.n*this.n));
		}
		
		this.t = (DenseVector)A.copy().mult(MatrixUtils.OnesVector(A.numColumns()), 
				new DenseVector(alpha.size())).scale(-1.0);
		
		int iter = 0;
		int maxIter = 10;
		boolean ready = false;
		DenseVector B = new DenseVector(this.n);
		DenseVector Z = new DenseVector(this.n);
		DenseMatrix N = new DenseMatrix(this.n, this.n+1);
		@SuppressWarnings("unused") double LHold;
        LHold = 0;
		double LHnew = 0;
		
		@SuppressWarnings("unused") double logLHold = 0;
		double logLHnew = 0;
		while(ready == false){
			iter ++;
			DenseMatrix newA = this.A.copy();
			DenseVector newAlpha = this.alpha.copy();
			DenseVector newT = this.t.copy();
			eStep(data, B, Z, N, newA, newAlpha, newT);
			LHold = LHnew;
			
			LHnew = mStep(data.length, B, Z, N, newA, newAlpha, newT);
			
			logLHold = logLHnew;
			logLHnew=0; 
			DenseContPhaseVar var = new DenseContPhaseVar(newAlpha, newA);
			for(int v=0; v < data.length; v++){
				logLHnew += Math.log(var.pdf(data[v]));
			}
			
			for(int i = 0; i< this.n; i++)System.out.print(B.get(i)+"\t");
			System.out.println();
			
			Matrix difMat = this.A.copy().add(-1, newA);
			double dif = Math.abs(difMat.norm(Matrix.Norm.Maxvalue));
						
			this.A.set(newA);
			this.alpha.set(newAlpha);
			this.t.set(newT);
			if(dif < logPrecision || iter > maxIter )ready = true;
		}
		return logLHnew;
	}

	/**
	 * E-Step of the EM method 
	 */
	private void eStep(double[] data, DenseVector B, DenseVector Z, DenseMatrix N, 
			DenseMatrix newA, DenseVector newAlpha, DenseVector newT) {
		double min = data[0];
		double max = data[data.length-1];
		int points = data.length*evalPoints; //number of points to evaluate
		double step = (max-min)/points;
		DenseVector c[] = new DenseVector[this.n];
		for(int i = 0; i < this.n; i++)c[i]=this.t.copy().zero();
		Solution[] res = RungeKutta(step, min, max,
							newAlpha.copy(), newT.copy(), c);

		Solution[] eval = new Solution[data.length]; //yv
		int k = 0;
		
		//Evaluate vector functions on sample points 
		for(int i = 0 ; i < data.length; i++){
			while(res[k].t < data[i] && k < points-1)k++;
			if(k==0){
				eval[i] = new Solution(data[i], res[k].a.copy(), res[k].b.copy(),res[k].c);
			}else if(res[k].t < data[i] && k == points-1){
				eval[i] = new Solution(data[i], res[res.length-1].a.copy(), 
						res[res.length-1].b.copy(),res[res.length-1].c);
			}else{
				double tTemp = data[i]; 
				double pondA = (data[i] - res[k-1].t)/step;
				double pondB = (res[k].t - data[i])/step;
				DenseVector aTemp = (DenseVector)res[k-1].a.copy().scale(pondA).add(pondB, res[k].a.copy() );
				DenseVector bTemp = (DenseVector)res[k-1].b.copy().scale(pondA).add(pondB, res[k].b.copy() );
				DenseVector cTemp[] = new DenseVector[this.n];
				for(int j = 0; j < this.n; j++)
					cTemp[j] = (DenseVector)res[k-1].c[j].copy().scale(pondA).add(pondB, res[k].c[j].copy() );
				eval[i] = new Solution(tTemp, aTemp, bTemp, cTemp);
			}
			 
		}
		
		//Computation of the expectation-based estimate(B,Z,N)
		double[][] Biv = new double[this.n][data.length];
		double[][] Ziv = new double[this.n][data.length];
		double[][][] Nijv = new double[this.n][this.n][data.length];
		double[][] Ni0v = new double[this.n][data.length];
		for(int v = 0; v < data.length; v++){
			double denom = newAlpha.dot(eval[v].b);
			for(int i =0; i< this.n; i++){
				Biv[i][v]=newAlpha.get(i)*eval[v].b.get(i)/denom;
				Ziv[i][v]=eval[v].c[i].get(i)/denom;
				Ni0v[i][v]=newT.get(i)*eval[v].a.get(i)/denom;
				for(int j = 0; j < this.n; j++)
					if(i!=j)Nijv[i][j][v]=newA.get(i,j)*eval[v].c[i].get(j)/denom;
			}
		}
		
		for(int i = 0; i < this.n; i++){
			double temp = 0;
			for(int v = 0; v < data.length; v++)temp+=Biv[i][v];
			temp = (Math.abs(temp) < precisionParam) ? 0.0 : temp; 
			B.set(i, temp);
						
			temp = 0;
			for(int v = 0; v < data.length; v++)temp+=Ziv[i][v];
			temp = (Math.abs(temp) < precisionParam) ? 0.0 : temp;
			Z.set(i, temp);
			
			temp = 0;
			for(int v = 0; v < data.length; v++)temp+=Ni0v[i][v];
			temp = (Math.abs(temp) < precisionParam) ? 0.0 : temp;
			N.set(i, 0, temp);
			
			for(int j=0; j< this.n ; j++){
				if(i!=j){
				temp = 0;
				for(int v = 0; v < data.length; v++)temp+=Nijv[i][j][v];
				temp = (Math.abs(temp) < precisionParam) ? 0.0 : temp;
				N.set(i, j+1, temp);
				}
			}
			
		}
	}
	
	/**
	 * M-step of the EM method
	 */
	private double mStep(int n, DenseVector B, DenseVector Z, DenseMatrix N, DenseMatrix newA, DenseVector newAlpha, DenseVector newT) {
		int p = B.size();
		for(int i = 0; i < p; i++)newAlpha.set(i, B.get(i)/n);
		for(int i = 0; i < p; i++)if(Z.get(i)!=0)newT.set(i, N.get(i,0)/Z.get(i));
		for(int i = 0; i < p; i++){
			for(int j = 1; j < p+1 ; j++){
				if((i+1) != j && Z.get(i)!=0)newA.set(i, j-1, N.get(i,j)/Z.get(i));
			}
		}
		for(int i = 0; i < p; i++){
			double temp = -newT.get(i);
			for(int j = 0; j < p; j++){
				if(i != j)temp -= newA.get(i, j);
			}
				
			newA.set(i, i, temp);
		}
		
		double LH=0;//log-likelihood
		for(int i = 0; i< this.n; i++){
			for(int j = 1; j< this.n+1; j++){
				if(i+1 != j && newA.get(i,j-1)>precision){
					LH += N.get(i,j) * Math.log( newA.get(i,j-1) );
				}
			}
			if(newT.get(i)>precision)LH += N.get(i,0) * Math.log( newT.get(i) );
		}
		for(int i =0; i< this.n; i++){
			LH += newA.get(i,i)*Z.get(i);
		}
		for(int i =0; i< this.n; i++){
			if(newAlpha.get(i)>precision)LH += B.get(i) * Math.log(newAlpha.get(i));
		}
		return LH;
	}

	/**
	 * Runge-Kutta estimation 
	 */
	private Solution[] RungeKutta(double step, double ini, double fin, 
			DenseVector a0, DenseVector b0, DenseVector[] c0){
		int N = (int)Math.round((fin-ini)/step);
		
		Solution[] res = new Solution[N+1];
		
		double t = ini;
		DenseVector a = a0.copy();
		DenseVector b = b0.copy();
		DenseVector[] c = new DenseVector[c0.length];
		copiarVector(c0,c,c0.length);
		
		res[0]=new Solution(t, a, b, c);
		
		DenseVector K1a, K2a, K3a, K4a;
		DenseVector K1b, K2b, K3b, K4b;
		DenseVector[] K1c, K2c, K3c, K4c;
		K1c = new DenseVector[this.n];
		K2c = new DenseVector[this.n];
		K3c = new DenseVector[this.n];
		K4c = new DenseVector[this.n];
		
		for(int i = 1; i <= N; i++){
			//K1 : h*function(t, w)
			DenseVector aTemp = a.copy();
			DenseVector bTemp = b.copy();
			DenseVector[] cTemp = new DenseVector[c.length];
			copiarVector(c,cTemp,c.length);
			function(aTemp, bTemp, cTemp);
			K1a = aTemp.copy().scale(step);
			K1b = bTemp.copy().scale(step);
			for(int j = 0; j < this.n;j++)K1c[j] = cTemp[j].copy().scale(step);
			
			//K2 : h*function(t+h/2, w+K1/2)
			aTemp = (DenseVector)a.copy().add(0.5,K1a);
			bTemp = (DenseVector)b.copy().add(0.5,K1b);
			cTemp = new DenseVector[c.length];
			copiarVector(c,cTemp,c.length);
			for(int j = 0; j < this.n;j++){
				cTemp[j] = (DenseVector)c[j].copy().add(0.5,K1c[j]);
			}
			function(aTemp, bTemp, cTemp);
			K2a = aTemp.copy().scale(step);
			K2b = bTemp.copy().scale(step);
			for(int j = 0; j < this.n;j++)K2c[j] = cTemp[j].copy().scale(step);
			
			//K3 : h*function(t+h/2, w+K2/2)
			aTemp = (DenseVector)a.copy().add(0.5,K2a);
			bTemp = (DenseVector)b.copy().add(0.5,K2b);
			cTemp = new DenseVector[c.length];
			copiarVector(c,cTemp,c.length);
			for(int j = 0; j < this.n;j++){
				cTemp[j] = (DenseVector)c[j].copy().add(0.5,K2c[j]);
			}
			function(aTemp, bTemp, cTemp);
			K3a = aTemp.copy().scale(step);
			K3b = bTemp.copy().scale(step);
			for(int j = 0; j < this.n;j++)K3c[j] = cTemp[j].copy().scale(step);
			
			//K4 : h*function(t+h, w+K3);
			aTemp = (DenseVector)a.copy().add(K3a);
			bTemp = (DenseVector)b.copy().add(K3b);
			cTemp = new DenseVector[c.length];
			copiarVector(c,cTemp,c.length);
			for(int j = 0; j < this.n;j++){
				cTemp[j] = (DenseVector)c[j].copy().add(K3c[j]);
			}
			function(aTemp, bTemp, cTemp);
			K4a = aTemp.copy().scale(step);
			K4b = bTemp.copy().scale(step);
			for(int j = 0; j < this.n;j++)K4c[j] = cTemp[j].copy().scale(step);
			
			a = (DenseVector)a.add(1.0/6.0,K1a).add(1.0/3.0,K2a)
							.add(1.0/3.0,K3a).add(1.0/6.0,K4a);
			b = (DenseVector)b.add(1.0/6.0,K1b).add(1.0/3.0,K2b)
							.add(1.0/3.0,K3b).add(1.0/6.0,K4b);
			for(int j = 0; j < this.n; j++){
				c[j] = (DenseVector)c[j].add(1.0/6.0,K1c[j]).add(1.0/3.0,K2c[j])
				.add(1.0/3.0,K3c[j]).add(1.0/6.0,K4c[j]);
			}
							
			t = ini + i*step;
			res[i]=new Solution(t, a, b, c);
		}
		return res;
	}
	
	/**
	 * @param c0 origin
	 * @param c destination
	 * @param length
	 */
	private void copiarVector(DenseVector[] c0, DenseVector[] c, int length) {
		for(int k=0; k < length; k++){
			c[k]=c0[k].copy();
			
		}
	}
	
	private void function(DenseVector a, DenseVector b, DenseVector[] c){
		for(int i = 0; i < this.n; i++){
			c[i].set( this.A.copy().mult( c[i].copy(),  c[i].copy().zero() ).add(a.get(i),this.t.copy()) );
			for(int j = 0; j < this.n; j++)if(Math.abs(c[i].get(j)) < precisionParam)c[i].set(j, 0.0);
		}
		a.set(this.A.copy().transMult(a.copy(),   a.copy().zero()));
		for(int j = 0; j < this.n; j++)if(Math.abs(a.get(j)) < precisionParam)a.set(j, 0.0);
		b.set(  this.A.copy().mult( b.copy(),      b.copy().zero())  );
		for(int j = 0; j < this.n; j++)if(Math.abs(b.get(j)) < precisionParam)b.set(j, 0.0);
	}
	
	
}



	/**
	 * Solution of the Runge-Kutta evaluation 
	 */
	class Solution{
		double t;
		DenseVector a;
		DenseVector b;
		DenseVector[] c;
		
		Solution(double t, DenseVector a, DenseVector b, DenseVector[] c){
			this.t = t;
			this.a = a.copy();
			this.b = b.copy();
			this.c = new DenseVector[c.length];
			copiarVector(c,this.c,c.length);
	}
	
    @Override
	public String toString(){
		String s = "";
		s+=String.format("%5.3f", t);
		for(int i = 0; i < a.size(); i++)s+=String.format("\t%7.6f", a.get(i));
		for(int i = 0; i < b.size(); i++)s+=String.format("\t%7.6f", b.get(i));
		for(int j = 0; j < c.length; j++){
			for(int i = 0; i < c[j].size(); i++)s+=String.format("\t%7.6f", c[j].get(i));
		}
		
		return s;
	}
	
	/**
	 * @param c0 origin
	 * @param c destination
	 * @param length
	 */
	private void copiarVector(DenseVector[] c0, DenseVector[] c, int length) {
		for(int k=0; k < length; k++){
			c[k]=c0[k].copy();
			
		}
	}
	
}
