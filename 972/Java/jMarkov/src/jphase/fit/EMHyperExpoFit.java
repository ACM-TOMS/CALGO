package jphase.fit;


import jphase.DenseContPhaseVar;


/**
 * This class implements the Maximum Likelihood method 
 * proposed by Khayari, Sadre and Haverkort in "Fitting
 * world-wide web request traces with the EM algorithm", 2003.
 * The method returns a Hyper-Exponential distribution, a 
 * subclass of Phase-type distributions.
 * @author Juan F. Pérez
 * @version 1.0 
 */
public class EMHyperExpoFit extends MLContPhaseFitter  {
	
	
	/**
	 * @see jphase.fit.MLContPhaseFitter#MLContPhaseFitter(double[])
	 */
	public EMHyperExpoFit(double[] data){
		super(data);
	}
	
	
	
	/**
	 * Precision for the convergence criterion
	 * in the algorithm
	 */
	public static double precision = 10E-5;
	
	/**
	 * Precision for the convergence criterion in
	 * the coefficient of variance
	 */
	public static double precisionParam = 10E-3;
	
	/**
	 * Number of exponential phases in the hyper-exponential
	 * distribution in evaluation 
	 */
	private int n;
	
	
	/**
	 * Probability of choose each exponential phase in the 
	 * hyper-exponential distribution in evaluation 
	 */
	private double[] probs;
	
	/**
	 * Rate at each exponential phase in the 
	 * hyper-exponential distribution in evaluation 
	 */
	private double[] rates;
	
	
	/**
	 * 
	 * @return HyperExponential variable with the best fit
	 */
	@Override
	public DenseContPhaseVar fit() {
		int N = 15;
		int best = 0;
		double[] LH = new double[N]; 
		DenseContPhaseVar[] vars = new DenseContPhaseVar[N];
		for(int k=0; k < N; k++){
			System.out.println("\n\nITERATION: "+(k+1));
			System.out.println("N: "+(k+1));
			this.n=k+1;
			vars[k] = new DenseContPhaseVar(k+1);
			LH[k] = doFitN(data);
			vars[k] = DenseContPhaseVar.HyperExpo(rates, probs);
			System.out.println("Likelihood: "+LH[k]);
			System.out.println("Mean: "+vars[k].expectedValue());
			System.out.println("Variance: "+vars[k].variance());
			System.out.println("CV: "+vars[k].CV());
			if(LH[k]>LH[best])best=k;
		}
		System.out.println("The best variable found is: "
				+vars[best].toString());
		System.out.println("Likelihood: "+LH[best]);
		return vars[best];
	}
	
	
	
	/**
	 * This method implements the EM algorithm from the data, with 
	 * the specified number of exponential phases
	 * @param data non-negative data trace from independent experiments
	 * @return The loglikelihood of the solution found. The values of the
	 * parameters are stored in the attributes probs and rates
	 */
	public double doFitN(double[] data){
		//Initialization
		probs = new double[n];
		rates = new double[n];
		for(int i = 0; i < this.n; i++){
			probs[i] = 1.0/this.n;
			probs[i] = (0.9 * Math.pow(10, -i ))/(1-Math.pow(10, -this.n));
			rates[i] = Math.pow(10, -i+1);
		}

		//EM algorithm
		int iter = 0;
		int maxIter = 300;
		boolean ready = false;
		double[][] p = new double[data.length][this.n];
		double LHold = 0;
		double LHnew = 0;
		while(ready == false){
			iter ++;
			eStep(data, p);
			LHold = LHnew;
			LHnew = mStep(data, p);
			double dif = Math.abs(LHold-LHnew);
			if(dif < precision || iter > maxIter )ready = true;
		}
		return LHnew;
	}

	
	/**
	 * Executes the Expectation step of the EM algorithm
	 * @param data non-negative data trace from independent experiments
	 * @param p likelihood estimate to be computed
	 */
	private void eStep(double[] data, double[][]p){
		for(int n = 0; n < data.length; n++){
			double denom =0.0;
			for(int i = 0; i < this.n; i++){
				 double num = probs[i]*rates[i]*Math.exp(-rates[i]*data[n]);
				 p[n][i]= num;
				 denom += num;
			}
			for(int i = 0; i < this.n; i++)p[n][i] /= denom;
		}

	}
	
	/**
	 * Executes the Maximization step of the EM algorithm
	 * @param data non-negative data trace from independent experiments
	 * @param p likelihood estimate calculated at the Expectation step
	 * @return log-likelihood of the resulting parameter set
	 */
	private double mStep(double[] data, double[][] p){
		int N = data.length;
		for(int i = 0; i < this.n; i++){
			double sum = 0;
			double sumX = 0;
			for(int n = 0; n < N; n++){
				sum += p[n][i];
				sumX += p[n][i]*data[n];
			}
			probs[i] = sum/N;
			rates[i] = sum/sumX;
		}
		double LH=0; 
		for(int n = 0; n < N; n++){
			double sumLH = 0;
			for(int i = 0; i < this.n; i++){
				sumLH += probs[i]*rates[i]*Math.exp(-rates[i]*data[n]);
			}
			LH += Math.log(sumLH); 
		}
			return LH;
	}
}
