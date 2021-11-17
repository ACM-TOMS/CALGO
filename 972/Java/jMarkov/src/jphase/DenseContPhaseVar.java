package jphase;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;

/** 
 * This class allows the creation and manipulation 
 * of Continuous Phase-type distributions represented by dense matrices.
 * @author German Riaño
 * @author Juan F. Perez
 * @version 1.0 
 */
public class DenseContPhaseVar extends AbstractContPhaseVar {
	
	/**
	 *  Rate Matrix
	 */
	protected DenseMatrix A;

	/**
	 *  Initial probability distribution vector
	 */
	protected DenseVector alpha;
	
	/**
	 * Constructs an empty continuous Phase-type Distribution 
	 * with dense representation 
	 */
    public DenseContPhaseVar() {
        this.A = new DenseMatrix(new double[][] { { -1.0 } });
        this.alpha = new DenseVector(new double[]  { 0 } );
    }
	
    /**
	 * Constructs an empty Continuous Phase-type Distribution 
	 * of size n with dense representation
	 * @param n size of the Continuous Phase-type Distribution  
	 */
    public DenseContPhaseVar(int n) {
        this.A = new DenseMatrix(n,n);
        this.alpha = new DenseVector(n);
    }
    
	/**
	 * Constructs a continuous Phase-type Distribution 
	 * with dense representation
	 * @param A rate matrix
	 * @param alpha initial probability distribution vector
	 */
	public DenseContPhaseVar(DenseVector alpha, DenseMatrix A){
		try{
        	if( checkContPhaseVar(alpha, A)){
        		this.A = A;
        		this.alpha = alpha;
        	}
        }catch (Exception e){
        	e.printStackTrace(); 
        }
		
	}

	/**
	 * Constructs a continuous Phase-type Distribution 
	 * with dense representation
	 * @param A rate matrix
	 * @param alpha initial probability distribution vector
	 */
	public DenseContPhaseVar(Vector alpha, Matrix A){
		try{
        	if( checkContPhaseVar(alpha, A)){
        		this.A = (DenseMatrix)A.copy();
        		this.alpha = (DenseVector)alpha.copy();
        	}
        }catch (Exception e){
        	e.printStackTrace(); 
        }
		
	}
	
	/**
	 * Constructs a continuous Phase-type Distribution 
	 * with dense representation
	 * @param A rate matrix
	 * @param alpha initial probability distribution vector
	 * @throws Exception 
	 */
    public DenseContPhaseVar(double[] alpha, double[][] A) {
        int n1, n2;
        n1 = alpha.length;
        n2 = A.length;
        if (n1 != n2) {
            throw new IndexOutOfBoundsException(
                    "Vector and matrix not of the same size");
        } else {
            DenseMatrix Am= new DenseMatrix(A);
            DenseVector Vm= new DenseVector(alpha);
            try{
            	if( checkContPhaseVar(Vm, Am)){
	            	this.A= Am;
	            	this.alpha = Vm;
            	}
            }catch (Exception e){
            	e.printStackTrace();
            }
        }
    }
    
    public boolean checkContPhaseVar(Vector alpha, Matrix A) throws Exception{
    	if (!MatrixUtils.checkSubStochasticVector((Vector)alpha)){
    		throw new Exception(
                    "Vector provided is not sub-stochastic");
    	}else if (!MatrixUtils.checkSubGeneratorMatrix(A)){
    		throw new Exception(
                    "Matrix provided is not a sub-generator");
        	
        }
        return true;
    }
    
    
	/**
	 * Constructs a Phase-Type representation of an Exponential 
	 * distribution with rate lambda
	 * @param lambda exponential distribution rate
	 * @return Dense Continuous Phase-Type Distribution
	 */
	public static DenseContPhaseVar expo(double lambda){
		double[][] matriz = new double[1][1];
		matriz[0][0] = -lambda;
		DenseMatrix A = new DenseMatrix(matriz);

		double[] vector = new double[1];
		vector[0] = 1;
		DenseVector alpha = new DenseVector(vector);
		try{
			return new DenseContPhaseVar(alpha, A);
		}catch(Exception e){
			e.printStackTrace();
			return new DenseContPhaseVar();
		}
	}

	/**
	 * Constructs a Phase-Type representation of an Erlang 
	 * distribution with rate lambda and n exponential phases
	 * @param lambda exponential rate in each phase
	 * @param n number of exponential phases
	 * @return Dense Continuous Phase-Type Distribution
	 */
	public static DenseContPhaseVar Erlang(double lambda, int n){
		double[][] matriz = new double[n][n];
		for (int i = 0; i < n; i++) {
			matriz[i][i] = -lambda;
			if (i < n - 1)
				matriz[i][i + 1] = lambda;
		}

		DenseMatrix A = new DenseMatrix(matriz);

		double[] vector = new double[n];
		vector[0] = 1;
		DenseVector alpha = new DenseVector(vector);
		try{
			return new DenseContPhaseVar(alpha, A);
		}catch(Exception e){
			e.printStackTrace();
			return new DenseContPhaseVar();
		}
	}
	
	public static DenseContPhaseVar HipoExponential(double[] lambdas){
		HypoExponentialVar var = new HypoExponentialVar(lambdas);
		Matrix A = var.getMatrix();
		Vector alpha = var.getVector();
		try{
			return new DenseContPhaseVar(alpha, A);
		}catch(Exception e){
			e.printStackTrace();
			return new DenseContPhaseVar();
		}
	}

	
	/**
	 * Constructs a Phase Distribution that represents
	 * a HyperExponential distribution with the especified parameters 
	 * @param lambdas ecah one of the exponential rates
	 * @param probs initial probability vector
	 * @return Dense Continuous Phase-Type Distribution
	 */
	public static DenseContPhaseVar HyperExpo(double[] lambdas, double[] probs){
		if(lambdas.length==probs.length){
			int n=lambdas.length;
			double[][] matrix = new double[n][n];
			double[] vector = new double[n];
			for (int i = 0; i < n; i++) {
				matrix[i][i] = -lambdas[i];
				vector[i] = probs[i];
			}
			DenseMatrix A = new DenseMatrix(matrix);
			DenseVector alpha = new DenseVector(vector);
			try{
				return new DenseContPhaseVar(alpha, A);
			}catch(Exception e){
				e.printStackTrace();
				return new DenseContPhaseVar();
			}
		}
        throw new RuntimeException("Rates and probability vectors have different length");
	}
	
	/**
	 * Constructs a Phase-Type representation of a Hyper-Erlang 
	 * distribution with k erlang branches, its k rates and n number
	 * of phases per branch
	 * @param k number of Erlang branches
	 * @param lambdas exponential rate of each phase in each branch
	 * @param n number of exponential phases in each branch
	 * @param probs probability distribution of taking each one of the branches 
	 * @return Dense Continuous Phase-Type Distribution
	 */
	public static DenseContPhaseVar HyperErlang(
		int k,
		double[] lambdas,
		int[] n,
		double[] probs){
		int N, l, r;
		N = l = r = 0;
		for (int i = 0; i < k; i++)
			N += n[i];
		double[][] matriz = new double[N][N];
		double[] vector = new double[N];

		for (int i = 0; i < k; i++) {
			vector[l] = probs[i];
			for (int j = 0; j < n[i]; j++) {
				matriz[l][l] = -lambdas[i];
				if (l < N - 1 ) {
					if (j < n[i] -1) {
					matriz[l][l + 1] = lambdas[i];
					}
				}
				l++;
			}
			r+=n[i];
		}

		DenseMatrix A = new DenseMatrix(matriz);
		DenseVector alpha = new DenseVector(vector);
		try{
			return new DenseContPhaseVar(alpha, A);
		}catch(Exception e){
			e.printStackTrace();
			return new DenseContPhaseVar();
		}
	}

	
	/**
	 * Constructs a Phase-Type representation of a Hyper-Erlang 
	 * distribution from a Dense representation of the same 
	 * distribution
	 * @param var HyperErlang variable from which the Dense 
	 * Continuous variable must be constructed
	 * @return Dense Continuous Phase-Type Distribution
	 */
	public static DenseContPhaseVar HyperErlang(HyperErlangVar var){
		DenseMatrix A = new DenseMatrix(var.getMatrix());
		DenseVector alpha = new DenseVector(var.getVector());
		try{
			return new DenseContPhaseVar(alpha, A);
		}catch(Exception e){
			e.printStackTrace();
			return new DenseContPhaseVar();
		}
	}
	
	
	/**
	 * Constructs a Phase-Type representation of a Coxian 
	 * distribution with n phases
	 * @param n number of phases
	 * @param lambdas exponential rates of each phase
	 * @param probs probability of going to the next phase 
	 * (no absorption) in each phase except the last one.  
	 * @return Dense Continuous Phase-Type Distribution
	 */
	public static DenseContPhaseVar Coxian(int n, double[] lambdas, double[] probs){
		if(n <= 0)throw new IllegalArgumentException("The number of phases in the" +
		"Coxian distribution must be greater than zero");
		
		if(lambdas.length != n) throw new IllegalArgumentException("" +
				"The number of rates is different from the number of phases");
		
		if(probs.length != n-1) throw new IllegalArgumentException("" +
		"The number of probabilities is different from the number of phases - 1");
		
		double[][] matriz = new double[n][n];
		double[] vector = new double[n];
		vector[0]=1.0;
		
		for (int i = 0; i < n; i++) {
			matriz[i][i] = -lambdas[i];
			if(i < n-1)matriz[i][i+1]=probs[i]*lambdas[i];
				
		}

		DenseMatrix A = new DenseMatrix(matriz);
		DenseVector alpha = new DenseVector(vector);
		try{
			return new DenseContPhaseVar(alpha, A);
		}catch(Exception e){
			e.printStackTrace();
			return new DenseContPhaseVar();
		}
	}
	
	
	/**
	 * Constructs a Phase-Type representation of an ErlangCoxian 
	 * distribution as defined by Osogami and Harchol in "Closed 
	 * form solutions for mapping general distributions to 
	 * quasi-minimal PH distributions", 2005.
	 * @param n total number of phases (Erlang degree: n-2)
	 * @param p probability of having a positive elapse time
	 * in the distribution. 1-p: mass probability at zero
	 * @param lambdaY rate of the Erlang distribution
	 * @param lambdaX1 rate of the first stage of the Coxian 
	 * distribution
	 * @param lambdaX2 rate of the second stage of the Coxian 
	 * distribution 
	 * @param px probability of going from the first to the second
	 * stage in the Coxian distribution. 1-p: probability of 
	 * absorption at the first stage of the Coxian distribution 
	 * @return Dense Continuous Phase-Type Distribution
	 */
	public static DenseContPhaseVar ErlangCoxian(
		int n,
		double p,
		double lambdaY,
		double lambdaX1,
		double lambdaX2,
		double px) {
		
		if(n <= 0)throw new IllegalArgumentException("The number of phases in the" +
				"Erlang-Coxian distribution must be greater than zero");
		else if(n == 1)return expo(lambdaX1);
		else if(n == 2)return Coxian(n, new double[]{lambdaX1, lambdaX2}, new double[]{px});
		else{
			double[][] matriz = new double[n][n];

			
			for (int i = 0; i < n - 2; i++) {
				matriz[i][i] = -lambdaY;
				matriz[i][i + 1] = lambdaY;
			}
			matriz[n-2][n-2] = -lambdaX1;
			matriz[n-2][n-1] = px*lambdaX1;
			matriz[n-1][n-1] = -lambdaX2;
			
			double[] vector = new double[n];
			vector[0]=p;
						
			DenseMatrix A = new DenseMatrix(matriz);
			DenseVector alpha = new DenseVector(vector);
			return new DenseContPhaseVar(alpha, A);
		}
	}
	
	public Matrix getMatrix() {
		return this.A;
	}

	public Vector getVector() {
		return this.alpha;
	}

	public void setMatrix(Matrix A) {
		this.A = (DenseMatrix) A.copy();
	}

	public void setVector(Vector alpha) {
		this.alpha = (DenseVector) alpha.copy();
	}
	
	public ContPhaseVar copy(){
		return new DenseContPhaseVar(
				this.getVector(), this.getMatrix());
	}
	
	public ContPhaseVar newVar(int n){
        return new DenseContPhaseVar(n);
    }

}
