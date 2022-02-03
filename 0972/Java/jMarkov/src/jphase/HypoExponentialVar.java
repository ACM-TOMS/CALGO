package jphase;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;


/**
 * This class allows the creation and manipulation 
 * of hypoexponential distributions. The associated matrix has a   
 * dense representation 
 * @author German Riaño
 * @author Juan F. Pérez
 * @version 1.0
 * 
 */
public class HypoExponentialVar  extends AbstractContPhaseVar implements PhaseVar {
    
    /**
     * Total number of phases in the hypoexponential distribution 
     */
    private int N;
    
    /**
     * Rates of each Erlang branch
     */
    private double[] lambdas;
    
    
    /**
     * Constructor of a hypoexponential variable in dense representation.
     * As default it has just one branch that is taken with probability
     * one. the unique branch has one phase with rate 1 per time unit. 
     */
    public HypoExponentialVar() {
        this.N=1;
        this.lambdas=new double[1];
        this.lambdas[0]=1;
    }
    
    /**
     * Constructor of a hypoexponential variable with n phases
     * in dense representation
     * @param n Total number of phases
     */
    public HypoExponentialVar(int n){
        this.N=n;
    }  
    
    
    /**
     * Constructor of a hypoexponential variable in dense representation
     * @param lambdas Rate associated to each branch
     */
    public HypoExponentialVar(double[] lambdas) {
		N = lambdas.length;
		this.lambdas = lambdas;
    }
    
    /**
     * Return the rates associated to each branch
     * @return Rate associated to each branch
     */
    public double[] getLambdas(){
        return this.lambdas;
    }
    
    /**
     * Sets the rates associated to each branch
     * @param lambdas Rates associated to each branch to set
     */
    public void setLambdas(double[] lambdas){
        this.lambdas=lambdas;
    }
    /**
     * 
     */
	public Matrix getMatrix() {
		return new DenseMatrix(getDMatrix());        
    }
    
    public double[][] getDMatrix(){
		double[][] matriz = new double[N][N];
		for (int i = 0; i < N; i++) {
			double lambda = lambdas[i];
			matriz[i][i] = -lambda;
			if (i < N - 1)
				matriz[i][i + 1] = lambda;
		}

        return matriz;
    }

	public void setMatrix(Matrix A) {
        throw new UnsupportedOperationException();
    }

	public Vector getVector() {
        return new DenseVector(getDVector());
    }
    
    public double[] getDVector(){
		double[] vector = new double[N];
		vector[0] = 1;
        return vector;
    }

	public void setVector(Vector alpha) {
        throw new UnsupportedOperationException();
    }

    public ContPhaseVar copy() {
        return new HypoExponentialVar(this.lambdas);
    }
    
	public ContPhaseVar newVar(int n){
        return new HypoExponentialVar(n);
    }
    
    @Override
    public int getNumPhases(){
        return N;
    }
   
        
    @Override
    public double expectedValue() {
    	double suma = 0;
    	for(int i = 0; i!= lambdas.length;i++)
    		suma+=(1/lambdas[i]);
    	return suma;
    }
}
