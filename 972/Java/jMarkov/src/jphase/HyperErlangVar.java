package jphase;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;


/**
 * This class allows the creation and manipulation 
 * of HyperErlang distributions. The associated matrix has  
 * dense representation 
 * @author German Riaño
 * @author Juan F. Perez
 * @version 1.0
 */
public class HyperErlangVar  extends AbstractContPhaseVar implements PhaseVar {
    
    /**
     * Total Number of phases in the HyperErlang distribution 
     */
    private int N;
    
    /**
     * Number of Erlang branches 
     */
    private int M;
    
    /**
     * Number of phases on each Erlang branch
     */
    private int[] r;
    
    /**
     * Probability mass function of choosing each one of the branches 
     */
    private double[] alphas;
    
    /**
     * Rates of each erlang branch
     */
    private double[] lambdas;
    
    
    /**
     * Constructor of a Hyper Erlang variable in dense representation.
     * As default it has just one branch that is taken with probability
     * one. the unique branch has one phase with rate 1 per time unit. 
     */
    public HyperErlangVar() {
        this.N=1;
        this.M=1;
        this.r=new int[M];
        this.r[0]=1;
        this.alphas=new double[M];
        this.alphas[0]=1;
        this.lambdas=new double[M];
        this.lambdas[0]=1;
    }
    
    /**
     * Constructor of a Hyper Erlang variable with n phases
     * in dense representation
     * @param n Total number of phases
     */
    public HyperErlangVar(int n){
        this.N=n;
    }
    
    /**
     * Constructor of a Hyper Erlang variable in dense representation
     * @param N Total number of phases
     * @param M Number of branches
     * @param r Number of phases in each branch
     * @param alphas Probability associated to each branch 
     * @param lambdas Rate associated to each branch
     * @param deep True if this is a deep constructor, false if not
     */
    public HyperErlangVar(int N, int M, int[] r, double[] alphas, double[] lambdas, boolean deep){
        if(M==r.length && M==alphas.length && M==lambdas.length){
            this.N=N;
            this.M=M;
            if(deep==true){
                this.r=new int[M];
                this.alphas=new double[M];
                this.lambdas=new double[M];
                System.arraycopy(r,0,this.r,0,M);
                System.arraycopy(alphas,0,this.alphas,0,M);
                System.arraycopy(lambdas,0,this.lambdas,0,M);
            }else{
                this.r=r;
                this.alphas=alphas;
                this.lambdas=lambdas;
            }
        }
    }
    
    
    
    /**
     * Constructor of a Hyper Erlang variable in dense representation
     * @param r Number of phases in each branch
     * @param alphas Probability associated to each branch 
     * @param lambdas Rate associated to each branch
     * @param deep True if this is a deep constructor, false if not
     */
    public HyperErlangVar(int[] r, double[] alphas, double[] lambdas, boolean deep) {
        if(r.length==alphas.length && alphas.length==lambdas.length){
            this.M=r.length;
            for(int i=0;i<M;i++){
                this.N+=r[i];
            }
            if(deep==true){
                this.r=new int[M];
                this.alphas=new double[M];
                this.lambdas=new double[M];
                System.arraycopy(r,0,this.r,0,M);
                System.arraycopy(alphas,0,this.alphas,0,M);
                System.arraycopy(lambdas,0,this.lambdas,0,M);
            }else{
                this.r=r;
                this.alphas=alphas;
                this.lambdas=lambdas;
            }
        }
    }
    
    
    /**
     * @return Total number of phases
     */
    public int getN(){
        return this.N;
    }
    
    /**
     * @param N Total number of phases to set
     */
    public void setN(int N){
        this.N=N;
    }
    
    /**
     * @return Nnumber of branches
     */
    public int getM(){
        return this.M;
    }
    
    /**
     * @param M Number of branches to set
     */
    public void setM(int M){
        this.M=M;
    }
    
    
    /**
     * @return Number of phases in each branch
     */
    public int[] getR(){
        return this.r;
    }
    
    /**
     * @param r Number of phases in each branch to set
     */
    public void setR(int[] r){
        this.r=r;
    }
    

    /**
     * @return Probability associated to each branch
     */
    public double[] getAlphas(){
        return this.alphas;
    }
    
    /**
     * @param alphas Probability associated to each branch to set
     */
    public void setAlphas(double[] alphas){
        this.alphas=alphas;
    }
    
    
    /**
     * @return Rate associated to each branch
     */
    public double[] getLambdas(){
        return this.lambdas;
    }
    
    /**
     * @param lambdas Rates associated to each branch to set
     */
    public void setLambdas(double[] lambdas){
        this.lambdas=lambdas;
    }

	public Matrix getMatrix() {
        return new DenseMatrix(getDMatrix());        
    }
    
    /**
     * Returns the Double MAtrix that represents the variable
     * @return Double Matrix
     */
    public double[][] getDMatrix(){
        int N = this.getN();
        int l =0;
        double[][] matriz = new double[N][N];
        for (int i = 0; i < this.getM(); i++) {
            for (int j = 0; j < this.getR()[i]; j++) {
                matriz[l][l] = -this.getLambdas()[i];
                if (l < N - 1 && j < this.getR()[i]-1) {
                    matriz[l][l + 1] = this.getLambdas()[i];
                }
                l++;
            }
        }

        return matriz;
    }

	public void setMatrix(Matrix A) {
        throw new UnsupportedOperationException();
    }

	public Vector getVector() {
        return new DenseVector(getDVector());
    }

    
    /**
     * Returns the Double MAtrix that represents the variable
     * @return Double Vector
     */
    public double[] getDVector(){
        int N = this.getN();
        int l =0;
        
        double[] vector = new double[N];
//        System.arraycopy(this.getAlphas(),0,vector,0,this.getM());
        for (int i = 0; i < this.getM(); i++) {
            for (int j = 0; j < this.getR()[i]; j++) {
            	if(j == 0){
                    vector[l] = this.getAlphas()[i];	
            	}
            	else 
            		vector[l] = 0;
                l++;
            }
        }
        return vector;

//        double[] vector = new double[N];
//        System.arraycopy(this.getAlphas(),0,vector,0,this.getM());
//        return vector
    }

	public void setVector(Vector alpha) {
        throw new UnsupportedOperationException();
    }

    
	public ContPhaseVar copy() {
        return new HyperErlangVar(
                this.N, this.M, this.r,
                this.alphas, this.lambdas, true
        );
    }
    
	public ContPhaseVar newVar(int n){
        return new HyperErlangVar(n);
    }
    
    @Override
    public int getNumPhases(){
        return this.getN();
    }
   
        
    @Override
    public double expectedValue() {
        double m = 0;
        for(int i=0;i<this.getM();i++)m+=this.getAlphas()[i]*this.getR()[i]/this.getLambdas()[i];
        return m;
    }

    
    @Override
    public double moment(int k) {
        double sum = 0;
        for(int i = 0; i< this.M; i++ ){
            double temp1 =this.alphas[i]/(Math.pow(this.lambdas[i],k));
            for(int j = 0; j <=k-1;j++)temp1*=(this.r[i]+j);
            sum+=temp1;
        }
        return sum;
    }
    
    /*
    @Override
    public String description(){
        String s = "__________________________________________________\n";
        s += "HyperErlang Distribution";
        s+= "\nNumber of Branches: "+this.M+"\n";
        s += "Number of Phases: "+this.N+"\n";
        s=s+"\tBranch\t Ri\t alphai\t\t lambdai\n";
        //for(int i = 0 ; i < this.M; i++)s+=""+(i+1)+"\t "+this.r[i]+"\t"+this.alphas[i]+"\t"+this.lambdas[i]+"\n";
        for(int i = 0 ; i < this.M; i++){
            s+="\t"+(i+1)+"\t "+this.r[i];
            s+=String.format("\t %5.4f",this.alphas[i]);
            s+=String.format("\t\t %6.4f",this.lambdas[i]);
            s+="\n";
        }
        s+="__________________________________________________\n";
        return s;
    }
    */
}
