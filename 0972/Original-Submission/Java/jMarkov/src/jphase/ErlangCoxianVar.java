package jphase;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;

/**
 * Phase-Type representation of an ErlangCoxian 
 * distribution as defined by Osogami and Harchol in "Closed 
 * form solutions for mapping general distributions to 
 * quasi-minimal PH distributions", 2005.
 * @author Juan F. Pérez
 * @version 1.0
 */
public class ErlangCoxianVar  extends AbstractContPhaseVar implements PhaseVar {

     /** Total number of phases (Erlang degree: n-2)
      */
     private int n;
     
     /** Probability of having a positive elapse time
      * in the distribution. 1-p: mass probability at zero
      */
     private double p;
     
     /** Erlang distribution rate*/
     private double lambdaY;
     
     /** Rate of the first stage at the Coxian 
     * distribution*/
     private double lambdaX1;
     
     /** Rate of the second stage af the Coxian 
     * distribution*/
     private double lambdaX2; 
     
     /** Probability of going from the first to the second
     * stage in the Coxian distribution. 1-p: probability of 
     * absorption at the first stage if the Coxian distribution*/
     private double px;
     
        
    /**
     * Constructor of an Erlang Coxian variable in dense representation.
     * As default it has just one phase in the Erlang that is taken with probability
     * one. This phase and those of the Coxian distribution has rate 1 per time unit. 
     */
    public ErlangCoxianVar() {
        this.n=3;
        this.p=1;
        this.lambdaY=1;
        this.lambdaX1=1;
        this.lambdaX2=1;
        this.px=1;
    }
    
    /**
     * Constructor of a ErlangCoxian variable in dense representation
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
     */
    public ErlangCoxianVar(int n, double p, double lambdaY, double lambdaX1, double lambdaX2, double px){
        if(n < 1)
            throw new IllegalArgumentException("The number of phases n must be >= 1");
        else if(p < 0.0 || p > 1.0)
            throw new IllegalArgumentException("The probability p must be 0 <= p <= 1");
        else if(lambdaY < 0)
            throw new IllegalArgumentException("The rate lambdaY must be >= 0");
        else if(lambdaX1 < 0)
            throw new IllegalArgumentException("The rate lambdaX1 must be >= 0");
        else if(lambdaX2 < 0)
            throw new IllegalArgumentException("The rate lambdaX2 must be >= 0");
        else if(px < 0.0 || px > 1.0)
            throw new IllegalArgumentException("The probability px must be 0 <= px <= 1");
        else{
            this.n=n;
            this.p=p;
            this.lambdaY=lambdaY;
            this.lambdaX1=lambdaX1;
            this.lambdaX2=lambdaX2;
            this.px=px;
        }
    }

    /**
     * Constructor of a ErlangCoxian variable in dense representation
     * @param n total number of phases (Erlang degree: n-2) 
     * absorption at the first stage of the Coxian distribution 
     */
    public ErlangCoxianVar(int n){
        if(n < 1)
            throw new IllegalArgumentException("The number of phases n must be >= 1");
        else{
            this.n=n;
        }
    }
    
    /**
     * @return Total number of phases
     */
    public int getN(){
        return this.n;
    }
    
    /**
     * @param n Total number of phases
     */
    public void setN(int n){
        this.n=n;
    }
     
    /**
     * @return Probability of having a positive elapse time
     * in the distribution. 1-p: mass probability at zero
     */
    public double getP(){
        return this.p;
    }
    
    /**
     * @param p Probability of having a positive elapse time
      * in the distribution. 1-p: mass probability at zero
     */
    public void setP(double p){
        this.p = p;
    }
    
    /**
     * @return Erlang distribution rate
     */
    public double getLambdaY(){
        return this.lambdaY;
    }
    
    /**
     * @param lambdaY Erlang distribution rate
     */
    public void setLambdaY(double lambdaY){
        this.lambdaY = lambdaY;
    }

    /**
     * @return Rate of the first stage at the Coxian 
     * distribution*/
    public double getLambdaX1(){
        return this.lambdaX1;
    }
    
    /**
     * @param lambdaX1 Rate of the first stage at the Coxian 
     * distribution*/
    public void setLambdaX1(double lambdaX1){
        this.lambdaX1=lambdaX1;
    }
    
    /**
     * @return Rate of the second stage at the Coxian 
     * distribution*/
    public double getLambdaX2(){
        return this.lambdaX2;
    }
    
    /**
     * @param lambdaX2 Rate of the second stage at the Coxian 
     * distribution*/
    public void setLambdaX2(double lambdaX2){
        this.lambdaX2=lambdaX2;
    }
    
    /**
     * @return Probability of going from the first to the second
     * stage in the Coxian distribution. 1-p: probability of 
     * absorption at the first stage if the Coxian distribution
     */
    public double getPx(){
        return this.px;
    }
    
    /**
     * @param px Probability of going from the first to the second
     * stage in the Coxian distribution. 1-p: probability of 
     * absorption at the first stage if the Coxian distribution
     */
    public void setPx(double px){
        this.px=px;
    }
    
        
    public Matrix getMatrix() {
        double[][] matriz = new double[n][n];
        if(n==1)matriz[0][0] = -lambdaX2;
        else if(n==2){
            matriz[0][0]= - lambdaX1;
            matriz[0][1]= px*lambdaX1;
            matriz[1][1]= - lambdaX2;
        }else{
            for (int i = 0; i < n - 2; i++) {
                matriz[i][i] = -lambdaY;
                matriz[i][i + 1] = lambdaY;
            }
            matriz[n-2][n-2] = -lambdaX1;
            matriz[n-2][n-1] = px*lambdaX1;
            matriz[n-1][n-1] = -lambdaX2;
        }

        return new DenseMatrix(matriz);
    }

    public void setMatrix(Matrix A) {
        // TODO Auto-generated method stub
    }

    public Vector getVector() {
        double[] vector = new double[n];
        vector[0]=p;
        return new DenseVector(vector);
    }

    public void setVector(Vector alpha) {
        // TODO Auto-generated method stub
    }

    
    public ContPhaseVar copy() {
        return new ErlangCoxianVar(
                this.n,  this.p, this.lambdaY,
                this.lambdaX1, this.lambdaX2, this.px
        );
    }
    
    public ContPhaseVar newVar(int n){
        return new ErlangCoxianVar(n);
    }
    
    /**
     * Creates a Dense Continuous Phase Variable that represents
     * the mixture of the original ErlangCoxian distribution (p) 
     * and an exponential phase with rate lambda (1-p) 
     * @param lambda rate of the exponential phase to be included
     * in the mixture
     * @param p probability mass of the ErlangCoxian distribution
     * in the mixture
     * @return Dense Continuous Phase Variable that represents
     * the mixture of the original ErlangCoxian distribution (p) 
     * and an exponential phase with rate lambda (1-p)
     */
    public PhaseVar mixtureExpo(double lambda, double p){
        Matrix oldMat = this.getMatrix();
        Vector oldVec = this.getVector();
        Matrix newMat = new DenseMatrix(n+1, n+1);
        Vector newVec = new DenseVector(n+1);
        double[][] matrizExpo = {{-lambda}};
        newMat = MatrixUtils.concatQuad(oldMat, new DenseMatrix(n,1), 
                new DenseMatrix(1,n), new DenseMatrix(matrizExpo),
                newMat);
        newVec.set(0  , oldVec.get(0)*p );
        newVec.set(n, 1-p );
        return new DenseContPhaseVar(newVec, newMat);
    }
    
    
    /**
     * Creates a Dense Continuous Phase Variable that represents
     * the convolution of the original ErlangCoxian distribution 
     * and an exponential phase with rate lambda  
     * @param lambda rate of the exponential phase to be included
     * in the convolution
     * @return Dense Continuous Phase Variable that represents
     * the convolution of the original ErlangCoxian distribution 
     * and an exponential phase with rate lambda
     */
    public PhaseVar convoExpo(double lambda){
        Matrix oldMat = this.getMatrix();
        Vector oldVec = this.getVector();
        Matrix newMat = new DenseMatrix(n+1, n+1);
        Vector newVec = new DenseVector(n+1);
        double[][] matrizRightUp = new double[n][1];
        matrizRightUp[n-1][0]=-lambdaX2;
        double[][] matrizExpo = {{-lambda}};
        newMat = MatrixUtils.concatQuad(oldMat, new DenseMatrix(matrizRightUp), 
                new DenseMatrix(1,n), new DenseMatrix(matrizExpo),
                newMat);
        newVec.set(0  , oldVec.get(0) );
        return new DenseContPhaseVar(newVec, newMat);
    }

    @Override
    public String description(){
        String s = "__________________________________________________\n";
        s += "ErlangCoxian Distribution";
        s+= "\nNumber of Phases: "+this.n+"\n";
        s += "Probability of Positive time : "+this.p+"\n";
        s += "Erlang rate : "+this.lambdaY+"\n";
        s += "Coxian First stage rate : "+this.lambdaX1+"\n";
        s += "Coxian Second stage rate : "+this.lambdaX2+"\n";
        s += "NON-absorption Probability at Coxian First Stage : "+this.px+"\n";
        s+="__________________________________________________\n";
        return s;
    }
}
