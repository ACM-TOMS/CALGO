package jphase;

import jphase.MatrixUtils;
import static jphase.Utils.*;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Matrices;


/**
 *  Abstract class that defines the behaviour of the Discrete
 *  Phase-Type distributions
 *  @author German Riaño
 *  @author Juan F. Perez
 *  @version 1.0
 *  
 */

public abstract class AbstractDiscPhaseVar implements DiscPhaseVar {

	/**
     * @see jphase.DiscPhaseVar#sumPH(jphase.DiscPhaseVar, jphase.DiscPhaseVar)
     */
	public DiscPhaseVar sumPH(DiscPhaseVar B, DiscPhaseVar res){
		int n1 = this.getNumPhases();
		int n2 = B.getNumPhases();
		int n3 = res.getNumPhases();
		if (n1*n2 != n3)
			throw new IndexOutOfBoundsException("The closure operation cannot be done:\n " +
					"this.Phases * B.Phases != res.Phases ("
					+ n1 + " * " + n2 +" != " + n3 + ")");
		//ISInv = (I-alpha_(m+1)S)^(-1)
		Matrix ISInv = Matrices.identity(n2).add(-this.getVec0(),B.getMatrix());
		ISInv = ISInv.solve(Matrices.identity(this.getNumPhases()), ISInv.copy());
		
		res.setVector( MatrixUtils.kroneckerVectors(this.getVector(), 
							ISInv.transMult(B.getVector(), B.getVector().copy()),
								res.getVector()));
		//L1 = T kron I
		Matrix L1 = MatrixUtils.kronecker(this.getMatrix(), Matrices.identity(n2), res.getMatrix());
		//L2 = (1-alpha_(m+1))t*alpha k ISInv*S
		Matrix L2 = MatrixUtils.kronecker( (MatrixUtils.multVector(
													this.getMat0(), this.getVector(), 
														this.getMatrix().copy()) ).scale(1-this.getVec0()),
											ISInv.mult(B.getMatrix(), B.getMatrix().copy()),
											res.getMatrix()
										); 
		res.setMatrix(L1.add(L2));
		return res;
	}
	
    /**
     * @see jphase.DiscPhaseVar#sumPH(jphase.DiscPhaseVar)
     */
    public DiscPhaseVar sumPH(DiscPhaseVar v2) {
        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();
        DiscPhaseVar res = this.newVar(n1*n2);
        return this.sumPH(v2, res);
    }
	
	/**
     * @see jphase.DiscPhaseVar#getNumPhases()
     */
	public int getNumPhases() {
		return this.getVector().size();
	}
	
	/**
     * @see jphase.DiscPhaseVar#getVec0()
     */
	public double getVec0(){
    	return 1.0 - this.getVector().dot(MatrixUtils.OnesVector(this.getNumPhases()));
    }
	
    /**
     * @see jphase.DiscPhaseVar#getMat0()
     */
    public Vector getMat0(){
    	return MatrixUtils.OnesVector(this.getNumPhases()).add(
    			-1, this.getMatrix().mult(
    					MatrixUtils.OnesVector(this.getNumPhases()),
    					this.getVector().copy()
	    			)
    			);
    }

	/**
     * @see jphase.DiscPhaseVar#getMatrixArray()
     */
	public double[][] getMatrixArray(){
		 return Matrices.getArray(this.getMatrix());
	}
	
	/**
     * @see jphase.DiscPhaseVar#getVectorArray()
     */
	public double[] getVectorArray(){
		return Matrices.getArray(this.getVector());
	}
	
	/**
     * @see jphase.DiscPhaseVar#getMat0Array()
     */
	public double[] getMat0Array(){
		return Matrices.getArray(this.getMat0());
	}
    
	/**
     * @see jphase.DiscPhaseVar#expectedValue()
     */
	public double expectedValue() {
		return moment(1);
	}
	
	/**
     * @see jphase.DiscPhaseVar#variance()
     */
	public double variance() {
		double m1 = moment(1);
		return moment(2) - (m1 * m1);
	}

	/**
     * @see jphase.DiscPhaseVar#stdDeviation()
     */
	public double stdDeviation() {
		return Math.sqrt(variance());
	}

	/**
     * @see jphase.DiscPhaseVar#CV()
     */
	public double CV() {
		double m=expectedValue();
		return moment(2)/(m*m)-1;
	}
	
	/**
     * @see jphase.DiscPhaseVar#moment(int)
     */
	public double moment(int k) {
		//T=I-A
		Matrix T = Matrices.identity(this.getNumPhases()).add(-1, this.getMatrix());
		//TInv = (I-A)^-1
		Matrix TInv = T.solve(Matrices.identity(this.getNumPhases()), this.getMatrix().copy());
		//TInv = (I-A)^-k
		TInv = MatrixUtils.matPower(TInv, k);
			
		double res =  
				this.getVector().dot(
				 TInv.mult(MatrixUtils.matPower(this.getMatrix(),k-1),this.getMatrix().copy().zero()).mult(
						 MatrixUtils.OnesVector(this.getVector()), this.getVector().copy().zero() 
				 )
				);
        
        res *= fact(k);
        return res;
	}
	
	
    /**
     * @see jphase.DiscPhaseVar#cdf(double)
     */
    public double cdf(double x) {
    	int k = (int) x;
    	if(k==0)return this.getVec0();
    	return 1 - MatrixUtils.matPower(this.getMatrix(), k ,this.getVector(),this.getMat0());
    }
    
    /**
     * @see jphase.DiscPhaseVar#cdf(int, double)
     */
    public double[] cdf(int n, double delta) {
    	double[] res = new double[n+1];
    	res[0] = this.getVec0();
    	for(int i =1; i <= n; i++){
    		res[i] = this.cdf(i*delta);
    	}
    	return res;
    }
    
	
    /**
     * @see jphase.DiscPhaseVar#pmf(int)
     */
	public double pmf(int k) {
		if(k==0)return this.getVec0();
		return MatrixUtils.matPower(this.getMatrix(), k-1 , this.getVector(), this.getMat0());
	}
	
    /**
     * @see jphase.DiscPhaseVar#pmf(int, int)
     */
	public double[] pmf(int n, int delta) {
    	double[] res = new double[n+1];
    	res[0] = this.getVec0();
    	Matrix mat = MatrixUtils.matPower(this.getMatrix(),delta-1);
    	Matrix step = MatrixUtils.matPower(this.getMatrix(),delta);
    	for(int i =1; i <= n; i++){
    		res[i] = this.getVector().dot( 
    				 mat.mult(
    						this.getMat0(), this.getMat0().copy().zero())  );
    		mat = mat.mult(step,mat.copy());
    	}
    	return res;
    }
	
	/**
     * @see jphase.DiscPhaseVar#prob(double, double)
     */
	public double prob(double a, double b) {
        return (b > a) ? cdf(b) - cdf(a) : 0.0;
	}
    
	/**
     * @see jphase.DiscPhaseVar#survival(double)
     */
	public double survival(double x) {
        return (1.0 - this.cdf(x));
	}

	/**
     * @see jphase.DiscPhaseVar#survival(int, double)
     */
	public double[] survival(int n, double delta) {
        double cdfs[] = cdf(n, delta);
        double res[] = new double[n];
        for (int i = 0; i < n; i++)
            res[i] = 1.0 - cdfs[i];
        return res;
	}
    
    /**
     * @see jphase.PhaseVar#lossFunction1(double)
     */
    public double lossFunction1(double x){
        //TODO
        return 0;
    }
    
    /**
     * @see jphase.PhaseVar#lossFunction2(double)
     */
    public double lossFunction2(double x){
        //TODO
        return 0;
    }
	
	/**
     * @see jphase.DiscPhaseVar#quantil(double)
     */
	public double quantil(double p) {
		double x, dif, derv;
		int ITMAX = 100;
		int cnt;
		x = this.expectedValue();
		for (cnt = 0; cnt < ITMAX; cnt++) {
			derv = this.pmf((int)x);
			dif = this.cdf(x) - p;
			if (Math.abs(dif) < 1.0e-10)
				return x;
			x = Math.max(x - dif / derv, 0);
		}
		return 0.0;
	}

	/**
     * @see jphase.DiscPhaseVar#median()
     */
	public double median() {
		return this.quantil(0.5);
	}
	
	/**
     * @see jphase.DiscPhaseVar#sum(jphase.DiscPhaseVar, jphase.DiscPhaseVar)
     */
	public DiscPhaseVar sum(DiscPhaseVar v2, DiscPhaseVar res) {
        if (Matrices.cardinality(this.getVector())==0)
            return v2;
        if (Matrices.cardinality(v2.getVector())==0)
            return this;
        
        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();
        
        Vector A0 = this.getMat0();
        Vector vec1 = this.getVector();
        double vec1_0 = this.getVec0();
        Vector vec2 = v2.getVector();
        
        Vector vec2Temp = vec2.copy();
        vec2Temp.scale(vec1_0);
        
        res.setVector(MatrixUtils.concatVectors(vec1, vec2Temp, res.getVector() ) );
        
        int[] zeroRowsLD = new int[n2];
        int[] zeroColsLD = new int[n1];
        for(int i = 0; i < n2; i++) zeroRowsLD[i] = n1 + i;
        for(int i = 0; i < n1; i++) zeroColsLD[i] = i;

        int[] zeroRowsRU = new int[n1];
        int[] zeroColsRU = new int[n2];
        for(int i = 0; i < n1; i++) zeroRowsRU[i] = i;
        for(int i = 0; i < n2; i++) zeroColsRU[i] = n1 + i;
        
       
        res.setMatrix(MatrixUtils.concatQuad(
				this.getMatrix(),  MatrixUtils.multVector(A0, vec2, Matrices.getSubMatrix(res.getMatrix(), zeroRowsRU , zeroColsRU))	, 
				Matrices.getSubMatrix(res.getMatrix(), zeroRowsLD , zeroColsLD) , v2.getMatrix(), 
				res.getMatrix()
				)
			);
        return res;
	}
	
    /**
     * @see jphase.DiscPhaseVar#sum(jphase.DiscPhaseVar)
     */
    public DiscPhaseVar sum(DiscPhaseVar v2) {
        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();
        DiscPhaseVar res = this.newVar(n1+n2);
        return this.sum(v2, res);
    }
    

    /**
     * @see jphase.DiscPhaseVar#sumGeom(double)
     */
    public DiscPhaseVar sumGeom(double p) {
        DiscPhaseVar res = this.copy();
        res.setVector(this.getVector());
        res.setMatrix(this.getMatrix().copy().add(
                1-p, MatrixUtils.multVector(
                        this.getMat0(),
                        this.getVector(),
                        res.getMatrix().copy().zero()))  );
        return res;
    }
    
	/**
     * @see jphase.DiscPhaseVar#mix(double, jphase.DiscPhaseVar, jphase.DiscPhaseVar)
     */
	public DiscPhaseVar mix(double p, DiscPhaseVar v2, DiscPhaseVar res) {
        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();
        Vector vec1Temp = this.getVector().copy().scale(p);
        Vector vec2Temp = v2.getVector().copy().scale(1-p);
        
        res.setVector( MatrixUtils.concatVectors(vec1Temp, vec2Temp, res.getVector()) );

        int[] zeroRowsLD = new int[n2];
        int[] zeroColsLD = new int[n1];
        for(int i = 0; i < n2; i++) zeroRowsLD[i] = n1 + i;
        for(int i = 0; i < n1; i++) zeroColsLD[i] = i;

        int[] zeroRowsRU = new int[n1];
        int[] zeroColsRU = new int[n2];
        for(int i = 0; i < n1; i++) zeroRowsRU[i] = i;
        for(int i = 0; i < n2; i++) zeroColsRU[i] = n1 + i;
        
        res.setMatrix( MatrixUtils.concatQuad(
        		this.getMatrix(), 	Matrices.getSubMatrix(res.getMatrix(), zeroRowsRU , zeroColsRU), 
        		Matrices.getSubMatrix(res.getMatrix(), zeroRowsLD , zeroColsLD), 	v2.getMatrix(),
        		res.getMatrix())
        		);
        
        return res;
	}
    
    /**
     * @see jphase.DiscPhaseVar#mix(double, jphase.DiscPhaseVar)
     */
    public DiscPhaseVar mix(double p, DiscPhaseVar v2) {
        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();
        DiscPhaseVar res = this.newVar(n1+n2);
        return this.mix(p, v2, res);
    }
    
	
	/**
     * @see jphase.DiscPhaseVar#min(jphase.DiscPhaseVar, jphase.DiscPhaseVar)
     */
	public DiscPhaseVar min(DiscPhaseVar v2, DiscPhaseVar res) {
		int n1 = this.getNumPhases();
		int n2 = v2.getNumPhases();
		int n3 = res.getNumPhases();
		if (n1*n2 != n3)
			throw new IndexOutOfBoundsException("The closure operation cannot be done:\n n1*n2 != n3("
					+ n1 + " * " + n2 + " != " + n3 + ")");
		
		res.setVector(MatrixUtils.kroneckerVectors(this.getVector(), v2.getVector(), res.getVector()));
		res.setMatrix(MatrixUtils.kronecker(this.getMatrix(), v2.getMatrix(), res.getMatrix()));
		return res;
	}
    
    /**
     * @see jphase.DiscPhaseVar#min(jphase.DiscPhaseVar)
     */
    public DiscPhaseVar min(DiscPhaseVar v2) {
        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();
        DiscPhaseVar res = this.newVar(n1*n2);
        return this.min(v2, res);
    }
	
	/**
     * @see jphase.DiscPhaseVar#max(jphase.DiscPhaseVar, jphase.DiscPhaseVar)
     */
	public DiscPhaseVar max(DiscPhaseVar v2, DiscPhaseVar res) {
		int n1 = this.getNumPhases();
		int n2 = v2.getNumPhases();
		int n3 = res.getNumPhases();
		if (n1*n2+n1+n2 != n3)
			throw new IndexOutOfBoundsException("The closure operation cannot be done:\n n1*n2 + n1 + n2 != n3("
					+ n1 + "*" + n2 +" + "+ n1 + " + "+ n2 +" != " + n3 + ")");
       
    
        Vector vec1Temp = this.getVector().copy().scale(v2.getVec0());
        Vector vec2Temp = v2.getVector().copy().scale(this.getVec0());
        
        int[] zeroN1N2 = new int[n1*n2];
        int[] zeroN1masN2 = new int[n1+n2];
        int[] zeroN1 = new int[n1];
        int[] zeroN2 = new int[n2];
        for(int i = 0; i < n1*n2; i++) zeroN1N2[i] = i;
        for(int i = 0; i < n1+n2; i++) zeroN1masN2[i] = n1*n2+i;
        for(int i = 0; i < n1; i++) zeroN1[i] = n1*n2+i;
        for(int i = 0; i < n2; i++) zeroN2[i] = n1*n2+n1+i;
        
        
        res.setVector(MatrixUtils.concatVectors(
        		MatrixUtils.kroneckerVectors(this.getVector(), v2.getVector(), 
        				Matrices.getSubVector(res.getVector(), zeroN1N2)),
        		MatrixUtils.concatVectors(vec1Temp, vec2Temp, Matrices.getSubVector(res.getVector(), zeroN1masN2)),
        		res.getVector()
        		)
        		);

        Matrix sumAB = MatrixUtils.kroneckerSum(this.getMatrix(), v2.getMatrix(), 
        		Matrices.getSubMatrix(res.getMatrix(),zeroN1N2,zeroN1N2));
        Matrix A00B = MatrixUtils.concatQuad(
        		this.getMatrix(), 	Matrices.getSubMatrix(res.getMatrix(), zeroN1 , zeroN2).copy().zero(), 
        		Matrices.getSubMatrix(res.getMatrix(), zeroN2 , zeroN1) , 	v2.getMatrix(),
        		Matrices.getSubMatrix(res.getMatrix(), zeroN1masN2, zeroN1masN2)
        		);
        Matrix matRU = MatrixUtils.concatCols(
        							MatrixUtils.kronecker(Matrices.identity(n1),v2.getMat0(),
        									Matrices.getSubMatrix(res.getMatrix(),zeroN1N2,zeroN1)),
        							MatrixUtils.kronecker(this.getMat0(),Matrices.identity(n2),
                							Matrices.getSubMatrix(res.getMatrix(),zeroN1N2,zeroN2)),
                					Matrices.getSubMatrix(res.getMatrix(),zeroN1N2,zeroN1masN2)
        				);
        
        res.setMatrix(MatrixUtils.concatQuad( sumAB, matRU, 
        				Matrices.getSubMatrix(res.getMatrix(),zeroN1masN2,zeroN1N2).copy().zero(), A00B,
        				res.getMatrix())
        			);
        return res;
	}
    
    /**
     * @see jphase.DiscPhaseVar#max(jphase.DiscPhaseVar)
     */
    public DiscPhaseVar max(DiscPhaseVar v2) {
        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();
        DiscPhaseVar res = this.newVar(n1+n2+n1*n2);
        return this.max(v2, res);
    }
	
	/**
     * @see jphase.DiscPhaseVar#toString()
     */
	@Override
    public final String toString() {
        return label();
    }

    public String label() {
        return "DPH - " + getNumPhases() + " Phases";
    }

    public String description() {
		String s = "\n\n__________________________________________________\n";
		s += "Phase-Type Distribution";
		s+="\nNumber of Phases: "+this.getNumPhases()+"\n";
		s=s+"Vector:\n\t";
		for(int i = 0 ; i < this.getNumPhases(); i++){
			s+=String.format("%5.4f",this.getVector().get(i));
			s+="\t";
		}
		s=s+"\nMatrix:\n";
		for(int i = 0 ; i < this.getNumPhases(); i++){
			s+="\t";
			for(int j = 0 ; j < this.getNumPhases(); j++){
				s+=String.format("%7.4f",this.getMatrix().get(i,j));
				s+="\t";
			}
			s+="\n";
		}
		s+="__________________________________________________\n";
		return s;
	}
}
