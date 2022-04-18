package jphase;

import static jphase.Utils.*;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.sparse.CG;
import no.uib.cipr.matrix.sparse.DefaultIterationMonitor;
import no.uib.cipr.matrix.sparse.IterativeSolver;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
/**
 * Abstract class that defines a  Continuous Phase-Type Distribution
 * @author German Riaño
 * @author Juan F. Pérez
 * @version 1.0
 */
public abstract class AbstractContPhaseVar implements ContPhaseVar {

    /**
     * @see jphase.ContPhaseVar#sumPH(jphase.DiscPhaseVar,
     *      jphase.ContPhaseVar)
     */
    public ContPhaseVar sumPH(DiscPhaseVar B, ContPhaseVar res) {
        int n1 = this.getNumPhases();
        int n2 = B.getNumPhases();
        int n3 = res.getNumPhases();
        if (n1 * n2 != n3)
            throw new IndexOutOfBoundsException(
                    "The closure operation cannot be done:\n "
                            + "this.Phases * B.Phases != res.Phases (" + n1
                            + " * " + n2 + " != " + n3 + ")");
        
        
        // ISInv = (I-alpha_(m+1)S)^(-1)
        Matrix ISInv = Matrices.identity(n2)
                .add(-this.getVec0(), B.getMatrix());
        //Inverse 
        ISInv = ISInv.solve(Matrices.identity(n2), ISInv.copy());

        res.setVector(MatrixUtils.kroneckerVectors(this.getVector(), ISInv
                .transMult(B.getVector(), B.getVector().copy()), res
                .getVector()));
        // L1 = T kron I
        Matrix L1 = MatrixUtils.kronecker(this.getMatrix(), Matrices
                .identity(n2), res.getMatrix().copy());
        // L2 = (1-alpha_(m+1))t*alpha k ISInv*S
        Matrix L2 = MatrixUtils.kronecker(MatrixUtils.multVector(
                this.getMat0(), this.getVector(), this.getMatrix().copy())
                .scale(1 - this.getVec0()), ISInv.mult(B.getMatrix(), B
                .getMatrix().copy()), res.getMatrix().copy());
        res.setMatrix(L1.add(L2));
        return res;
    }
    
    /**
     * @see jphase.ContPhaseVar#sumPH(jphase.DiscPhaseVar)
     */
    public ContPhaseVar sumPH(DiscPhaseVar v2) {
        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();
        ContPhaseVar res = this.newVar(n1*n2);
        return this.sumPH(v2, res);
    }

    /**
     * @see jphase.ContPhaseVar#getNumPhases()
     */
    public int getNumPhases() {
        return this.getVector().size();
    }

    /**
     * @see jphase.ContPhaseVar#getVec0()
     */
    public double getVec0() {
        return 1.0 - this.getVector().dot(
                MatrixUtils.OnesVector(this.getNumPhases()));
    }

    /**
     * @see jphase.ContPhaseVar#getMat0()
     */
    public Vector getMat0() {
        return this.getMatrix().mult(-1,
                MatrixUtils.OnesVector(this.getNumPhases()),
                this.getVector().copy());
    }

    /**
     * @see jphase.ContPhaseVar#getMatrixArray()
     */
    public double[][] getMatrixArray() {
        return Matrices.getArray(this.getMatrix());
    }

    /**
     * @see jphase.ContPhaseVar#getVectorArray()
     */
    public double[] getVectorArray() {
        return Matrices.getArray(this.getVector());
    }

    /**
     * @see jphase.ContPhaseVar#getMat0Array()
     */
    public double[] getMat0Array() {
        return Matrices.getArray(this.getMat0());
    }

    /**
     * @see jphase.ContPhaseVar#expectedValue()
     */
    public double expectedValue() {
        return moment(1);
    }

    /**
     * @see jphase.ContPhaseVar#variance()
     */
    public double variance() {
        double m1 = moment(1);
        return moment(2) - (m1 * m1);
    }

    /**
     * @see jphase.ContPhaseVar#stdDeviation()
     */
    public double stdDeviation() {
        return Math.sqrt(variance());
    }

    /**
     * @see jphase.ContPhaseVar#CV()
     */
    public double CV() {
        double m = expectedValue();
        return moment(2) / (m * m) - 1;
    }

    /**
     * @see jphase.ContPhaseVar#moment(int)
     */
    public double moment(int k) {
        double res = 0;
        
        Vector x = this.getVector().copy().zero();
        IterativeSolver solv = new CG(x);
        Matrix A1 = MatrixUtils.matPower( this.getMatrix(), k).transpose();
        Matrix A = A1.copy().zero();
        A = A1.transAmult(A1, A);
        Vector b = A1.transMult(this.getVector(), this.getVector().copy());
        solv.setIterationMonitor(new DefaultIterationMonitor(100000, 1e-10,  1e-50, 1e+5));
        
        try{
            solv.solve(A, b, x);
            res = x.dot(MatrixUtils.OnesVector(this.getVector().size()));
        }catch(IterativeSolverNotConvergedException e){
            System.out.println("Imposible to calculate moment "+k+"("+e.getReason()+")");
        }
      
        if (k % 2 != 0)
            res = -res;
        res *= fact(k);
        return res;
    }

    /**
     * @see jphase.ContPhaseVar#cdf(double)
     */
    public double cdf(double x) {
        return 1.0 - MatrixUtils.expTimesOnes(this.getMatrix(), x, this
                .getVector());
    }

    /**
     * @see jphase.ContPhaseVar#cdf(int, double)
     */
    public double[] cdf(int n, double delta) {
        double[] result = new double[n];
        double survivals[] = MatrixUtils.expTimesOnes(this.getMatrix(), n,
                delta, this.getVector());
        for (int i = 0; i < n; i++) {
            result[i] = 1.0 - survivals[i];
        }
        return result;
    }

    /**
     * @see jphase.ContPhaseVar#pdf(double)
     */
    public double pdf(double x) {
        return MatrixUtils.exp(this.getMatrix(), x, this.getVector(), this
                .getMat0());
    }

    /**
     * @see jphase.ContPhaseVar#pdf(int, double)
     */
    public double[] pdf(int n, double delta) {
        return MatrixUtils.exp(this.getMatrix(), n, delta, this.getVector(),
                this.getMat0(), true);
    }

    /**
     * @see jphase.ContPhaseVar#prob(double, double)
     */
    public double prob(double a, double b) {
        return (b > a) ? cdf(b) - cdf(a) : 0.0;
    }

    /**
     * @see jphase.ContPhaseVar#survival(double)
     */
    public double survival(double x) {
        return (1.0 - this.cdf(x));
    }

    /**
     * @see jphase.ContPhaseVar#survival(int, double)
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
    public double lossFunction1(double t){
        double res = 0;
   
        Vector x = this.getVector().copy().zero();
        IterativeSolver solv = new CG(x);
        
        Matrix A = this.getMatrix().copy();
        Matrix At = A.copy().zero();
        A.transpose(At);
        Matrix B = A.copy().zero();
        
        B = A.mult(At, B);
        Vector temp = x.copy(); 
        temp = MatrixUtils.exp(this.getMatrix(), t).mult(
                MatrixUtils.OnesVector(this.getNumPhases()), temp.copy());
                        
        Vector b = A.mult(this.getVector(), x.copy().zero());
        
        solv.setIterationMonitor(new DefaultIterationMonitor(100000, 1e-10,  1e-50, 1e+5));
        
        try{
            solv.solve(B, b, x);
            res = x.dot(temp);
            
            
        }catch(IterativeSolverNotConvergedException e){
            //e.printStackTrace();
            System.out.println("Imposible to calculate order-1 loss Function at "+t+"("+e.getReason()+")");
        }
        return -res;
    }
    
    /**
     * @see jphase.PhaseVar#lossFunction2(double)
     */
    public double lossFunction2(double t){
        double res = 0;
           
        Vector x = this.getVector().copy().zero();
        IterativeSolver solv = new CG(x);
        
        Matrix A = MatrixUtils.matPower(this.getMatrix().copy(), 2);
        Matrix At = A.copy();
        A.transpose(At);
        Matrix B = A.copy().zero();
        B = At.mult(A, B);
        Vector temp = x.copy(); 
        temp = MatrixUtils.exp(this.getMatrix(), t).mult(
                MatrixUtils.OnesVector(this.getNumPhases()), temp);
                        
        Vector b = At.mult(temp, x.copy());
        
        solv.setIterationMonitor(new DefaultIterationMonitor(100000, 1e-10,  1e-50, 1e+5));
        
        try{
            solv.solve(B, b, x);
            res = this.getVector().dot(x);
            
            
        }catch(IterativeSolverNotConvergedException e){
            //e.printStackTrace();
            System.out.println("Imposible to calculate order-2 loss Function at "+t+"("+e.getReason()+")");
        }
        return res;
    }
    
    /**
     * @see jphase.ContPhaseVar#quantil(double)
     */
    public double quantil(double p) {
        double x, dif, derv;
        int ITMAX = 100;
        int cnt;
        x = this.expectedValue();
        for (cnt = 0; cnt < ITMAX; cnt++) {
            derv = this.pdf(x);
            dif = this.cdf(x) - p;
            if (Math.abs(dif) < 1.0e-10)
                return x;
            x = Math.max(x - dif / derv, 0);
        }
        return 0.0;
    }

    /**
     * @see jphase.ContPhaseVar#median()
     */
    public double median() {
        return this.quantil(0.5);
    }

    /**
     * @see jphase.ContPhaseVar#sum(jphase.ContPhaseVar, jphase.ContPhaseVar)
     */
    public ContPhaseVar sum(ContPhaseVar v2, ContPhaseVar res) {
        if (Matrices.cardinality(this.getVector()) == 0)
            return v2;
        if (Matrices.cardinality(v2.getVector()) == 0)
            return this;

        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();

        Vector A0 = this.getMat0();

        Vector vec1 = this.getVector();
        double vec1_0 = this.getVec0();
        Vector vec2 = v2.getVector();

        Vector vec2Temp = vec2.copy();
        vec2Temp.scale(vec1_0);

        res.setVector(MatrixUtils
                .concatVectors(vec1, vec2Temp, res.getVector()));

        int[] zeroRowsLD = new int[n2];
        int[] zeroColsLD = new int[n1];
        for (int i = 0; i < n2; i++)
            zeroRowsLD[i] = n1 + i;
        for (int i = 0; i < n1; i++)
            zeroColsLD[i] = i;

        int[] zeroRowsRU = new int[n1];
        int[] zeroColsRU = new int[n2];
        for (int i = 0; i < n1; i++)
            zeroRowsRU[i] = i;
        for (int i = 0; i < n2; i++)
            zeroColsRU[i] = n1 + i;

        res.setMatrix(MatrixUtils.concatQuad(this.getMatrix(), MatrixUtils
                .multVector(A0, vec2, Matrices.getSubMatrix(res.getMatrix(),
                        zeroRowsRU, zeroColsRU)), Matrices.getSubMatrix(res
                .getMatrix(), zeroRowsLD, zeroColsLD), v2.getMatrix(), res
                .getMatrix()));
        return res;
    }
    
    /**
     * @see jphase.ContPhaseVar#sum(jphase.ContPhaseVar)
     */
    public ContPhaseVar sum(ContPhaseVar v2) {
        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();
        ContPhaseVar res = this.newVar(n1+n2);
        return this.sum(v2, res);
    }

    
    /**
     * @see jphase.ContPhaseVar#sumGeom(double)
     */
    public ContPhaseVar sumGeom(double p) {
        ContPhaseVar res = this.copy();
        res.setVector(this.getVector());
        res.setMatrix(this.getMatrix().copy().add(
                1 - p,
                MatrixUtils.multVector(this.getMat0(), this.getVector(), res
                        .getMatrix().copy().zero())));
        return res;
    }
    
    /**
     * @see jphase.ContPhaseVar#mix(double, jphase.ContPhaseVar, jphase.ContPhaseVar)
     */
    public ContPhaseVar mix(double p, ContPhaseVar v2, ContPhaseVar res) {
        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();
        Vector vec1Temp = this.getVector().copy().scale(p);
        Vector vec2Temp = v2.getVector().copy().scale(1 - p);

        res.setVector(MatrixUtils.concatVectors(vec1Temp, vec2Temp, res
                .getVector()));

        int[] zeroRowsLD = new int[n2];
        int[] zeroColsLD = new int[n1];
        for (int i = 0; i < n2; i++)
            zeroRowsLD[i] = n1 + i;
        for (int i = 0; i < n1; i++)
            zeroColsLD[i] = i;

        int[] zeroRowsRU = new int[n1];
        int[] zeroColsRU = new int[n2];
        for (int i = 0; i < n1; i++)
            zeroRowsRU[i] = i;
        for (int i = 0; i < n2; i++)
            zeroColsRU[i] = n1 + i;

        res.setMatrix(MatrixUtils.concatQuad(this.getMatrix(), Matrices
                .getSubMatrix(res.getMatrix(), zeroRowsRU, zeroColsRU),
                Matrices.getSubMatrix(res.getMatrix(), zeroRowsLD, zeroColsLD),
                v2.getMatrix(), res.getMatrix()));

        return res;
    }
    
    /**
     * @see jphase.ContPhaseVar#mix(double, jphase.ContPhaseVar)
     */
    public ContPhaseVar mix(double p, ContPhaseVar v2) {
        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();
        ContPhaseVar res = this.newVar(n1+n2);
        return this.mix(p, v2, res);
    }
    

    
    /**
     * @see jphase.ContPhaseVar#min(jphase.ContPhaseVar, jphase.ContPhaseVar)
     */
    public ContPhaseVar min(ContPhaseVar v2, ContPhaseVar res) {
        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();
        int n3 = res.getNumPhases();
        if (n1 * n2 != n3)
            throw new IndexOutOfBoundsException(
                    "The closure operation cannot be done:\n n1*n2 != n3(" + n1
                            + " * " + n2 + " != " + n3 + ")");

        res.setVector(MatrixUtils.kroneckerVectors(this.getVector(), v2
                .getVector(), res.getVector()));
        res.setMatrix(MatrixUtils.kroneckerSum(this.getMatrix(),
                v2.getMatrix(), res.getMatrix()));
        return res;
    }
    
    /**
     * @see jphase.ContPhaseVar#min(jphase.ContPhaseVar)
     */
    public ContPhaseVar min(ContPhaseVar v2) {
        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();
        ContPhaseVar res = this.newVar(n1*n2);
        return this.min(v2, res);
    }

    /**
     * @see jphase.ContPhaseVar#max(jphase.ContPhaseVar, jphase.ContPhaseVar)
     */
    public ContPhaseVar max(ContPhaseVar v2, ContPhaseVar res) {
        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();
        int n3 = res.getNumPhases();
        if (n1 * n2 + n1 + n2 != n3)
            throw new IndexOutOfBoundsException(
                    "The closure operation cannot be done:\n n1*n2 + n1 + n2 != n3("
                            + n1 + "*" + n2 + " + " + n1 + " + " + n2 + " != "
                            + n3 + ")");

        Vector vec1Temp = this.getVector().copy().scale(v2.getVec0());
        Vector vec2Temp = v2.getVector().copy().scale(this.getVec0());

        int[] zeroN1N2 = new int[n1 * n2];
        int[] zeroN1masN2 = new int[n1 + n2];
        int[] zeroN1 = new int[n1];
        int[] zeroN2 = new int[n2];
        for (int i = 0; i < n1 * n2; i++)
            zeroN1N2[i] = i;
        for (int i = 0; i < n1 + n2; i++)
            zeroN1masN2[i] = n1 * n2 + i;
        for (int i = 0; i < n1; i++)
            zeroN1[i] = n1 * n2 + i;
        for (int i = 0; i < n2; i++)
            zeroN2[i] = n1 * n2 + n1 + i;

        res.setVector(MatrixUtils.concatVectors(MatrixUtils.kroneckerVectors(
                this.getVector(), v2.getVector(), Matrices.getSubVector(res
                        .getVector(), zeroN1N2)), MatrixUtils.concatVectors(
                vec1Temp, vec2Temp, Matrices.getSubVector(res.getVector(),
                        zeroN1masN2)), res.getVector()));

        Matrix sumAB = MatrixUtils.kroneckerSum(this.getMatrix(), v2
                .getMatrix(), Matrices.getSubMatrix(res.getMatrix(), zeroN1N2,
                zeroN1N2));
        Matrix A00B = MatrixUtils.concatQuad(this.getMatrix(), Matrices
                .getSubMatrix(res.getMatrix(), zeroN1, zeroN2).copy().zero(),
                Matrices.getSubMatrix(res.getMatrix(), zeroN2, zeroN1), v2
                        .getMatrix(), Matrices.getSubMatrix(res.getMatrix(),
                        zeroN1masN2, zeroN1masN2));
        Matrix matRU = MatrixUtils.concatCols(MatrixUtils.kronecker(Matrices
                .identity(n1), v2.getMat0(), Matrices.getSubMatrix(res
                .getMatrix(), zeroN1N2, zeroN1)), MatrixUtils.kronecker(this
                .getMat0(), Matrices.identity(n2), Matrices.getSubMatrix(res
                .getMatrix(), zeroN1N2, zeroN2)), Matrices.getSubMatrix(res
                .getMatrix(), zeroN1N2, zeroN1masN2));
        res.setMatrix(MatrixUtils.concatQuad(sumAB, matRU, Matrices
                .getSubMatrix(res.getMatrix(), zeroN1masN2, zeroN1N2).copy()
                .zero(), A00B, res.getMatrix()));
        return res;
    }

    /**
     * @see jphase.ContPhaseVar#max(jphase.ContPhaseVar)
     */
    public ContPhaseVar max(ContPhaseVar v2) {
        int n1 = this.getNumPhases();
        int n2 = v2.getNumPhases();
        ContPhaseVar res = this.newVar(n1+n2+n1*n2);
        return this.max(v2, res);
    }
    
    /**
     * @see jphase.ContPhaseVar#times(double)
     */
    public ContPhaseVar times(double c) {
        ContPhaseVar res = this.copy();
        res.setVector(this.getVector());
        res.setMatrix(this.getMatrix());
        res.getMatrix().scale(1.0 / c);
        return res;
    }

    /**
     * @see jphase.ContPhaseVar#residualTime(double)
     */
    public ContPhaseVar residualTime(double x) {
        ContPhaseVar res = this.copy();
        double denom = this.survival(x);
        Vector vecTemp = this.getVector().copy().scale(1.0 / denom);
        res.setVector(MatrixUtils.exp(this.getMatrix(), x).transMult(vecTemp,
                res.getVector()));

        res.setMatrix(this.getMatrix());

        return res;
    }

    /**
     * @see jphase.ContPhaseVar#eqResidualTime()
     */
    public ContPhaseVar eqResidualTime() {
        ContPhaseVar res = this.copy();
        double mu = this.expectedValue();
        Vector pi = this.getVector().copy().zero();
        IterativeSolver solv = new CG(pi);

        Matrix A1 = this.getMatrix().copy().transpose();
        Matrix A = A1.copy();
        A = A1.transAmult(A1, A);
        
        Vector b = A1.transMult(this.getVector(), this.getVector().copy());
        try{
            solv.solve(A, b, pi);
                
        }catch(IterativeSolverNotConvergedException e){
            e.printStackTrace();
        }
        
        res.setVector(pi.scale(-1.0/mu));
        res.setMatrix(this.getMatrix());
        return res;
    }


    /**
     * @see jphase.ContPhaseVar#waitingQ(double)
     */
    public ContPhaseVar waitingQ(double rho) {
        ContPhaseVar res = this.copy();
        if (rho >= 1.0)
            throw new IllegalArgumentException("rho (" + rho + ") >= 1\n"
                    + "non-stable queue");

        double mu = this.expectedValue();
        Vector pi = this.getVector().copy().zero();
        IterativeSolver solv = new CG(pi);
        
        Matrix At = this.getMatrix().copy().transpose();
        Matrix A = At.copy();
        A = At.transAmult(At, A);
        Vector b = At.transMult(this.getVector(), this.getVector().copy());

        try{
            solv.solve(A, b, pi);
                
        }catch(IterativeSolverNotConvergedException e){
            e.printStackTrace();
        }
        
        res.setVector(pi.scale(-rho/mu));
        res.setMatrix(this.getMatrix().copy().add(
                MatrixUtils.multVector(this.getMat0(), res.getVector(), res
                        .getMatrix().copy())));

        return res;
    }
    
    /**
     * @see jphase.ContPhaseVar#residualVar(double)
     */    
    public ContPhaseVar residualVar(double a){
        ContPhaseVar res = this.copy();
        res.setVector( (MatrixUtils.exp(this.getMatrix(), a)).transMult(this.getVector(),this.getVector().copy())  );
        res.setMatrix(this.getMatrix());
        return res;
    }
    

    /**
     * @see jphase.ContPhaseVar#toString()
     */
    @Override
    public final String toString() {
        return description();
    }

    public String label() {
        return "CPH - " + getNumPhases() + " Phases";
    }

    public String description() {
        String s = "__________________________________________________\n";
        s += "Phase-Type Distribution";
        s += "\nNumber of Phases: " + this.getNumPhases() + "\n";
        s = s + "Vector:\n\t";
        for (int i = 0; i < this.getNumPhases(); i++) {
            s += String.format("%6.4f", this.getVector().get(i));
            s += "\t";
        }
        s = s + "\nMatrix:\n";
        for (int i = 0; i < this.getNumPhases(); i++) {
            s += "\t";
            for (int j = 0; j < this.getNumPhases(); j++) {
                s += String.format("%6.4f", this.getMatrix().get(i, j));
                s += "\t";
            }
            s += "\n";
        }
        s += "__________________________________________________\n";
        return s;
    }
}
