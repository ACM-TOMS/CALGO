package jphase;

import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;

/**
 * This class allows the creation and manipulation of Discrete
 * Phase-type distributions represented by dense matrices.
 * @author German Riaño
 * @author Juan F. Perez
 * @version 0.1 
 */
public class DenseDiscPhaseVar extends AbstractDiscPhaseVar implements
        DiscPhaseVar {

    /**
     * Probability Transition Matrix
     */
    protected DenseMatrix A;

    /**
     * Initial Probability distribution vector
     */
    protected DenseVector alpha;

    /**
     * Constructs an empty Discrete Phase-type Distribution with dense
     * representation
     */
    public DenseDiscPhaseVar() {
        this.A = new DenseMatrix(new double[][] { { 0.0 } });
        this.alpha = new DenseVector(new double[] { 0 });
    }

    /**
     * Constructs an empty Discrete Phase-type Distribution of size n
     * with dense representation
     * @param n size of the Discrete Phase-type Distribution
     */
    public DenseDiscPhaseVar(int n) {
        this.A = new DenseMatrix(n, n);
        this.alpha = new DenseVector(n);
    }

    /**
     * Constructs a Discrete Phase-type Distribution with dense
     * representation
     * @param A transition probability matrix
     * @param alpha initial probability distribution vector
     */
    public DenseDiscPhaseVar(DenseVector alpha, DenseMatrix A) {
        this.A = A;
        this.alpha = alpha;
    }

    /**
     * Constructs a Discrete Phase-type Distribution with dense
     * representation
     * @param A transition probability matrix
     * @param alpha initial probability distribution vector
     */
    public DenseDiscPhaseVar(Vector alpha, Matrix A) {
        this.A = (DenseMatrix) A.copy();
        this.alpha = (DenseVector) alpha.copy();
    }

    /**
     * Constructs a Discrete Phase-type Distribution with dense
     * representation
     * @param A transition probability matrix
     * @param alpha initial probability distribution vector
     */
    public DenseDiscPhaseVar(double[] alpha, double[][] A) {
        int n1, n2;
        n1 = alpha.length;
        n2 = A.length;
        if (n1 != n2) {
            throw new IndexOutOfBoundsException(
                    "Vector and matrix not of the same size");
        } else {
            DenseMatrix Am = new DenseMatrix(A);
            DenseVector Vm = new DenseVector(alpha);
            this.A = Am;
            this.alpha = Vm;
        }
    }

    /**
     * Discrete Phase distribution that represents a geometric
     * distribution with probability of success p
     * @param p probability of success
     * @return Dense Discrete Phase-Type Distribution
     */
    public static DenseDiscPhaseVar Geom(double p) {
        double[][] matriz = new double[1][1];
        matriz[0][0] = 1 - p;
        DenseMatrix A = new DenseMatrix(matriz);

        double[] vector = new double[1];
        vector[0] = 1;
        DenseVector alpha = new DenseVector(vector);
        return new DenseDiscPhaseVar(alpha, A);
    }

    /**
     * Discrete Phase Distribution that represents a Negative
     * Binomial distribution with parameters p and r
     * @param p probability of success in one trial
     * @param r number of successes until absorption
     * @return Dense Discrete Phase-Type distribution
     */
    public static DenseDiscPhaseVar NegativeBinomial(double p, int r) {
        double[][] matriz = new double[r][r];
        for (int i = 0; i < r; i++) {
            matriz[i][i] = 1 - p;
            if (i < r - 1)
                matriz[i][i + 1] = p;
        }

        DenseMatrix A = new DenseMatrix(matriz);

        double[] vector = new double[r];
        vector[0] = 1;
        DenseVector alpha = new DenseVector(vector);
        return new DenseDiscPhaseVar(alpha, A);
    }

    /**
     * @see jphase.PhaseVar#getMatrix()
     */
    public Matrix getMatrix() {
        return this.A;
    }

    /**
     * @see jphase.PhaseVar#getVector()
     */
    public Vector getVector() {
        return this.alpha;
    }

    /**
     * @see jphase.PhaseVar#setMatrix(no.uib.cipr.matrix.Matrix)
     */
    public void setMatrix(Matrix A) {
        this.A = (DenseMatrix) A;
    }

    /**
     * @see jphase.PhaseVar#setVector(no.uib.cipr.matrix.Vector)
     */
    public void setVector(Vector alpha) {
        this.alpha = (DenseVector) alpha;
    }

    /**
     * @see jphase.ContPhaseVar#copy()
     */
    public DiscPhaseVar copy() {
        return new DenseDiscPhaseVar(this.getVector(), this.getMatrix());
    }

    public DiscPhaseVar newVar(int n) {
        return new DenseDiscPhaseVar(n);
    }

}
