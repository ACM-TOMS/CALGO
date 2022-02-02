package jphase;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;
import no.uib.cipr.matrix.sparse.SparseVector;

/**
 * This class allows the creation and manipulation 
 * of Continuous Phase-type distributions represented by 
 * sparse (Flexible Compressed Row) matrices. 
 * @author German Riaño
 * @author Juan F. Perez
 * @version 1.0
 * 
 */
public class SparseContPhaseVar extends AbstractContPhaseVar{
	
	
	/**
	 *  Rate Matrix in Sparse representation (CompRowMatrix)
	 */
	protected FlexCompRowMatrix A;

	/**
	 *  Initial Probability distribution vector
	 */
	protected SparseVector alpha;
	
	
    /**
	 * Constructs an empty Continuous Phase-type Distribution 
	 * of size n with sparse representation (CompRowMatrix)
	 * @param n size of the Continuous Phase-type Distribution  
	 */
    public SparseContPhaseVar(int n) {
        this.A = new FlexCompRowMatrix(n,n);
        this.alpha = new SparseVector(n);
    }
    
	/**
	 * Constructs a continuous Phase-type Distribution 
	 * with sparse representation (CompRowMatrix)
	 * @param A rate matrix
	 * @param alpha initial probability distribution vector
	 */
	public SparseContPhaseVar(SparseVector alpha, FlexCompRowMatrix A) {
		this.A = A;
		this.alpha = alpha;
	}
	
	/**
	 * Constructs a continuous Phase-type Distribution 
	 * with sparse representation (CompRowMatrix)
	 * @param A rate matrix
	 * @param alpha initial probability distribution vector
	 */
	public SparseContPhaseVar(Vector alpha, Matrix A) {
		this.A = (FlexCompRowMatrix)A.copy();
		this.alpha = (SparseVector)alpha.copy();
	}
	
	/**
	 * Constructs a continuous Phase-type Distribution 
	 * with sparse representation (CompRowMatrix)
	 * @param A rate matrix
	 * @param alpha initial probability distribution vector
	 */
	public SparseContPhaseVar(double[] alpha, double[][] A) {
		this.A = new FlexCompRowMatrix( new DenseMatrix(A) );
		this.alpha = new SparseVector( new DenseVector(alpha));
	}

	public Matrix getMatrix() {
		return this.A;
	}

	public void setMatrix(Matrix A) {
		this.A = (FlexCompRowMatrix)A.copy();
	}

	public Vector getVector() {
		return this.alpha;
	}

	public void setVector(Vector alpha) {
		this.alpha = (SparseVector)alpha.copy();
	}

	public ContPhaseVar copy() {
		return new SparseContPhaseVar(
				this.getVector(), this.getMatrix());
	}

	public ContPhaseVar newVar(int n){
        return new SparseContPhaseVar(n);
    }

}
