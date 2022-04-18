package jphase;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;
import no.uib.cipr.matrix.sparse.SparseVector;

/**
 * This class allows the creation and manipulation 
 * of Discrete Phase-type distributions represented by 
 * sparse (Flexible Compressed Row) matrices. 
 * @author German Riaño
 * @author Juan F. Perez
 * @version 1.0
 * 
 */
public class SparseDiscPhaseVar extends AbstractDiscPhaseVar{
	
	
	/**
	 *  Transition Matrix in Sparse representation (FlexCompRowMatrix)
	 */
	protected FlexCompRowMatrix A;

	/**
	 *  Initial Probability distribution vector
	 */
	protected SparseVector alpha;
	
	
    /**
	 * Constructs an empty Discrete Phase-type Distribution 
	 * of size n with sparse representation (FlexCompRowMatrix)
	 * @param n size of the Discrete Phase-type Distribution  
	 */
    public SparseDiscPhaseVar(int n) {
    	this.A = new FlexCompRowMatrix(n,n);
        this.alpha = new SparseVector(n);
    }
    
	/**
	 * Constructs a discrete Phase-type Distribution 
	 * with sparse representation (FlexCompRowMatrix)
	 * @param A transition probability matrix
	 * @param alpha initial probability distribution vector
	 */
	public SparseDiscPhaseVar(SparseVector alpha, FlexCompRowMatrix A) {
		this.A = A;
		this.alpha = alpha;
	}
	
	/**
	 * Constructs a discrete Phase-type Distribution 
	 * with sparse representation (FlexCompRowMatrix)
	 * @param A transition probability matrix
	 * @param alpha initial probability distribution vector
	 */
	public SparseDiscPhaseVar(Vector alpha, Matrix A) {
		this.A = (FlexCompRowMatrix)A.copy();
		this.alpha = (SparseVector)alpha.copy();
	}
	
	/**
	 * Constructs a discrete Phase-type Distribution 
	 * with sparse representation (FlexCompRowMatrix)
	 * @param A transition probability matrix
	 * @param alpha initial probability distribution vector
	 */
	public SparseDiscPhaseVar(double[] alpha, double[][] A) {
		this.A = new FlexCompRowMatrix( new DenseMatrix(A) );
		this.alpha = new SparseVector( new DenseVector(alpha));
	}
	
    /**
     * @see jphase.PhaseVar#getMatrix()
     */
	public Matrix getMatrix() {
		return this.A;
	}

    /**
     * @see jphase.PhaseVar#setMatrix(no.uib.cipr.matrix.Matrix)
     */
	public void setMatrix(Matrix A) {
		this.A = (FlexCompRowMatrix)A.copy();
	}

    /**
     * @see jphase.PhaseVar#getVector()
     */
	public Vector getVector() {
		return this.alpha;
	}

    /**
     * @see jphase.PhaseVar#setVector(no.uib.cipr.matrix.Vector)
     */
	public void setVector(Vector alpha) {
		this.alpha = (SparseVector)alpha.copy();
	}
	
    /**
     * @see jphase.ContPhaseVar#copy()
     */
	public DiscPhaseVar copy() {
		return new SparseDiscPhaseVar(
				this.getVector(), this.getMatrix());
	}
   
	public DiscPhaseVar newVar(int n){
        return new SparseDiscPhaseVar(n);
    }

}
