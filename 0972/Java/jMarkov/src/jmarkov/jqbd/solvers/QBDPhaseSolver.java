package jmarkov.jqbd.solvers;

import java.util.ArrayList;

import jphase.MatrixUtils;
import jphase.PhaseVar;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

/**
 * This class is abstrac and inplements the mathods necesaries
 * to compute the performance measures for a PHPH1 model 
 * and generate the caracteristic matrix from it.
 * @author Andrés Sarmiento Romero
 */
public abstract class QBDPhaseSolver {

	/**
	 * R Matrix
	 */
	protected DenseMatrix R;

	/**
	 * A_0 Matrix For a NON-Homogeneous QBD process 
	 */
	protected Matrix MA0;

	/**
	 * A_1 Matrix For a NON-Homogeneous QBD process 
	 */
	protected Matrix MA1;

	/**
	 * A_2 Matrix For a NON-Homogeneous QBD process 
	 */
	protected Matrix MA2;

	/**
	 * B_0 Matrix For a NON-Homogeneous QBD process 
	 */
	protected Matrix MB0;

	/**
	 * B_1 Matrix For a NON-Homogeneous QBD process 
	 */
	protected Matrix MB1;

	/**
	 * B_2 Matrix For a NON-Homogeneous QBD process 
	 */
	protected Matrix MB2;

	/**
	 * N Matrix For a QBD process 
	 */
	protected DenseMatrix N;

	/**
	 * U Matrix For a QBD process 
	 */
	protected DenseMatrix U;

	/**
	 * G Matrix For a QBD process 
	 */
	protected DenseMatrix G;

	/**
	 * Steady state probabilities for the level
	 */    
	protected ArrayList<Double> pis;

	/**
	 * Steady state probabilities specifically for each level
	 */    
	protected ArrayList<DenseMatrix> pisP;

	/**
	 * Phase Type of the inter-arrival times
	 */
	protected PhaseVar arr;

	/**
	 * Phase Type of the service times
	 */
	protected PhaseVar serv;

	/**
	 * Arrival rate
	 */
	private double lambdaP;

	/**
	 * This method computes the 6 QBD caracteristic matrices
	 * using the interarrival and service PH distributions
	 * @param arrivals PH that represents the inter-arrival times
	 * @param service PH that represents the service times
	 */
	public void initilicePHPH1(PhaseVar arrivals, PhaseVar service) {

		arr = arrivals;
		serv = service;

		Vector alpha = arrivals.getVector();
		Vector beta = service.getVector();

		Matrix A = arrivals.getMatrix();
		Matrix M = service.getMatrix();

		int tamA = A.numRows();
		int tamM = M.numRows();

		double[] era = new double[tamA];
		double[] erm = new double[tamM];
		for(int i = 0; i != tamA ; i++)
			era[i] = 1;
		for(int i = 0; i != tamM ; i++)
			erm[i] = 1;

		DenseVector lambda = negativeMatVecMult(A, era, new DenseVector(tamA));
		DenseVector muV = negativeMatVecMult(M, erm, new DenseVector(tamM));		

		Matrix IA = new DenseMatrix(Matrices.identity(tamA));
		Matrix IM = new DenseMatrix(Matrices.identity(tamM));

		MA0 = new DenseMatrix(tamA*tamM, tamA*tamM);
		MA2 = new DenseMatrix(tamA*tamM, tamA*tamM);

		MB0 = new DenseMatrix(tamA, tamA*tamM);
		MB2 = new DenseMatrix(tamA*tamM, tamA);

		MA0 = MatrixUtils.kronecker(multVector(lambda, new DenseVector(alpha), new DenseMatrix(lambda.size(), alpha.size())), IM, MA0);
		MA1 = MatrixUtils.kroneckerSum(A, M);
		MA2 = MatrixUtils.kronecker(IA , multVector(muV, new DenseVector(beta), new DenseMatrix(muV.size(), beta.size())), MA2);

		MB0 = MatrixUtils.kronecker(multVector(lambda, new DenseVector(alpha), new DenseMatrix(lambda.size(), alpha.size())), rowVectoraMat(new DenseVector(beta)), MB0);
		MB1 = A.copy();
		MB2 = MatrixUtils.kronecker(IA , muV, MB2);
	}

	/**
	 * Computes the performance measures of a PH/PH/1 model
	 * @param R1 R Matrix for the QBD process
	 * @return A string with the performance masures
	 */
	public String performanceMeasures (double[][] R1){
		// Estructures that will kee the steady state probabilities
		pis = new ArrayList<Double>( ); // Just the probability per level
		pisP = new ArrayList<DenseMatrix>( ); // The array of probabilities per level

		R = new DenseMatrix(R1);

		// Auxiliar Matrices
		DenseMatrix BON = new DenseMatrix(MB0.numRows(), N.numColumns()); //P1 = B0*N
		DenseMatrix BAL = new DenseMatrix(MB1.numRows(), MB1.numColumns()); //P2 = B1 + B0*N*B2
		DenseMatrix NORM = new DenseMatrix(MB0.numRows(), 1); //P3 = 1 + B0N*(I-R)^-1*1
		DenseMatrix I = new DenseMatrix(Matrices.identity(R.numColumns()));

		MB0.mult(N, BON); //P1 = B0*N
		BON.mult(MB2, BAL); // P2 = B0*N*B2
		BAL.add(MB1); //P2 = B1 + B0*N*B2

		// e Vector of size R
		double[][] ones2 = new double[R.numRows()][1];
		for (int i = 0; i < R.numRows(); i++) {
			ones2[i][0] = 1;
		}        
		DenseMatrix oneR = new DenseMatrix(ones2);

		// e vector of size B1
		double[][] ones1 = new double[MB0.numRows()][1];
		for (int i = 0; i < MB0.numRows(); i++) {
			ones1[i][0] = 1;
		}
		DenseMatrix oneB = new DenseMatrix(ones1);

		DenseMatrix IR = new DenseMatrix(R.copy().scale(-1)); // IR = -R
		DenseMatrix IR1 = new DenseMatrix(R.numRows(),1);
		IR.add(I); // IR = I-R
		IR.copy().solve(I, IR); // IR = (I-R)^1
		IR.mult(oneR, IR1); //P3 = (I-R)^-1*1
		BON.mult(IR1, NORM); //P3 = B0*N*(I-R)^-1*1
		NORM.add(oneB); //P3 =1 + (I-R)^-1*1

		Matrix P3T = new DenseMatrix(NORM.numColumns(), NORM.numRows());
		Matrix P2T = new DenseMatrix(BAL.numColumns(), BAL.numRows());

		BAL.copy().transpose(P2T); //P2 = (B1 + B0*N*B2)T
		NORM.copy().transpose(P3T); //P3 = (1 + B0*N*(I-R)^-1*1)T

		double[][] A = new double[P2T.numRows()][P3T.numColumns()];

		for(int i = 0; i != P2T.numRows() ; i++){
			for(int j = 0; j != P3T.numColumns() ; j++){
				if(i != P2T.numRows()-1){
					A[i][j] = P2T.get(i, j); //Balance
				}
				else{
					A[i][j] = P3T.get(0, j); //Normalization  		
				}
			}        	
		}

		DenseMatrix MF = new DenseMatrix(A);  //Balance and normalization
		DenseMatrix IF = new DenseMatrix(Matrices.identity(MF.numRows()));  
		MF.copy().solve(IF, MF); // (B-N)^1

		// 0 Vector with and 1 at the end
		double[][] oones = new double[MF.numColumns()][1];
		for (int i = 0; i < MF.numColumns(); i++) {
			if(i == MF.numColumns() - 1)
				oones[i][0] = 1;
			else
				oones[i][0] =0;
		}
		DenseMatrix oone = new DenseMatrix(oones);

		DenseMatrix PI0T = new DenseMatrix(MF.numColumns(),1);
		MF.mult(oone, PI0T);

		DenseMatrix pi0 = new DenseMatrix(1,MF.numColumns()); //pi0
		PI0T.transpose(pi0);

		DenseMatrix pi1 = new DenseMatrix(1,BON.numColumns());
		pi0.copy().mult(BON.copy(),pi1); //pi1

		double sum = 0;
		double Lq = 0;
		double piCero = 0;

		double pi = 0;
		for(int i = 0 ; i != pi0.numColumns(); i++){
			pi += pi0.get(0, i);
		}

		piCero = pi;
		sum+= pi;
		pis.add(pi);
		pisP.add(pi0);

		pi = 0;
		for(int i = 0 ; i != pi1.numColumns(); i++){
			pi += pi1.get(0, i);
		}
		sum+= pi;
		pis.add(pi); 
		pisP.add(pi1);   

		DenseMatrix pin = pi1.copy();

		for(int i = 1; 1 - sum > 1e-8; i++){
			DenseMatrix piN = new DenseMatrix(1,R.numColumns());
			pin.copy().mult(R.copy(),piN);
			pi = 0;
			for(int j = 0 ; j != pi1.numColumns(); j++){
				pi += piN.get(0, j);
			}
			pin = piN.copy();
			sum+= pi;
			Lq += pi*i;
			pis.add(pi);
			pisP.add(piN); 
		}

		double lam = getLambda();

		String resp = "";

		resp += "WIP (Queue): " + String.format("%6.4f", Lq) + "\n";
		resp += "WIP (Service): " + String.format("%6.4f", (1 - piCero)) + "\n";
		resp += "WIP (Total): " + String.format("%6.4f", (Lq + 1 - piCero)) + "\n";
		resp += "____________________________________\n";
		resp += "Expected time in Queue: " + String.format("%6.4f", Lq*lam) + "\n";
		resp += "Expected time in service: " + String.format("%6.4f", (1 - piCero)*lam) + "\n";
		resp += "Expected time in the system: " + String.format("%6.4f", (Lq + 1 - piCero)*lam) + "\n";
		resp += "____________________________________\n";

		return resp;
	}

	/**
	 * @return the level steady state probabilities
	 */
	public ArrayList<Double> getLevelSteadyStateProbs() {
		return pis;
	}

	/**
	 * @return the steady state probabilities per level
	 */
	public ArrayList<DenseMatrix> getSteadyStateProbsPerLevel() {
		return pisP;
	}

	/**
	 * Prints in console the Caracteristic matrices
	 */
	public void printMatrices(){
		printVector(new DenseMatrix (MA0), "A0");
		printVector(new DenseMatrix (MA1), "A1");
		printVector(new DenseMatrix (MA2), "A2");
		printVector(new DenseMatrix (MB0), "B0");
		printVector(new DenseMatrix (MB1), "B1");
		printVector(new DenseMatrix (MB2), "B2");
		printVector(U, "U");
		printVector(G, "G");
		printVector(R, "R");		
	}

	/**
	 * Prints a Matrix with its name
	 * @param ma: Matrix
	 * @param name: Name of the matrix
	 */
	protected void printVector(DenseMatrix matrix, String matrixName) {
		System.out.println("Matrix: " + matrixName);
		double[][] matriz = Matrices.getArray(matrix);

		for(int i = 0 ; i != matrix.numRows(); i++){
			for(int j = 0 ; j != matrix.numColumns(); j++){
				System.out.print(matriz[i][j]+"\t");
			}
			System.out.println();
		}
	}

	/**
	 * Gets the arrival rate;
	 * @return arrival rate
	 */
	private double getLambda(){
		if(lambdaP == 0){
			Matrix A = arr.getMatrix();
			int tamA = A.numColumns();
			double[] era = new double[tamA];
			for(int i = 0; i != tamA ; i++)
				era[i] = 1;

			DenseMatrix AI = new DenseMatrix(tamA, tamA);
			DenseMatrix IA = new DenseMatrix(Matrices.identity(tamA));

			A.solve(IA, AI);

			DenseVector lambda1 = negativeMatVecMult(AI, era, new DenseVector(tamA));

			this.lambdaP = rowCol(new DenseVector(arr.getVector()), lambda1);
		}
		return lambdaP;
	}

	/**
	 * Computes the product of - a Matrix and a Vector A x B
	 * @param A Vector
	 * @param era Matrix
	 * @param denseVector Vector to store the resulting matrix
	 * @return res = A x B^T
	 */
	public static DenseVector negativeMatVecMult(Matrix A, double[] era, DenseVector denseVector) {
		for (int i = 0; i != A.numRows(); i++) {
			double sum = 0;
			for (int j = 0; j < A.numColumns(); j++) {
				sum += A.get(i,j) * era[j];
			}
			denseVector.set(i, -sum);
		}
		return denseVector;
	}

	/**
	 * Computes the product of a rowvetor with a column vector
	 * @param fila: row vector
	 * @param colu: column vector
	 * @return The value of the multiplication
	 */
	private double rowCol(DenseVector fila, DenseVector colu) {
		double[] row = fila.getData();
		double[] col = colu.getData();
		if (row.length != col.length)
			throw new IndexOutOfBoundsException(
					"row.length != col.length " + 
					row.length +" != " + col.length);
		double ret = 0;
		for(int i  = 0; i != row.length ; i++){
			ret += row[i]*col[i];
		}
		return ret;
	}

	/**
	 * Computes the product of two vectors A x B^T
	 * @param A Vector
	 * @param B Vector
	 * @param res Vector to store the resulting matrix
	 * @return res = A x B^T
	 */
	public static Matrix multVector(DenseVector A,
			DenseVector B, DenseMatrix res) {
		int n = res.numRows();
		int m = res.numColumns();

		if (A.size() != n)
			throw new IndexOutOfBoundsException("A.size() != res.numRows() ("
					+ A.size() + " != " + n + ")");
		if (B.size() != m)
			throw new IndexOutOfBoundsException(
					"B.size() != res.numColumns() (" + B.size() + " != " + m
					+ ")");
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				res.set(i, j, A.get(i) * B.get(j));
			}
		}
		return res;
	}

	/**
	 * Creates a dense matrix from a dense vector
	 * @param dv: dense vector
	 * @return the dense matrix
	 */
	private static DenseMatrix rowVectoraMat(DenseVector dv){
		double[] data = dv.getData();
		double[][] mat = new double[1][data.length];
		for(int i = 0; i != data.length; i++){
			mat[0][i] = data[i];
		}

		return new DenseMatrix(mat);
	}
	
	/**
	 * Verifies if the system is stable, this means that the arrival rate
	 * is lower than the service rate
	 * @return True if is stable, False if not
	 */
	public boolean unstableSystem(){
		if(arr == null || serv == null)
			return true;
		return arr.expectedValue()<=serv.expectedValue();
	}

	/**
	 * Sets the arrival rate
	 * @param lambda arrival rate
	 */
	public void setLambda(double lambda) {
		// TODO Auto-generated method stub
		this.lambdaP = lambda;
	}

	/**
	 * Computs the R matrix for the QBD process
	 * @return the R Matrix
	 */
	public abstract double[][] getRmatrix();
}
