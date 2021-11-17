package jmarkov;

import java.io.PrintWriter;
import java.util.SortedSet;
import java.util.TreeSet;

import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.NotUnichainException;
import jmarkov.solvers.GeometricSolver;
import jmarkov.solvers.MtjLogRedSolver;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.VectorEntry;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;

/**
 * The class GeomProcess represents a continuous or discrete Quasi
 * Birth and Death process. This class extends the class
 * SimpleMarkovProcess. The class generate the G matrix through the
 * Logarithmic Reduction algorithm. The user should extend State to
 * generate Sub-States.
 * @param <Sub> The sub-states.
 * @param <E> The events class
 * @author Julio Goez.
 * @version 1.0
 */
public abstract class GeomProcess<Sub extends State, E extends Event> extends
        SimpleMarkovProcess<GeomState<Sub>, E> {

    /** Integer to identify the final State of rLevel i. */
    private int level1Idx = -1;
    private int level2Idx = -1;
    private int boundaryIdx = -1;

    /** Current Solver */
    protected GeometricSolver GeometrixSolver = null;
    /** Default Solver */
    protected GeometricSolver defaultGeometrixSolver = null;

    /**
     * Builds a GeomProcess
     * @param i0 Initial state. MUST be a boundary state!
     * @param eSet the event set.
     */

    public GeomProcess(Sub i0, EventsSet<E> eSet) {
        super(new GeomState<Sub>(i0, 0), eSet);
        addMOP("Expected Level");
    }

    /**
     * @return an array with all the states in levels: 0 (boundary),
     *         1, and 2.
     */
    @Override
    public StatesSet<GeomState<Sub>> getStates() {
        StatesSet<GeomState<Sub>> states = super.getStates();
        getLevelsIndices();
        return states;
    }

    /**
     * Computes the max index for levels 0,1 and 2.
     * @return The max index for levels 0,1 and 2.
     */
    private int[] getLevelsIndices() {
        GeomState<Sub>[] theGeomStates = null;
        /*
         * It initializes the theGeomStates with the States fully
         * analyzed and stored in objs and initializes the integers
         * level1, level2 and boundary that are identifying the last
         * State in every rLevel.
         */
        theGeomStates = super.getStates().toStateArray();
        int maxIdx[] = new int[3];
        for (GeomState<Sub> s : theGeomStates) {
            int level = s.getLevel();
            maxIdx[level] = s.getIndex();
        }
        boundaryIdx = maxIdx[0];
        level1Idx = maxIdx[1];
        level2Idx = maxIdx[2];
        return maxIdx;
    }

    // TODO: theBoundaryStates y theTypicalStates deben ser
    // construidos al tiempo
    private Sub[] theBoundaryStates = null;

    /**
     * Returns an array with the States in the boundary level.
     * @return an array with the Sub-States
     */
    public Sub[] getBoundaryStates() {
        if (theBoundaryStates == null) {
            GeomState<Sub>[] theStates = getStates().toStateArray();
            StatesSet<Sub> bdrySet = new StatesSet<Sub>();
            for (int i = 0; i <= boundaryIdx; i++) {
                bdrySet.add(theStates[i].getSubState());
            }
            theBoundaryStates = bdrySet.toStateArray();
        }
        return theBoundaryStates;
    }

    private Sub[] theTypicalStates = null;

    /**
     * Returns an array with the States in the typical levels.
     * @return an array with the Sub-States
     */
    public Sub[] getTypicalStates() {
        if (theTypicalStates == null) {
            GeomState<Sub>[] theStates = getStates().toStateArray();
            StatesSet<Sub> typicStates = new StatesSet<Sub>();
            for (int i = boundaryIdx + 1; i <= level1Idx; i++) {
                typicStates.add(theStates[i].getSubState());
            }
            theTypicalStates = typicStates.toStateArray();
        }
        return theTypicalStates;
    }

    /**
     * The Number of States in the boundary level.
     * @return the Number of States in the boundary level.
     */
    public int getNumBoundaryStates() {
        if (!isGenerated())
            generate();
        getLevelsIndices();
        return (boundaryIdx + 1);
    }

    /**
     * The number of states in the typical levels.
     * @return the number of states in the typical levels.
     */
    public int getNumTypicalStates() {
        if (!isGenerated())
            generate();
        getLevelsIndices();
        return (level1Idx - boundaryIdx);

    }

    /**
     * This method constructs any A_n or B_ij matrix existing in the
     * process and that are necessary to calculated R matrix.
     * @param rowMax The upper limit row index.
     * @param rowMin The lower limit row index.
     * @param colMax The upper limit column index.
     * @param colMin The lower limit column index.
     * @return any A matrix in the process.
     */
    public Matrix getSubMatrices(int rowMin, int rowMax, int colMin, int colMax) {
        int[] rows = new int[rowMax - rowMin + 1];
        int[] rowsSub = new int[rowMax - rowMin + 1];
        int[] columns = new int[colMax - colMin + 1];
        int[] columnsSub = new int[colMax - colMin + 1];
        //long startTimeIdx = System.currentTimeMillis();
    	
        for (int i = rowMin; i < rowMax + 1; i++) {
            rows[i - rowMin] = i;
            rowsSub[i - rowMin] = i - rowMin;
        }

        for (int i = colMin; i < colMax + 1; i++) {
            columns[i - colMin] = i;
            columnsSub[i - colMin] = i - colMin;
        }

        /*long stopTimeIdx = System.currentTimeMillis();
        long elapsedTimeIdx = stopTimeIdx - startTimeIdx;
        System.out.println("getSubMatrices - Index building time: "+elapsedTimeIdx+" ms");*/
        
        //long startTimeGetGen = System.currentTimeMillis();
        Matrix generator = this.getMtjGenerator();
        /*long stopTimeGetGen = System.currentTimeMillis();
        long elapsedTimeGetGen = stopTimeGetGen - startTimeGetGen;
        System.out.println("getSubMatrices - getMtjGeneartor: "+elapsedTimeGetGen+" ms");*/
        
        
        // int rDim = rows.length;
        // int cDim = columns.length;
        // double[][] tempM = new double[rDim][cDim];

        // generator.get(rows, columns, tempM);
        // FlexCompRowMatrix SubMatrix = new FlexCompRowMatrix(rDim,
        // cDim);
        // SubMatrix.set(rowsSub, columnsSub, tempM);
        //long startTimeGetSubMat = System.currentTimeMillis();
        
        Matrix SubMatrix = Matrices.getSubMatrix(generator, rows, columns);
        /*long stopTimeGetSubMat = System.currentTimeMillis();
        long elapsedTimeGetSubMat = stopTimeGetSubMat - startTimeGetSubMat;
        System.out.println("getSubMatrices - Matrices.getSubmatrix: "+elapsedTimeGetSubMat+" ms");*/
        return SubMatrix;
    }

    private Matrix R = null;

    /**
     * The R Matrix of the Geometric solution. This matrix solves
     * <p>
     * <B>A<sub>0</sub>+RA<sub>1</sub> + R<sup>2</sup>A<sub>2</sub> =
     * 0</B>
     * @return The <b>R</b> Matrix of the Geometric solution. If the
     *         system is not stable it returns a zero matrix.
     * @throws NotUnichainException
     */

    public Matrix getRmatrix() throws NotUnichainException {

        if (R == null) {
            int n = getNumTypicalStates();
            if (!stabilityChecked)
                isStable = isStable();
            if (!isStable)
                return new DenseMatrix(n, n);
            //long startTime = System.currentTimeMillis();
            R = new DenseMatrix(getGeometrixSolver().getRmatrix());
	    	/*long stopTime = System.currentTimeMillis();
	        long elapsedTime = stopTime - startTime;
	        System.out.println("Compute R exec time: "+elapsedTime+" ms");*/
        }
        return R;
    }

    /**
     * @return Matrix R in an array of doubles.
     * @throws NotUnichainException
     */
    public double[][] matrixRtoArray() throws NotUnichainException {
        if (R == null)
            R = getRmatrix();
        int row;
        int column;
        int n = R.numRows();
        double[][] matrixR = new double[n][n];
        for (MatrixEntry e : R) {
            row = e.row();
            column = e.column();
            matrixR[row][column] = e.get();
        }
        return matrixR;
    }

    private Matrix[] A; // The A Matrices computed

    /**
     * Returns the matrices of the repeating levels, A0, A1, and A2.
     * If the model has not been generated it will be.
     * @return [A0, A1, A2]
     */
    public Matrix[] getAMatrices() {
        if (A == null) {
        	//long startTimeGetIdx = System.currentTimeMillis();
        	getLevelsIndices();
            //long stopTimeGetIdx = System.currentTimeMillis();
            //long elapsedTimeGetIdx = stopTimeGetIdx - startTimeGetIdx;
	        //System.out.println("getMatrices - GetLevelIndices time: "+elapsedTimeGetIdx+" ms");
            
	        //long startTimeGetSubMat = System.currentTimeMillis();
            Matrix A0 = getSubMatrices((boundaryIdx + 1), level1Idx,
                    (level1Idx + 1), level2Idx);
            Matrix A1 = getSubMatrices((boundaryIdx + 1), level1Idx,
                    (boundaryIdx + 1), level1Idx);
            Matrix A2 = getSubMatrices((level1Idx + 1), level2Idx,
                    (boundaryIdx + 1), level1Idx);
            //long stopTimeGetSubMat = System.currentTimeMillis();
            //long elapsedTimeGetSubMat = stopTimeGetSubMat - startTimeGetSubMat;
	        //System.out.println("getMatrices - GetSubMatrices time: "+elapsedTimeGetSubMat+" ms");
	        
	        
            A = new Matrix[] { A0, A1, A2 };
            debug(3, "A0 is:\n" + A0.toString());
            debug(3, "A1 is:\n" + A1.toString());
            debug(3, "A2 is:\n" + A2.toString());
        }
        return A;
    }

    private Matrix[] B = null;

    /**
     * Returns the matrices B00, B01 and B10. It causes the generation
     * of the model if it has not been generated.
     * @return an array with {B00,B01,B10} in that order.
     */
    public Matrix[] getBMatrices() {
        if (B == null) {
            getLevelsIndices();
            Matrix B00 = getSubMatrices(0, boundaryIdx, 0, boundaryIdx);
            Matrix B01 = getSubMatrices(0, boundaryIdx, (boundaryIdx + 1),
                    level1Idx);
            Matrix B10 = getSubMatrices((boundaryIdx + 1), level1Idx, 0,
                    boundaryIdx);
            B = new Matrix[] { B00, B01, B10 };
            debug(3, "B00 is:\n" + B00.toString());
            debug(3, "B01 is:\n" + B01.toString());
            debug(3, "B10 is:\n" + B10.toString());
        }
        return B;
    }

    private double[] pis = null;

    /**
     * Computes and returns the initial solution [ pi(0), pi(1) ].
     * @return an array with the initial solution [ pi(0), pi(1) ]
     * @throws NotUnichainException
     */
    public double[] getInitialSol() throws NotUnichainException {

        if (pis == null) {
            /** Stability test */
        	//long startTimeStable = System.currentTimeMillis();
            boolean stable = isStable();
            /*long stopTimeStable = System.currentTimeMillis();
	        long elapsedTimeStable = stopTimeStable - startTimeStable;
	        System.out.println("\nisStable exec time: "+elapsedTimeStable+" ms\n");*/
	        
            if (!stable) {
                debug(0, "The system is unstable");
                return new double[getNumBoundaryStates()
                        + getNumTypicalStates()];
            }
            getLevelsIndices();
            Matrix A[] = getAMatrices();
            Matrix A1 = A[1], A2 = A[2];
            Matrix B[] = getBMatrices();
            Matrix B00 = new FlexCompRowMatrix(B[0], true), B01 = B[1], B10 = new FlexCompRowMatrix(
                    B[2], true);

            for (int i = 0; i < B00.numRows(); i++) {
                B00.set(i, 0, 1);
            }

            //long startTimeR = System.currentTimeMillis();
            Matrix R = getRmatrix();
            /*long stopTimeR = System.currentTimeMillis();
	        long elapsedTimeR = stopTimeR - startTimeR;
	        System.out.println("Compute R exec time (from getInitialSol() ): "+elapsedTimeR+" ms");*/
	        
	        /*long startTimePi0 = System.currentTimeMillis();
            long startTime2 = System.currentTimeMillis();*/
	        DenseVector e = new DenseVector(R.numRows());
            
            //for (int i = 0; i < e.numRows(); i++) {
	        for (VectorEntry i : e){
	        	i.set(1);
            }
            Matrix I = Matrices.identity(R.numRows());
            I.add(-1, R);
            DenseVector x = new DenseVector(R.numRows());
            I.solve(e, x);
            for (int i = 0; i < B10.numRows(); i++) {
                B10.set(i, 0, x.get(i));
            }
            /*
            long stopTime2 = System.currentTimeMillis();
	        long elapsedTime2 = stopTime2 - startTime2;
	        System.out.println("Time NEW pi solve I-R: "+elapsedTime2+" ms"); */
	        
            

            int columDimension = B00.numColumns() + B01.numColumns();
            int rowDimension = B00.numRows() + B10.numRows();
            //long startTime3 = System.currentTimeMillis();
            R.multAdd(A2, A1);
            /*long stopTime3 = System.currentTimeMillis();
            long elapsedTime3 = stopTime3 - startTime3;
	        System.out.println("Time pi multAdd: "+elapsedTime3+" ms");*/
	        
            DenseVector Pis = new DenseVector(rowDimension);
            //DenseVector Pis2 = new DenseVector(rowDimension);
            DenseVector zeros = new DenseVector(columDimension);
            zeros.zero();
            zeros.set(0, 1);

            //startTime3 = System.currentTimeMillis();
            DenseMatrix MTotal = assmbleMatrix(B00, B01, B10, A1);
            /*stopTime3 = System.currentTimeMillis();
            elapsedTime3 = stopTime3 - startTime3;
	        System.out.println("Time pi assemble: "+elapsedTime3+" ms");*/
	        
            //startTime3 = System.currentTimeMillis();
            MTotal.transSolve(zeros, Pis);
            /*stopTime3 = System.currentTimeMillis();
	        elapsedTime3 = stopTime3 - startTime3;
	        System.out.println("Time pi solve pi0: "+elapsedTime3+" ms");*/
	        
            /*long stopTimePi0 = System.currentTimeMillis();
	        long elapsedTimePi0 = stopTimePi0 - startTimePi0;
	        System.out.println("Compute Pi0 exec time: "+elapsedTimePi0+" ms");*/
            pis = Pis.getData();
        }
        return pis;
    }

    private DenseMatrix assmbleMatrix(Matrix B00, Matrix B01, Matrix B10,
            Matrix B11) {
        int columDimension = B00.numColumns() + B01.numColumns();
        int rowDimension = B00.numRows() + B10.numRows();
        DenseMatrix MTotal = new DenseMatrix(rowDimension, columDimension);

        setSubMatrix(MTotal, B00, 0, 0);
        setSubMatrix(MTotal, B01, 0, B00.numColumns());
        setSubMatrix(MTotal, B10, B00.numRows(), 0);
        setSubMatrix(MTotal, B11, B00.numRows(), B00.numColumns());
        return MTotal;
    }

    private void setSubMatrix(Matrix mat, Matrix submatrix, int rowOffset,
            int colOffset) {
        for (MatrixEntry entry : submatrix) {
            mat.set(entry.row() + rowOffset, entry.column() + colOffset, entry
                    .get());
        }
    }

    // private Matrix assmbleMatrix(Matrix B00, Matrix B01, Matrix
    // B10, Matrix M22){
    // int columDimension = B00.numColumns() + B01.numColumns();
    // int rowDimension = B00.numRows() + B10.numRows();
    // DenseMatrix MTotal = new DenseMatrix(rowDimension,
    // columDimension);
    // double[][] Mtotal11 = Matrices.getArray(B00);
    // double[][] Mtotal21 = Matrices.getArray(B10);
    // double[][] Mtotal12 = Matrices.getArray(B01);
    // double[][] Mtotal22 = Matrices.getArray(M22);
    //
    // int[] row1 = new int[B00.numRows()];
    // int[] row2 = new int[B10.numRows()];
    // int[] col1 = new int[B00.numColumns()];
    // int[] col2 = new int[B01.numColumns()];
    //
    // for (int i = 0; i < row1.length; i++) {
    // row1[i] = i;
    // }
    //
    // for (int i = 0; i < row2.length; i++) {
    // row2[i] = i + row1.length;
    // }
    //
    // for (int i = 0; i < col1.length; i++) {
    // col1[i] = i;
    // }
    //
    // for (int i = 0; i < col2.length; i++) {
    // col2[i] = i + col1.length;
    // }
    //
    // MTotal.set(row1, col1, Mtotal11);
    // MTotal.set(row1, col2, Mtotal12);
    // MTotal.set(row2, col1, Mtotal21);
    // MTotal.set(row2, col2, Mtotal22);
    // return Mtotal;
    // }

    /** Used for pi(1) */
    private double[] pi1 = null;

    /**
     * Returns the steady State probabilities for level 1.
     * @return the array pi(1)
     * @throws NotUnichainException
     */
    public double[] getVectorPi1() throws NotUnichainException {
        if (pi1 == null) {
            pis = getInitialSol();
            int n0 = getNumBoundaryStates();
            int n = getNumTypicalStates();
            pi1 = new double[n];
            System.arraycopy(pis, n0, pi1, 0, n);
        }
        return pi1;
    }

    /** Used for pi(1) */
    private double[] pi1mod = null;

    /**
     * Returns the steady State probabilities for level 1.
     * @return the array pi(1)(I-R)^(-1)
     * @throws NotUnichainException
     */
    public double[] getVectorPi1Mod() throws NotUnichainException {
        if (pi1mod == null) {
        	//long startTimePi1 = System.currentTimeMillis();
            pi1 = getVectorPi1();
            /*long stopTimePi1 = System.currentTimeMillis();
            long elapsedTimePi1 = stopTimePi1 - startTimePi1;
	        System.out.println("Compute Pi1 exec time: "+elapsedTimePi1+" ms");*/
	        
            int n = getNumTypicalStates();
            R = getRmatrix();
            
            //long startTimePi1Mod = System.currentTimeMillis();
            DenseVector sol = new DenseVector(n);
            Matrix ImR = new DenseMatrix(Matrices.identity(n));
            ImR.add(-1, R);
            ImR.transSolve(new DenseVector(pi1), sol);
            pi1mod = sol.getData();
            /*long stopTimePi1Mod = System.currentTimeMillis();
	        long elapsedTimePi1Mod = stopTimePi1Mod - startTimePi1Mod;
	        System.out.println("Compute Pi1Mod exec time: "+elapsedTimePi1Mod+" ms");*/

        }
        return pi1mod;
    }

    /** To keep the computed pi0 */
    private double pi0[] = null;

    /**
     * Returns the steady State probabilities for boundary level.
     * @return the array pi(0)
     * @throws NotUnichainException
     */
    public double[] getVectorPi0() throws NotUnichainException {
        if (pi0 == null) {
            pis = getInitialSol();
            pi0 = new double[boundaryIdx + 1];
            for (int i = 0; i < boundaryIdx + 1; i++) {
                pi0[i] = pis[i];
            }
        }
        return pi0;
    }

    /**
     * Return an array with the probabilities for the given level.
     * @param level
     * @return probabilities array pi(k).
     * @throws NotUnichainException
     */
    public double[] getSteadyState(int level) throws NotUnichainException {

        DenseVector piLevel;// pi(k)
        if (level == 0) {
            piLevel = new DenseVector(getVectorPi0());
        } else {
            piLevel = new DenseVector(getVectorPi1());
            Matrix R = getRmatrix();
            int size = piLevel.size();
            for (int i = 1; i < level; i++) {
                DenseVector newPiLevel = new DenseVector(size);
                R.transMult(piLevel, newPiLevel); // pi(k) =
                // pi(k-1)* R
                piLevel = newPiLevel;
            }
        }
        return piLevel.getData();
    }

    /**
     * Computes the steady state probabilities for the generated
     * States (up to level 2).
     * @return (pi(0), pi(1), pi(2)).
     * @throws NotUnichainException
     */
    public double[] steadyProbabilities() throws NotUnichainException {
        int N = getNumStates();
        double steadyProb[] = new double[N];
        if (isStable()) {
            double initialPis[] = getInitialSol();
            int n0 = initialPis.length;
            System.arraycopy(initialPis, 0, steadyProb, 0, n0);

            double pi2[] = getSteadyState(2);
            int n2 = pi2.length;
            System.arraycopy(pi2, 0, steadyProb, n0, n2);
        }
        return steadyProb;

    }

    /*
     * L = <b>1pi + 2 pi(2) + 3 pi(3) + ... = pi<sub>1</sub>(I-R)<sup>-2</sup>1
     * </B>
     */

    /**
     * Returns the Expected Value for the Level. That is
     * <P>
     * <p class=MsoNormal>
     * <span lang=ES-CO style='font-family:Times;mso-ansi-language:
     * ES-CO;mso-no-proof:yes'>L </span><span lang=ES-CO
     * style='font-family:Symbol;
     * mso-ansi-language:ES-CO;mso-no-proof:yes'>=<b
     * style='mso-bidi-font-weight:normal'> (p</b><sub>1 </sub>+ 2<b
     * style='mso-bidi-font-weight:normal'>p</b><sub>2</sub> +3<b
     * style='mso-bidi-font-weight:normal'>p</b><sub>3</sub> + ... )<b
     * style='mso-bidi-font-weight:normal'>1</b>= <b
     * style='mso-bidi-font-weight:normal'>p</b><sub>1</sub></span><span
     * lang=ES-CO
     * style='font-family:Times;mso-ansi-language:ES-CO;mso-no-proof:yes'>(<b
     * style='mso-bidi-font-weight:normal'>I-R</b>)</span><sup><span
     * lang=ES-CO
     * style='font-family:Symbol;mso-ansi-language:ES-CO;mso-no-proof:yes'>-2</span></sup><b
     * style='mso-bidi-font-weight:normal'><span lang=ES-CO
     * style='font-family:Symbol;
     * mso-ansi-language:ES-CO;mso-no-proof:yes'>1</span></b><span
     * lang=ES-CO
     * style='font-family:Symbol;mso-ansi-language:ES-CO;mso-no-proof:yes'><o:p></o:p></span>,
     * where <b style='mso-bidi-font-weight:normal'><span lang=ES-CO
     * style='font-family:Symbol;
     * mso-ansi-language:ES-CO;mso-no-proof:yes'>1</span></b> is a
     * column vector of ones.
     * </p>
     * @return L.
     * @throws NotUnichainException
     */
    public double getExpectedLevel() throws NotUnichainException {
        double expLevel = -1;
        if (expLevel == -1) {
        	//long startTimeEL = System.currentTimeMillis();
        	
            double pi1mod[] = getVectorPi1Mod();
            int n = getNumTypicalStates();
            
            //long startTimeR = System.currentTimeMillis();
            R = getRmatrix();
            /*long stopTimeR = System.currentTimeMillis();
            long elapsedTimeR = stopTimeR - startTimeR;
	        System.out.println("Compute R exec time (from getExpectedLevel): "+elapsedTimeR+" ms");*/
	        
            DenseVector sol = new DenseVector(n);
            Matrix ImR = new DenseMatrix(Matrices.identity(n));
            ImR.add(-1, R);
            ImR.transSolve(new DenseVector(pi1mod), sol);
            expLevel = sol.norm(Vector.Norm.One);
            
            /*long stopTimeEL = System.currentTimeMillis();
            long elapsedTimeEL = stopTimeEL - startTimeEL;
	        System.out.println("Compute Exp Level exec time: "+elapsedTimeEL+" ms");*/

        }
        return expLevel;
    }

    private double getExpectedLevelOld() throws NotUnichainException {
        double result = 0.0;
        Matrix R = getRmatrix();
        DenseVector Pi1 = new DenseVector(getVectorPi1());

        DenseMatrix TMP = new DenseMatrix(Matrices.identity(R.numRows()));
        DenseVector Vtmp = new DenseVector(R.numRows());
        DenseMatrix Square = new DenseMatrix(R.numRows(), R.numRows());
        DenseMatrix Inv = new DenseMatrix(R.numRows(), R.numRows());

        TMP.add(-1, R);
        TMP.mult(TMP, Square);
        Square.solve(Matrices.identity(Square.numRows()), Inv);
        TMP = new DenseMatrix(Matrices.identity(R.numRows()));
        R.mult(Square, TMP);
        TMP.transpose();
        TMP.mult(Pi1, Vtmp);

        for (int i = 0; i < R.numRows(); i++) {
            result = result + Vtmp.get(i);
        }
        return result;
    }

    private boolean stabilityChecked = false;
    private boolean isStable = false;

    /**
     * Determines if the system is stable.
     * @return true if the system is stable.
     */
    public boolean isStable() {
        if (!stabilityChecked) {
        	Matrix A[] = getAMatrices();
            Matrix A0 = A[0], A1 = A[1], A2 = A[2];

            int dim = A0.numRows();
            DenseMatrix sum = new DenseMatrix(dim, dim);
            DenseMatrix b = new DenseMatrix(dim, 1);
            DenseMatrix f = new DenseMatrix(dim, 1);
            DenseMatrix ft = new DenseMatrix(1, dim);
            DenseMatrix right = new DenseMatrix(1, dim);
            DenseMatrix left = new DenseMatrix(1, dim);
            double R = 0;
            double L = 0;

            b.zero();
            b.set((dim - 1), 0, 1);
            
            sum.add(A0);
            sum.add(A1);
            sum.add(A2);
	        
            for (int i = 0; i < dim; i++) {
                sum.set(i, (dim - 1), 1);
            }
            sum.transSolve(b, f);
	        
            f.transpose(ft);

            ft.mult(A2, right);
            ft.mult(A0, left);

            for (int i = 0; i < dim; i++) {
                R += right.get(0, i);
                L += left.get(0, i);
            }
            
            isStable = (L < R);
            stabilityChecked = true;
        }
        return isStable;
    }

    @Override
    public void printAll(PrintWriter out) {
        int N = getNumStates();
        out.println(description() + "\n");
        out.println("System has: \n" + //
                getNumBoundaryStates() + " Boundary States.\n"
                + getNumTypicalStates() + " States in  typical levels.\n");
        boolean isStable = isStable();
        if (!isStable)
            out.println("The system is NOT stable !");
        if (N < 100) {
            printStates(out);
            printDenseMatrix(out, 10, 4, false, false, getLevelsIndices());
        }
        out.println();
        try {
            Matrix A[] = getAMatrices();
            Matrix A0 = A[0], A1 = A[1], A2 = A[2];
            Matrix R = getRmatrix();
            out.println();
            printDenseAMatrix(out, A0, getTypicalStates(), "A0");
            out.println();
            printDenseAMatrix(out, A1, getTypicalStates(), "A1");
            out.println();
            printDenseAMatrix(out, A2, getTypicalStates(), "A2");
            out.println();

            if (isStable) {
                printDenseAMatrix(out, R, getTypicalStates(), "R");
                out.println();
                printMOPs(out);
                printEventsRates(out);
            }
        } catch (NotUnichainException e) {
            out.print(e);
        }
    }

    @Override
    public void printStates(PrintWriter out, int width, int probDecimals) {
        State[] stts = getStates().toStateArray();
        debug(1, "Printing States...");
        int n = stts.length;
        int maxL = 9; // Lets find max label width

        try {

            double[] pi = steadyProbabilities();

            for (int i = 0; i < n; i++) {
                maxL = Math.max(maxL, stts[i].label().length());
            }
            int w = maxL + 3;
            int w2 = width + 3;
            // Print headers
            out.println(pad("", w) + pad("EQUILIBRUM", w2, false));
            out.println(pad("STATE", w, false) + pad("PROBAB.", w2, false)
                    + "DESCRIPTION");
            out.println();
            String hLine = hLine(w + w2 + w + w2 + 15);

            int idx[] = getLevelsIndices();
            int k = 0;
            out.println(hLine);
            for (int i = 0; i < n; i++) {
                out.println(pad(stts[i].label(), w, false)
                        + pad(pi[i], w, probDecimals, false)
                        + ((stts[i].description() != "") ? ""
                                + stts[i].description() : ""));
                if (i == idx[k]) {
                    out.println(hLine);
                    k++;
                }
            }

            for (int i = 0; i < 3; i++) {
                out.println(pad(".", w, false) + pad(".", w, false)
                        + pad(".", w, false));
            }
        } catch (NotUnichainException e) {
            out.print(e);
        }

        out.println();
    }

    private void printDenseAMatrix(PrintWriter out, Matrix A, Sub[] states,
            String Matrix) {
        printDenseAMatrix(out, 10, 4, false, A, states, Matrix);
    }

    private void printDenseAMatrix(PrintWriter out, int width,
            int rateDecimals, boolean printZeros, Matrix A, Sub[] subStates,
            String matrixName) {
        // setStatus(Status.WRITING);
        int n = subStates.length;
        int maxL = 0; // Lets find max label width
        for (int i = 0; i < n; i++) {
            maxL = Math.max(maxL, subStates[i].label().length());
        }
        int w = Math.max(maxL + 1, width);
        out.println(matrixName + " MATRIX: ");
        out.print(pad(" ", w));
        for (int i = 0; i < n; i++)
            out.print(pad(subStates[i].toString(), w));
        for (cnt = 0; cnt < n; cnt++) {
            int i = cnt;
            out.println();
            out.print(pad(subStates[i].toString(), w));
            for (int j = 0; j < n; j++) {
                if (A.get(i, j) == 0 && !printZeros) {
                    out.print(pad("", w));
                } else {
                    out.print(pad(A.get(i, j), w, rateDecimals));
                }
            }
        }
        out.flush();
    }

    private void printDenseRMatrix(PrintWriter out, Matrix R, String Matrix) {
        printDenseRMatrix(out, 10, 3, false, R, Matrix);
    }

    private void printDenseRMatrix(PrintWriter out, int width,
            int rateDecimals, boolean printZeros, Matrix R, String Matrix) {
        // setStatus(WRITING);
        int n = R.numColumns();
        int maxL = 0; // Lets find max label width
        // for (int i = 0; i < n; i++) {
        // maxL = Math.max(maxL, stts[i].label().length());
        // }
        int w = Math.max(maxL + 1, width);
        out.println(Matrix + " MATRIX: ");
        // double[][] theQ = getGenerator();
        out.print(pad(" ", w));
        // for (int i = 0; i < n; i++)
        // out.print(pad(stts[i].toString(), w));
        for (cnt = 0; cnt < n; cnt++) {
            int i = cnt;
            out.println();
            // out.print(pad(stts[i].toString(), w));
            for (int j = 0; j < n; j++) {
                if (R.get(i, j) == 0 && !printZeros) {
                    out.print(pad("", w));
                } else {
                    out.print(pad(R.get(i, j), w, rateDecimals));
                }
            }
        }
        out.flush();
    }

    /**
     * Determines the destination set of States when events e occurs.
     * It has to be implemented by the subclass.
     * @param i current State.
     * @param iLevel absolute level of current State. For QBD this is
     *        0, 1 or 2. Anything above 2 should report the same
     *        result.
     * @param e The Event that ocurred.
     * @return The destination States
     */

    public abstract GeomRelState[] dests(Sub i, int iLevel, E e);

    /**
     * Overrides SimpleMarkovProcess' method
     * @param i current state
     * @param e Event
     * @return the destinations states
     */
    @Override
    public final States<GeomState<Sub>> dests(GeomState<Sub> i, E e) {
        GeomRelState[] dests = dests(i.getSubState(), i.getLevel(), e);
        SortedSet<GeomState<Sub>> newDest = new TreeSet<GeomState<Sub>>();
        for (GeomRelState<Sub> s : dests) {
            GeomState<Sub> j;
            if (s.isBoundary()) {
                j = new GeomState<Sub>(s.getSubState(), 0);
                newDest.add(j);
            } else {
                int newLevel = i.getLevel() + s.getRelLevel();
                if (newLevel <= 2) {
                    j = new GeomState<Sub>(s.getSubState(), newLevel);
                    newDest.add(j);
                }
            }
        }
        int size = newDest.size();
        if (size == 0)
            return null;
        // GeomState<Sub> elem = newDest.first();
        // GeomState<Sub>[] result = (GeomState<Sub>[])
        // java.lang.reflect.Array
        // .newInstance(elem.getClass(), size);
        return new StatesSet<GeomState<Sub>>(newDest);
    }

    /**
     * The user must extend this method to determine which events are
     * active.
     * @param substate the current sub state
     * @param iLevel Absolute level of current State i. You should
     *        test only whether it is 0 (boundary), 1 or greater than
     *        1. Your code should not behave any different if the
     *        level is 2, or 3, etc
     * @param e The event being tested.
     * @return tru if this event occurs
     */
    public abstract boolean active(Sub substate, int iLevel, E e);

    /**
     * The user cannot extend this method. GeomProcess detemines this
     * based on <code>active</code>
     * @see #active(State, int, Event)
     * @see jmarkov.SimpleMarkovProcess#active(State, Event)
     */
    @Override
    public final boolean active(GeomState<Sub> i, E e) {
        return active(i.getSubState(), i.getLevel(), e);
    }

    /**
     * @param i current sub state
     * @param ilevel current state's absolute level
     * @param j destination sub state
     * @param jLevel destination level
     * @param e Event
     * @return rate of occurrance
     */
    public abstract double rate(Sub i, int ilevel, Sub j, int jLevel, E e);

    /*
     * (non-Javadoc)
     * @see jmarkov.SimpleMarkovProcess#rate(Stte, Stte, Evt)
     */
    @Override
    public double rate(GeomState<Sub> i, GeomState<Sub> j, E e) {
        Sub iSub = i.getSubState();
        Sub jSub = j.getSubState();

        return rate(iSub, i.getLevel(), jSub, j.getLevel(), e);
    }

    /*
     * (non-Javadoc)
     * @see jmarkov.SimpleMarkovProcess#reset()
     */
    @Override
    public void reset() {
        super.reset();
        this.A = null;
        this.B = null;
        this.boundaryIdx = -1;
        this.level1Idx = -1;
        this.level2Idx = -1;
        this.pi0 = null;
        this.pi1 = null;
        this.pis = null;
        this.R = null;
        this.pi1mod = null;
        this.stabilityChecked = false;
    }

    /**
     * @throws NotUnichainException
     * @see jmarkov.SimpleMarkovProcess#getMOPsMoment(int, int)
     */
    @Override
    public double getMOPsMoment(int mopNum, int m) throws NotUnichainException {
        if (getMOPNames(mopNum) == "Expected Level") {
            if (m == 1)
                return getExpectedLevel();
            else
                return Double.NaN;
        }
        double[] pi0 = getVectorPi0();
        double[] pi1m = getVectorPi1Mod();
        String mopName = getMOPNames(mopNum);
        debug(1, "Computing moment " + m + " for MOP " + mopName);
        int N0 = getNumBoundaryStates();
        int N1 = getNumTypicalStates();
        double sum = 0;
        Sub theBStates[] = getBoundaryStates();
        for (int i = 0; i < N0; i++) {
            sum += pi0[i] * Math.pow(theBStates[i].getMOP(mopNum), m);
        }
        Sub theTStates[] = getTypicalStates();
        // sum = pi(0) * r0, so far.
        for (int i = 0; i < N1; i++) {
            sum += pi1m[i] * Math.pow(theTStates[i].getMOP(mopNum), m);
        }
        return sum;
    }

    /**
     * @throws NotUnichainException
     * @see jmarkov.SimpleMarkovProcess#getEventRate(int)
     */
    @Override
    public double getEventRate(int eNum) throws NotUnichainException {
        Sub i;
        E e;
        e = theEvents[eNum];
        debug(1, "Computing Events Rate for Event " + e);
        double sumRate = 0.0;
        double pi0[] = getVectorPi0();
        Sub bdryStates[] = getBoundaryStates();
        int N0 = pi0.length;
        for (int s = 0; s < N0; s++) {
            i = bdryStates[s];
            if (active(i, 0, e)) {
                for (GeomRelState<Sub> j : dests(i, 0, e)) {
                    sumRate += pi0[s] * rate(i, 0, j.subState, j.rLevel, e);
                }
            }
        }
        // so far sum is rates for pi(0)
        double pi1m[] = getVectorPi1Mod();
        Sub typStates[] = getTypicalStates();
        int N1 = pi1m.length;
        for (int s = 0; s < N1; s++) {
            i = typStates[s];
            if (active(i, 1, e)) {
                for (GeomRelState<Sub> j : dests(i, 1, e)) {
                    sumRate += pi1m[s]
                            * rate(i, 1, j.subState, 1 + j.rLevel, e);
                }
            }
        }
        return sumRate;
    }

    /**
     * This return the Sub-States class, rather than GeomState.
     * @see jmarkov.SimpleMarkovProcess#getStateClass()
     */
    @Override
    public Class getStateClass() {
        return getInitialState().getSubState().getClass();
    }

    /*
     * SOLVERS MANAGEMENT
     */

    /**
     * Returns the default GeometrixSolver.
     * @return the default GeometrixSolver.
     */
    protected final GeometricSolver getDefaultGeometrixSolver() {
        if (defaultGeometrixSolver == null) {
            defaultGeometrixSolver = new MtjLogRedSolver(this);
        }
        return defaultGeometrixSolver;
    }

    /**
     * The currently defined solver.
     * @see jmarkov.solvers.GeometricSolver
     * @return Returns the GeometrixSolver.
     */
    public GeometricSolver getGeometrixSolver() {
        if (GeometrixSolver == null)
            GeometrixSolver = getDefaultGeometrixSolver();
        return GeometrixSolver;
    }

    /**
     * Allows the user to set an alternate solver.
     * @see jmarkov.solvers.GeometricSolver
     * @param geometrixSolver The GeometrixSolver to set.
     */
    public void setGeometrixSolver(GeometricSolver geometrixSolver) {
        this.GeometrixSolver = geometrixSolver;
        resetResults();
        debug(1, "New geometrix solver set: " + geometrixSolver);
    }

}// End of class
