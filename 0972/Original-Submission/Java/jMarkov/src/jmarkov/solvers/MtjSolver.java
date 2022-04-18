/**
 * 
 */
package jmarkov.solvers;

import jmarkov.MarkovProcess;
import jmarkov.basic.exceptions.NotUnichainException;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.BiCG;
import no.uib.cipr.matrix.sparse.BiCGstab;
import no.uib.cipr.matrix.sparse.CGS;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import no.uib.cipr.matrix.sparse.DefaultIterationMonitor;
import no.uib.cipr.matrix.sparse.DiagonalPreconditioner;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;
import no.uib.cipr.matrix.sparse.GMRES;
import no.uib.cipr.matrix.sparse.ILU;
import no.uib.cipr.matrix.sparse.IterationMonitor;
import no.uib.cipr.matrix.sparse.IterativeSolver;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
import no.uib.cipr.matrix.sparse.Preconditioner;
import no.uib.cipr.matrix.sparse.QMR;
import no.uib.cipr.matrix.sparse.SSOR;
import no.uib.cipr.matrix.sparse.SparseVector;

/**
 * This class uses MTJ to solve Steady State probabilities.
 * 
 * @author German Riano. Universidad de los Andes.
 */
public class MtjSolver extends SteadyStateSolver {

    /**
     * 
     * This is the list of preconditioner offered in MTJ.
     * 
     * @author German Riano. Universidad de los Andes.
     * 
     */
    public enum EnumPrecond {
        /** */
        IDEN("IdentityPreconditioner"),
        /** */
        SSOR("SuccessiveOverrelaxation"),
        /** */
        DIA("DiagonalPreconditioner"),
        /** */
        ILU("ILU");

        private String name = "";

        private EnumPrecond(String name) {
            this.name = name;
        }
        /**
         * Returns a label with the solver name.
         * @return The solver name
         * @see jmarkov.solvers.Solver#label()
         */
        public String label() {
            return name;
        }
    }

    /**
     * This is the list of solvers provided by MTJ.
     * 
     * @author German Riano. Universidad de los Andes.
     * 
     */
    public enum EnumSolver {
        /** */
        BiCG("BiConjugateGradients"),
        /** */
        BiCGstab("BiConjugateGradientsStabilized"),
        // Chebyshev( "Chebyshev"),
        // CG("ConjugateGradients"),
        /** */
        CGS("ConjugateGradientsSquared"), /** */
        GMRES("GeneralizedMinimalResiduals"),
        // IR( "IterativeRefinement"),
        /** */
        QMR("QuasiMinimalResiduals");
        private String name = "";

        private EnumSolver() {
        }

        private EnumSolver(String name) {
            this.name = name;
        }

        /**
         * @return The name.
         */
        public String getName() {
            return name;
        }

    };

    private long processTime = -1;
    private EnumPrecond currentPreConditioner = EnumPrecond.IDEN;
    private EnumSolver currentIterSolver = EnumSolver.BiCGstab;
    private IterativeSolver iterativeSolver = null;
    private Preconditioner preConditioner = null;
    private boolean tryOthers = true;
    private Matrix A; // linear system matrix
    private Vector b; // linear system RHS vector
    private Vector pi0; // linear system initial Vector
    private double maxLambda = 0.0;

    /**
     * Default constructor. Uses as default solver BiCGStab first. If it fails,
     * it tries other solvers.
     * 
     * @param mp
     */

    public MtjSolver(MarkovProcess mp) {
        this(mp, EnumSolver.BiCGstab, true);
    }

    /**
     * Construct a solver for the given SimpleMarkovProcess. It uses the given solver
     * enumeration to determine what solver to use.
     * 
     * @param mp
     *            Tha Markov process solved
     * @param solver
     *            The solver to use.
     */
    public MtjSolver(MarkovProcess mp, EnumSolver solver) {
        this(mp, solver, false);
    }

    /**
     * Construct a solver for the given SimpleMarkovProcess. It uses the given solver
     * enumeration to determine what solver to use. If it fails it tries other
     * solvers.
     * 
     * @param mp
     *            the Markov Process to solve.
     * @param solver
     *            the iterative solver form the enumeration IterSolver.
     * @param tryOthers
     *            whether a different solver should be tried if the first one
     *            fails.
     */
    public MtjSolver(MarkovProcess mp, EnumSolver solver, boolean tryOthers) {
        super(mp);
        this.currentIterSolver = solver;
        this.tryOthers = tryOthers;
    }

    /**
     * Generator Matrix
     */
    private Matrix genMatrix = null;

    /**
     * Returns the Generator matrix.
     * 
     * @return the Generator Matrix <B>G</B>.
     */
    public Matrix getGenerator() {
        if (genMatrix != null)
            return genMatrix;
        mp.debug(1, "Building MTJ Generator Matrix ...");
        //TODO: se podra hacer mejor???
        
        //long startTimeGetRates = System.currentTimeMillis();
        Matrix R = new FlexCompRowMatrix(mp.getMtjRates());
        /*long stopTimeGetRates = System.currentTimeMillis();
        long elapsedTimeGetRates = stopTimeGetRates - startTimeGetRates;
        System.out.println("getGenerator - getMtjRates: "+elapsedTimeGetRates+" ms");*/
        
        //long startTimeQ = System.currentTimeMillis();
        int n = R.numRows();
        double sum[] = new double[n];
        Matrix Q = R.copy();
        // First we compute diagonals
        for (MatrixEntry e : Q) {
            int i = e.row();
            int j = e.column();
            if (i != j) {
                sum[i] += R.get(i, j);
            }
        }
        // set the diagonal
        for (int i = 0; i < n; i++) {
            Q.set(i, i, -sum[i]);
            maxLambda = Math.max(sum[i], maxLambda);
        }
        /*long stopTimeQ = System.currentTimeMillis();
        long elapsedTimeQ = stopTimeQ - startTimeQ;
        System.out.println("getGenerator - Q from R: "+elapsedTimeQ+" ms");*/
        
        mp.debug(3, "The Generator matrix is: \n" + Q);
        return (genMatrix = Q);
    }

    /**
     * Returns the A matrix, used to solve pi*A = (0,...,1)
     * 
     * @return The matrix used in the linear system.
     */
    private Matrix getSystemMatrix() {
        if (A == null) {
            Matrix Q = getGenerator();
            A = Q.copy();
            int n = A.numRows();
            // Replace last row by ones..
            for (int i = 0; i < n; i++) {
                A.set(i, n - 1, 1.0);
            }
            mp.debug(4, "The A system matrix is: \n" + Q);
        }
        return A;
    }

    /**
     * Create unit vector b = (0,0,...,0,1)
     * 
     * @return Vector (0,0,...,0,1) used in the linear system.
     */
    private Vector getRightVector() {
        if (b == null) {
            int n = mp.getNumStates();
            b = new SparseVector(n, 1);// one non-zero
            b.set(n - 1, 1.0); // set one at the last position.
            mp.debug(4, "The b vector is: \n" + b);
        }
        return b;
    }

    /**
     * Returns an initial probability guess vector.
     * 
     * @return pi0 = (1/n) (1,1,...,1)
     */
    private Vector getInitialVector() {
        int n = mp.getNumStates();
        if (pi0 == null) {
            pi0 = new DenseVector(n); // to hold solution
            double pr = (1.0 / n);
            for (int i = 0; i < n; i++)
                pi0.set(i, pr);
        }
        return pi0;
    }

    @Override
    public double[] getSteadyState() throws NotUnichainException {

        mp.debug(1, "Computing  Steady State Probabilities with solver " + this
                + " ...");
        A = getSystemMatrix();
        b = getRightVector();
        pi0 = getInitialVector();

        mp.debug(4, "Solving: pi A = b ..." + b);
        DenseVector pi = null;
        try {
            pi = solve(A, b, pi0, currentIterSolver);
        } catch (IterativeSolverNotConvergedException e) {
            if (tryOthers) {
                mp.debug(1, "Solver " + currentIterSolver + " failed ...");
                for (EnumSolver solv : EnumSolver.values()) {
                    if (solv != currentIterSolver) {
                        mp.debug(1, "Trying solver " + solv + "...");
                        try {
                            pi = solve(A, b, pi0, solv);
                        } catch (IterativeSolverNotConvergedException e2) {
                            pi = null;
                        }
                        if (pi != null)
                            return pi.getData();
                    }
                }
            }
            if (pi == null) {
                throw new NotUnichainException(
                        "Exception solving Steady State probabilitis with solver "
                                + this + "\n(Reason :" + e.getReason()
                                + ", Iterations: " + e.getIterations() + ")");
            }
        }
        return pi.getData();

    }

    /**
     * Solves pi * A = b for pi, using the given solver.
     * 
     * @param A
     * @param b
     * @param pi0
     * @param iterSolver
     * @return pi
     */
    private DenseVector solve(Matrix A, Vector b, Vector pi0,
            EnumSolver iterSolver) throws IterativeSolverNotConvergedException {
        iterativeSolver = getIterativeSolver(pi0, iterSolver);
        // IterationMonitor im = iteartiveSolver.getIterationMonitor();

        preConditioner = getPreConditioner(A);
        DenseVector sol = null;

        if (currentPreConditioner != EnumPrecond.IDEN)
            iterativeSolver.setPreconditioner(preConditioner);
        long initialTime = System.currentTimeMillis();
        // solve it!
        //TODO: que pasó con solveTranspose?
        sol = new DenseVector(iterativeSolver.solve(A.transpose(), b, pi0), false);
        processTime = System.currentTimeMillis() - initialTime;
        mp.debug(1, "Solve time : " + processTime + " milliseconds.");

        return sol;
    }

    /**
     * @see jmarkov.solvers.Solver#label()
     */
    @Override
    public String label() {
        return "MTJ Solver: PreCond =  " + currentPreConditioner
                + " Iter Solver = " + currentIterSolver;
    }

    /**
     * @return Returns the currentIterSolver.
     */
    public EnumSolver getCurrentIterSolver() {
        return currentIterSolver;
    }

    /**
     * Sets the solver to use.
     * 
     * @param iterSolver
     *            The currentIterSolver to set.
     * @param tryOthers
     *            whether other solvers should be tryed if this fails.
     */
    public void setIterSolver(EnumSolver iterSolver, boolean tryOthers) {
        this.currentIterSolver = iterSolver;
        this.tryOthers = tryOthers;
        mp.resetResults();
    }

    /**
     * Sets the solver to use. It will not try other solvers if this one fails.
     * 
     * @param iterSolver
     */
    public void setCurrentIterSolver(EnumSolver iterSolver) {
        setIterSolver(iterSolver, false);
    }

    /**
     * @return Returns the currentPreConditioner.
     */
    public EnumPrecond getCurrentPreConditioner() {
        return currentPreConditioner;
    }

    /**
     * @param preConditioner
     *            The currentPreConditioner to set.
     */
    public void setCurrentPreConditioner(EnumPrecond preConditioner) {
        this.currentPreConditioner = preConditioner;
        mp.resetResults();
    }

    /**
     * @return Returns true if the solver shall try other solvers when it fails.
     */
    public final boolean isTryOthers() {
        return tryOthers;
    }

    /**
     * Sets whether the solver shall try other solvers when it fails.
     * 
     * @param tryOthers
     *            true if the solver shall try other solvers when it fails.
     */
    public final void setTryOthers(boolean tryOthers) {
        this.tryOthers = tryOthers;
    }

    /**
     * @return Returns the generator Matrix.
     */
    public final Matrix getGenMatrix() {
        return genMatrix;
    }

    /**
     * @return Returns the Process Time of the last solved problem.
     */
    public final long getProcessTime() {
        return processTime;
    }

    /**
     * @param pi The probability vector.
     * @return Returns the iteartiveSolver.
     */
    public IterativeSolver getIteartiveSolver(Vector pi) {
        return getIterativeSolver(pi, currentIterSolver);
    }

    /**
     * @param pi0 Initial guess value.
     * @param solver The solver used.
     * @return Returns the iterativeSolver.
     */
    public IterativeSolver getIterativeSolver(Vector pi0, EnumSolver solver) {
        if (this.iterativeSolver == null) {
            switch (solver) {
            case BiCG:
                iterativeSolver = new BiCG(pi0);
                break;
            case BiCGstab:
                iterativeSolver = new BiCGstab(pi0);
                break;
            // case (CG):
            // solver = new CG(pi);
            // break;
            case CGS:
                iterativeSolver = new CGS(pi0);
                break;
            // case (Chebyshev):
            // solver = new Chebyshev(pi, 0.0, maxLambda);
            // break;
            case GMRES:
                iterativeSolver = new GMRES(pi0);
                break;
            // case (IR):
            // solver = new IR(pi);
            // break;
            case QMR:
                iterativeSolver = new QMR(pi0);
                break;
            }
            IterationMonitor mon = new DefaultIterationMonitor(100000, 1e-6,
                    1e-50, 1e+5);
            iterativeSolver.setIterationMonitor(mon);
        }
        return iterativeSolver;
    }

    /**
     * @param A The system matrix.
     * @return Returns the preConditioner.
     */
    private Preconditioner getPreConditioner(Matrix A) {
        if (preConditioner == null) {

            switch (currentPreConditioner) {
            case IDEN:
                preConditioner = null;
                break;
            case SSOR:
                preConditioner = new SSOR((CompRowMatrix)A);
                break;
            case DIA:
                preConditioner = new DiagonalPreconditioner(A.numRows());
                break;
            case ILU:
                preConditioner = new ILU((CompRowMatrix)A);
                break;
            }
        }
        return preConditioner;
    }

    public String description() {
        return "This Solvers use the methods implemented in MTJ";
    }

}
