package jmarkov.jmdp.solvers;

import java.io.File;

import jmarkov.basic.Action;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDP;

/**
 * This class builds the Dual Linear Program for an average infinite
 * horizon MDP in a MPS file. It is an abstract class and needs to be implemented.
 * @author Diego Bello, Germán Riaño - Universidad de Los Andes (C)
 *         2005
 * @param <S> States class.
 * @param <A> Action class.
 */

public abstract class MpsLpAverageSolver<S extends State, A extends Action>
        extends AbstractAverageSolver<S, A> implements MpsLpSolver<S, A> {

    private static String fs = System.getProperty("file.separator");
    private MpsLpDiscountedSolver<S, A> discSolver = null;

    /**
     * The constructor method receives a problem of the type infinite
     * DTMDP, the working directory where the MPS file will be stored,
     * and the name that the user wants for the MPS File.
     * @param problem The problem to be solved.
     * @param workingDir Where the file will be created.
     * @param fileName Label for the MPS File.
     */

    public MpsLpAverageSolver(DTMDP<S, A> problem, String workingDir,
            String fileName) {
        super(problem);
        discSolver = new MpsLpDiscountedSolver<S, A>(problem, 0.0, workingDir,
                fileName, true) {

            @Override
            public void solveLP() throws SolverException {
                solveLPI();
            }

            @Override
            public long getIterations() {
                return getIterationsI();
            }

            @Override
            public String label() {
                return "MPS Xpress Solver (Avg)";
            }

            @Override
            public Solution<S, A> buildSolution() throws SolverException {
                return buildSolutionI();
            }
        };
    }

    /**
     * This cosntructor creates a solver for this problem. The created
     * mps file is stored in a temp folder.
     * @param problem The structure of the problem of type infinite
     *        DTMDP.
     */

    public MpsLpAverageSolver(DTMDP<S, A> problem) {
        this(problem, System.getProperty("user.dir"), "MDP.mps");
    }

    public final String getMpsFileName() {
        return discSolver.getMpsFileName();
    }

    @Override
    public final long getIterations() {
        return discSolver.getIterations();
    }

    @Override
    public final Solution<S, A> solve() throws SolverException {
        return discSolver.solve();
    }

    /**
     * @see jmarkov.jmdp.solvers.MpsLpSolver#getWorkingDir()
     */
    public final File getWorkingDir() {
        return discSolver.getWorkingDir();
    }

    public final File getMpsFile() {
        return discSolver.getMpsFile();
    }

    /* Internal Bridge methods */

    /**
     * @return Returns the discSolver.
     */
    protected MpsLpDiscountedSolver<S, A> getDiscSolver() {
        return discSolver;
    }

    private void solveLPI() throws SolverException {
        solveLP();
    }

    private Solution<S, A> buildSolutionI() throws SolverException {
        return buildSolution();
    }

    private long getIterationsI() {
        return getIterations();
    }

    public long getBuildTime() {
        return discSolver.getBuildTime();
    }

    public long getLpSolveTime() {
        return discSolver.getLpSolveTime();
    }

    public long getSolBuildTime() {
        return discSolver.getSolBuildTime();
    }

    @Override
    public final long getProcessTime() {
        return discSolver.getProcessTime();
    }

}
