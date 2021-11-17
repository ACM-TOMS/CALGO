package jmarkov.jmdp.solvers;

import java.io.File;
import java.io.IOException;
import java.util.Formatter;

import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDP;

/**
 * This class builds a Linear Program for a discounted infinite
 * horizon MDP in a MPS file. A extending class must code solveLP
 * method in order to solve the problem.
 * This is an abstract class, so one must implement it to give the MPS file to a solver.
 * @see #solveLP()
 * @author Diego Bello, Germán Riaño - Universidad de Los Andes (C)
 *         2005
 * @param <S> States class.
 * @param <A> Actions class.
 */

public abstract class MpsLpDiscountedSolver<S extends State, A extends Action>
        extends AbstractDiscountedSolver<S, A> implements MpsLpSolver<S, A> {
//TODO: One could write a simple implementation of this class that allows the user to only build the MPS file without linking to a Solver. 
    private Formatter output;
    private File mpsFile = null;
    private static String fs = File.separator;
    // private Map<S, String> statesNamesMap = new TreeMap<S,
    // String>();
    private boolean isAvg = false;
    private long processTime = -1;
    private long buildTime = -1;
    private long lpSolveTime = -1;
    private long solBuildTime = -1;

    /**
     * The constructor method exclusively receives a problem of the
     * type infinite DTMDP , an interest rate that is modified for
     * being used as discount factor and the name that the user wants
     * for the MPS File.
     * @param problem The structure of the problem of type infinite
     *        DTMDP.
     * @param interestRate A rate which is paid for the use of a
     *        resource.
     * @param workingDir Where the MPS file will be created.
     * @param fileName Name for the MPS File (with no path).
     */
    public MpsLpDiscountedSolver(DTMDP<S, A> problem, double interestRate,
            String workingDir, String fileName) {
        super(problem, interestRate);
        mpsFile = new File(workingDir, fileName);
    }

    /**
     * This is the default constructor for MpsLpDiscountedSolver
     * class, and defines the label MDP for the MPS File. The
     * constructor method exclusively receives a problem of the type
     * infinite DTMDP , an interest rate that is modified for being
     * used as discount factor.
     * @param problem The structure of the problem of type infinite
     *        DTMDP.
     * @param interestRate A rate which is paid for the use of a
     *        resource.
     */

    public MpsLpDiscountedSolver(DTMDP<S, A> problem, double interestRate) {
        super(problem, interestRate);
        try {
            File folder = new File(System.getProperty("user.dir"));
            mpsFile = null;
            if (folder.exists() && folder.isDirectory())
                mpsFile = File.createTempFile("MDP", ".mps", folder);
            else
                mpsFile = File.createTempFile("MDP", ".mps");
        } catch (IOException e) {
            // Exception will be caught at build time. (hopefully)
        }
        mpsFile.deleteOnExit();
    }

    /**
     * The constructor is used by the partenr average solver.
     * @param problem The structure of the problem of type infinite
     *        DTMDP.
     * @param interestRate A rate which is paid for the use of a
     *        resource.
     * @param workingDir Where the MPS file will be created.
     * @param fileName Name for the MPS File (with no path).
     * @param isAverage True if an average model is being built.
     */
    protected MpsLpDiscountedSolver(DTMDP<S, A> problem, double interestRate,
            String workingDir, String fileName, boolean isAverage) {
        super(problem, interestRate);
        mpsFile = new File(workingDir, fileName);
        isAvg = isAverage;
    }

    /**
     * The constructor is used by the partner average solver.
     * @param problem The structure of the problem of type infinite
     *        DTMDP.
     * @param interestRate A rate which is paid for the use of a
     *        resource.
     * @param isAverage True if an average model is being built.
     */

    protected MpsLpDiscountedSolver(DTMDP<S, A> problem, double interestRate,
            boolean isAverage) {
        super(problem, interestRate);
        isAvg = isAverage;
        try {
            File folder = new File(System.getProperty("user.dir"));
            mpsFile = null;
            if (folder.exists() && folder.isDirectory())
                mpsFile = File.createTempFile("MDP", ".mps", folder);
            else
                mpsFile = File.createTempFile("MDP", ".mps");
            mpsFile.deleteOnExit();
        } catch (IOException e) {
            // Exception will be caught at build time. (hopefully)
        }
    }

    /**
     * @return Returns the MPS File Name.
     */
    public final String getMpsFileName() {
        return mpsFile.getPath();
    }

    /**
     * @return Returns the mpsFile.
     */
    public File getMpsFile() {
        return mpsFile;
    }

    /**
     * Returns the working directory (where the MPS file is located)
     * @return Returns the MPS File folder.
     */
    public final File getWorkingDir() {
        return mpsFile.getParentFile();
    }

    /**
     * @return Returns the true if an Average problem is being solved.
     */
    public final boolean isAvg() {
        return isAvg;
    }

    /**
     * The implementator classes should override this class to solve
     * the problem using the mpsFile that has been created.
     */
    public abstract void solveLP() throws SolverException;

    /**
     * The implementator classes should override this class to build
     * the solution after the model has been solved.
     * @return The solution to the problem.
     */
    public abstract Solution<S, A> buildSolution() throws SolverException;

    @Override
    public final Solution<S, A> solve() throws SolverException {
        long startTime = System.currentTimeMillis();
        build();
        long endBuild = System.currentTimeMillis();
        buildTime = endBuild - startTime;
        solveLP();
        long endLP = System.currentTimeMillis();
        lpSolveTime = endLP - endBuild;
        Solution<S, A> sol = buildSolution();
        valueFunction = sol.getValueFunction();
        policy = sol.getPolicy();
        solBuildTime = System.currentTimeMillis() - endLP;
        solved = true;
        processTime = System.currentTimeMillis() - startTime;
        return sol;
    }

    /**
     * Builds the MPS file.
     * @throws SolverException
     */
    private void build() throws SolverException {
        DTMDP<S, A> problem = getDiscreteProblem();
        int numStates = problem.getAllStates().size();
        double alpha = (discountFactor == 1.0) ? 0.0 : 1.0 / numStates;
        // int indexStates = 0;
        // States<S> states = problem.getAllStates();
        // for (S state : states) {
        // statesNamesMap.put(state, "S" + indexStates);
        // indexStates++;
        // }
        openFile();
        writeName();
        writeRows();
        writeColumns();
        writeRHS(alpha);
        writeEnData();
        closeFile();
    }

    public final long getBuildTime() {
        return buildTime;
    }

    public final long getLpSolveTime() {
        return lpSolveTime;
    }

    public final long getSolBuildTime() {
        return solBuildTime;
    }

    @Override
    public final long getProcessTime() {
        return processTime;
    }

    /**
     * Open a new file if is posible.
     * @throws SolverException
     */
    private void openFile() throws SolverException {
        String fileName = mpsFile.getPath();
        try {
            File dir = mpsFile.getParentFile();
            if (!dir.exists()) {
                dir.mkdirs();
            }
            fileName = mpsFile.getPath();
            output = new Formatter(mpsFile);
        } // end try
        catch (SecurityException securityException) {
            throw new SolverException("You do not have write access to  file: "
                    + fileName);
        } // end catch
        catch (IOException fileNotFoundException) {
            throw new SolverException("Error creating file: " + fileName);
        } // end catch
    } // end method openFile

    private void writeName() {
        output.format("NAME          " + mpsFile.getName() + "\n");
    }

    /**
     * This method writes the ROWS section, and define the objective
     * function with the label OBJ, the normalization equation
     * (average case) with the label NORM, the other constraints are
     * labeled with statesMPS mask.
     */
    private void writeRows() {
        output.format("ROWS\n");
        output.format(" N  OBJ\n");
        States<S> states = getDiscreteProblem().getAllStates();
        for (S s : states) {
            String constraintLabel = "S" + s.getIndex();
            // don't include first row on AVG problems.
            boolean notFirstRow = !(constraintLabel.equals("S0"));
            if (!isAvg || notFirstRow)
                output.format(" E  %-8.8s\n", constraintLabel);
        }

        if (isAvg)
            output.format(" E  NORM\n");

    }

    private void writeColumns() {
        DTMDP<S, A> problem = getDiscreteProblem();
        output.format("COLUMNS\n");
        int indexVariables = 0;
        StatesSet<S> states = problem.getAllStates();
        for (S stateI : states) {
            Actions<A> actions = problem.feasibleActions(stateI);
            for (A action : actions) {
                double cost = problem.immediateCost(stateI, action);
                String variableName = "V" + String.valueOf(indexVariables);
                output.format("    %-8.8s  %-8.8s  %-10.8g\n", variableName,
                        "OBJ", cost);
                States<S> statesJ = problem.reachable(stateI, action);
                for (S stateJ : statesJ) {
                    double probability = problem.prob(stateI, stateJ, action);
                    double coefficient = ((stateI.equals(stateJ)) ? 1.0 : 0.0)
                            - ((isAvg) ? probability : discountFactor
                                    * probability);
                    stateJ = states.get(stateJ);
                    String constraintLabel = "S" + stateJ.getIndex();
                    // don't include first row on AVG problems.
                    boolean notFirstRow = !(constraintLabel.equals("S0"));
                    if (!isAvg || notFirstRow)
                        output.format("    %-8.8s  %-8.8s  %-10.8g\n",
                                variableName, constraintLabel, coefficient);

                }
                if (isAvg) {
                    output.format("    %-8.8s  %-8.8s  %-10.8g\n",
                            variableName, "NORM", 1.0);
                }
                indexVariables++;
            }
        }

    }

    private void writeRHS(double alpha) {
        output.format("RHS\n");
        for (S state : getDiscreteProblem().getAllStates()) {
            String constraintLabel = "S" + state.getIndex();
            // don't include first row on AVG problems.
            boolean notFirstRow = !(constraintLabel.equals("S0"));
            if (!isAvg || notFirstRow)
                output.format("    %-8.8s  %-8.8s  %-10.5g\n", "RHS",
                        constraintLabel, alpha);
        }
        if (isAvg)
            output.format("    %-8.8s  %-8.8s  %-10.5g\n", "RHS", "NORM", 1.0);
    }

    /**
     * Writes ENDATA in MPS File.
     */
    private void writeEnData() {
        output.format("ENDATA");
    }

    /**
     * Is closed the file.
     */
    private void closeFile() {
        if (output != null)
            output.close();
    }

}
