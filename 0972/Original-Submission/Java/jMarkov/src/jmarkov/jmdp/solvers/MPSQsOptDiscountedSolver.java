/**
 * MPSQsOptDiscountedSolver.java
 * Created: Dec 24, 2005
 */
package jmarkov.jmdp.solvers;

import java.io.IOException;

import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.DecisionRule;
import jmarkov.basic.Policy;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.CT2DTConverter;
import jmarkov.jmdp.CTMDP;
import jmarkov.jmdp.DTMDP;
import qs.Problem;
import qs.QS;
import qs.QSException;

/**
 * This solver solves an discounted infinite horizon MDP by building and solving a
 * linear problem using as interface QSopt-Optimizer. It implements MpsLpDiscountedSolver
 * creates an MPS file and feeds it to the solver. 
 * 
 * @author Germán Riaño. Universidad de los Andes. (C) 2005
 * @param <S>
 *            States class.
 * @param <A>
 *            Actions class.
 * 
 */
public class MPSQsOptDiscountedSolver<S extends State, A extends Action>
        extends MpsLpDiscountedSolver<S, A> {
    private Problem lpProb = null;

    /**
     * @param problem
     *            The problem to solve
     * @param interestRate
     *            The interest rate
     * @param workingDir
     *            The working directory where MPS file will be stored.
     * @param fileName
     *            The name for the generated MPD file.
     */
    public MPSQsOptDiscountedSolver(CTMDP<S, A> problem, double interestRate,
            String workingDir, String fileName) {
        this(problem, interestRate, workingDir, fileName, false);
    }

    /**
     * @param problem
     *            The problem to solve
     * @param interestRate
     *            The interest rate
     * @param workingDir
     *            The working directory where MPS file will be stored.
     * @param fileName
     *            The name for the generated MPD file.
     */
    public MPSQsOptDiscountedSolver(DTMDP<S, A> problem, double interestRate,
            String workingDir, String fileName) {
        super(problem, interestRate, workingDir, fileName);
    }

    /**
     * @param problem
     *            Probelem to solve
     * @param interestRate
     *            The interest rate.
     */
    public MPSQsOptDiscountedSolver(DTMDP<S, A> problem, double interestRate) {
        this(problem, interestRate, System.getProperty("user.dir"), "MDP.mps",
                false);
    }

    /**
     * 
     * @param problem
     *            problem to solve
     * @param interestRate
     *            interest rate
     */
    public MPSQsOptDiscountedSolver(CTMDP<S, A> problem, double interestRate) {
        this(problem, interestRate, System.getProperty("user.dir"), "MDP.mps",
                false);
    }

    /**
     * @param problem
     *            The problem to solve
     * @param interestRate
     *            The interest rate
     * @param workingDir
     *            The working directory where MPS file will be stored.
     * @param fileName
     *            The name for the generated MPD file.
     * @param isAvg
     *            true if problem is treated as Average cost problem
     */
    public MPSQsOptDiscountedSolver(DTMDP<S, A> problem, double interestRate,
            String workingDir, String fileName, boolean isAvg) {
        super(problem, interestRate, workingDir, fileName, isAvg);
    }

    /**
     * @param problem
     *            The problem to solve
     * @param interestRate
     *            The interest rate
     * @param workingDir
     *            The working directory where MPS file will be stored.
     * @param fileName
     *            The name for the generated MPD file.
     * @param isAvg
     *            true if problem is treated as Average cost problem
     */
    public MPSQsOptDiscountedSolver(CTMDP<S, A> problem, double interestRate,
            String workingDir, String fileName, boolean isAvg) {
        super(new CT2DTConverter<S,A>(problem), interestRate, workingDir, fileName,
                isAvg);
    }

    /**
     * @param problem
     *            The problem to solve
     * @param interestRate
     *            The interest rate
     * @param isAvg
     *            true if problem is treated as Average cost problem
     */
    public MPSQsOptDiscountedSolver(DTMDP<S, A> problem, double interestRate,
            boolean isAvg) {
        super(problem, interestRate, isAvg);
    }
    /**
     * @param problem
     *            The problem to solve
     * @param interestRate
     *            The interest rate
     * @param isAvg
     *            true if problem is treated as Average cost problem
     */
    public MPSQsOptDiscountedSolver(CTMDP<S, A> problem, double interestRate,
            boolean isAvg) {
        super(new CT2DTConverter<S,A>(problem), interestRate, isAvg);
    }

    /** Loads the LP problem */
    void loadProblem() throws SolverException {
        String mpsFile = getMpsFileName();
        try {
            boolean isMpsFormat = true;
            lpProb = Problem.read(mpsFile, isMpsFormat);
            if (lpProb == null) {
                throw new SolverException("Could not parse problem with file "
                        + mpsFile);
            }
        } catch (IOException e) {
            throw new SolverException("Problems reading file " + mpsFile
                    + ". Exception: " + e);
        }
    }

    /**
     * @see jmarkov.jmdp.solvers.AbstractInfiniteSolver#getIterations()
     */
    @Override
    public long getIterations() {
        return 0;
    }

    @Override
    public void solveLP() throws SolverException {
        loadProblem();
        try {
            lpProb.opt_primal();
        } catch (QSException e) {
            throw new SolverException("Error solving with QS:" + e);
        }
       //return getSolution();
    }

    @Override
    public Solution<S, A> buildSolution() throws SolverException {
        DecisionRule<S, A> dr = new DecisionRule<S, A>();
        ValueFunction<S> vf = new ValueFunction<S>();
        try {
            int status = lpProb.get_status();
            switch (status) {
            case (QS.LP_OPTIMAL):
                int nrows = lpProb.get_rowcount();
                double pi[] = new double[nrows];
                double slack[] = new double[nrows];
                int ncols = lpProb.get_colcount();
                double x[] = new double[ncols];
                double rc[] = new double[ncols];
                // now read the solution
                double val = lpProb.get_solution(x, pi, slack, rc);
                States<S> states = getProblem().getAllStates();
                int i = 0,
                var = 0;
                for (S s : states) {
                    // set value function and increase counter
                    vf.set(s, (isAvg()) ? val : pi[i++]);
                    Actions<A> acts = getProblem().feasibleActions(s);
                    for (A a : acts) {
                        if (x[var++] > 0) {
                            dr.set(s, a);// this is the optimal action
                        }
                    }
                }
                break;
            case (QS.LP_PRIMAL_INFEASIBLE):
                throw new SolverException(
                        "QSopt found the problem to be unfeasible.");
            case (QS.LP_PRIMAL_UNBOUNDED):
                throw new SolverException(
                        "QSopt found the problem to be unbounded.");
            default:
                throw new SolverException(
                        "QSopt did not find an optimal solution. Return value: "
                                + lpProb.get_status());
            }

        } catch (QSException e) {
            throw new SolverException("Error reading the solution from QS.", e);
        }
        return new Solution<S, A>(vf, new Policy<S, A>(dr));
    }

    /**
     * @see jmarkov.jmdp.solvers.Solver#label()
     */
    @Override
    public String label() {
        return "QSOpt Solver (disc)";
    }

}
