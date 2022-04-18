package jmarkov.jmdp.solvers;

import jmarkov.basic.Action;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDP;

/**
 * This solver solves a average-cost infinite horizon MDP by building
 * and solving a linear problem using as interface Xpress BCL. It
 * requires the professional version of XpressMP and the JAVA build
 * path must include the xprb.jar library, for running the
 * applications.
 * @author Diego Bello, Germán Riaño - Universidad de Los Andes (C)
 *         2005
 * @param <S> States Class
 * @param <A> Action class
 */

public class LPBCLAverageSolver<S extends State, A extends Action> extends
        AbstractAverageSolver<S, A> {

    /** Used to store the process time */
    private long processTime = -1;
    /** Used to store the build time */
    private long buildTime = -1;
    /** Used to store the Linear Programming solve time */
    protected long lpSolveTime = -1;
    /** used to store the iterations */
    private long solBuildTime = -1;

    /**
     * The constructor method exclusively receives a problem of the
     * type DTMDP because this solver is only designed to work on
     * infinite discrete horizon problems. This solver solves an
     * average DTMDP.
     * @param problem the structure of the problem of type DTMDP
     */
    public LPBCLAverageSolver(DTMDP<S, A> problem) {
        super(problem);
    }

    @Override
    public long getIterations() {
        return 0;
    }

    /**
     * Linear Programming Average Solver is a tool that builds the
     * solution based on the MDP's mathematical background given by
     * Puterman and the software provided by XpressMP (BCL libraries).
     * Is mandatory for the use to have a Xpress professional version.
     * It actually just builds a LPBCLDiscountedSolver object without a 
     * discount factor, all the heavy lifting is done by that class.
     * @throws SolverException
     */
    @Override
    public Solution<S, A> solve() throws SolverException {
        LPBCLDiscountedSolver<S, A> discSolver = new LPBCLDiscountedSolver<S, A>(
                getDiscreteProblem());
        discSolver.printValueFunction = this.printValueFunction;
        discSolver.printProcessTime = this.printProcessTime;
        Solution<S, A> answer = discSolver.solve();
        solved = true;
        buildTime = discSolver.getBuildTime();
        lpSolveTime = discSolver.getLpSolveTime();
        solBuildTime = discSolver.getBuildTime();
        processTime = discSolver.getProcessTime();
        return answer;
    }

    @Override
    public String label() {
        return "BCL Solver (avg)";
    }

    /**
     * @return Returns the processTime.
     */
    @Override
    public long getProcessTime() {
        return processTime;
    }

    /**
     * @return Returns the build time.
     */
    public long getBuildTime() {
        return buildTime;
    }

    /**
     * @return Returns the Linear Programming Solve Time
     */
    public long getLpSolveTime() {
        return lpSolveTime;
    }

    
    /**
     * Returns the time needed to build the Solution after the LP was
     * solved.
     * @return Returns the solBuildTime.
     */
    public long getSolBuildTime() {
        return solBuildTime;
    }
    
}
