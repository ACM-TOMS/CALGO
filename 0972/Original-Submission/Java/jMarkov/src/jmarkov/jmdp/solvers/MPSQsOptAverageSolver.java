/**
 * MPSQsOptAverageSolver.java
 * Created: Mar 13, 2006
 */
package jmarkov.jmdp.solvers;

import jmarkov.basic.Action;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.CT2DTConverter;
import jmarkov.jmdp.CTMDP;
import jmarkov.jmdp.DTMDP;

/**
 * This solver solves an average infinite horizon MDP by building and solving a
 * linear problem using as interface QSopt-Optimizer. It implements MpsLpDiscountedSolver
 * creates an MPS file and feeds it to the solver.
 *  
 * @author German Riano. Universidad de los Andes. (C) 2006
 * @param <S> Stete class
 * @param <A> Action class
 */
public class MPSQsOptAverageSolver<S extends State, A extends Action> extends
        MpsLpAverageSolver<S, A> {

    private MPSQsOptDiscountedSolver<S, A> discSolver = null;

    /**
     * @param problem
     */
    public MPSQsOptAverageSolver(DTMDP<S, A> problem) {
        super(problem);
        discSolver = new MPSQsOptDiscountedSolver<S, A>(problem, 0.0, true);
    }

    /**
     * Creates a solver for average cost problems.
     * @param problem The problem to solve.
     * @param workingDir working directory where files will be saved
     * @param fileName MDP file name.
     */
    public MPSQsOptAverageSolver(DTMDP<S, A> problem, String workingDir,
            String fileName) {
        super(problem, workingDir, fileName);
        discSolver = new MPSQsOptDiscountedSolver<S, A>(problem, 0.0,
                workingDir, fileName, true);
    }

    /**
     * @param problem
     */
    public MPSQsOptAverageSolver(CTMDP<S, A> problem) {
        super(new CT2DTConverter<S, A>(problem));
        discSolver = new MPSQsOptDiscountedSolver<S, A>(problem, 0.0, true);
    }

    /**
     * Creates a solver for average cost problems.
     * @param problem The problem to solve.
     * @param workingDir working directory where files will be saved
     * @param fileName MDP file name.
     */
    public MPSQsOptAverageSolver(CTMDP<S, A> problem, String workingDir,
            String fileName) {
        super(new CT2DTConverter<S, A>(problem), workingDir, fileName);
        discSolver = new MPSQsOptDiscountedSolver<S, A>(problem, 0.0,
                workingDir, fileName, true);
    }

    /**
     * @see jmarkov.jmdp.solvers.Solver#label()
     */
    @Override
    public String label() {
        return "QsOpt Solver (AVG Cost)";
    }

    public void solveLP() throws SolverException {
        discSolver.solveLP();
    }

    /**
     * @see jmarkov.jmdp.solvers.MPSQsOptDiscountedSolver#buildSolution()
     */
    public Solution<S, A> buildSolution() throws SolverException {
        return discSolver.buildSolution();
    }

}
