package jmarkov.jmdp.solvers;

import jmarkov.basic.Action;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDP;

/**
 * This solver solves an average infinite horizon MDP by building and solving a
 * linear problem using as interface Xpress-Optimizer. It requires the
 * professional version of XpressMP.
 * 
 * @author Diego Bello, Germán Riaño - Universidad de Los Andes (C) 2005
 * @param <S>
 *            States class.
 * @param <A>
 *            Actions class.
 * 
 */
public class MPSXpressAverage<S extends State, A extends Action> extends
        MpsLpAverageSolver<S, A> {

    private long processTime = 0;
    private MPSXpressDiscounted<S, A> discSolver = null;

    /**
     * 
     * The constructor method exclusively receives a problem of the type
     * infinite DTMDP and the name that the user wants as path for MPS Files and
     * solutions from Xpress-Optimizer.
     * 
     * @param problem
     *            the structure of the problem of type infinite DTMDP.
     * @param workingDir
     *            The path where will be saved and read the MPS file.
     * @param fileName
     *            The MPS file name.
     */
    public MPSXpressAverage(DTMDP<S, A> problem, String workingDir,
            String fileName) {
        super(problem, workingDir, fileName);
        discSolver = new MPSXpressDiscounted<S, A>(problem, 0.0, workingDir,
                fileName);
        problem.setSolver(this);
    }

    /**
     * This is the default constructor for MPSXpressDiscounted class, and
     * defines the path for the working directory equals to the MPS folder tied
     * to the project. The constructor method exclusively receives a problem of
     * the type infinite DTMDP.
     * 
     * @param problem
     *            the structure of the problem of type infinite DTMDP.
     */
    public MPSXpressAverage(DTMDP<S, A> problem) {
        super(problem);
        discSolver = new MPSXpressDiscounted<S, A>(problem, 0.0);
        problem.setSolver(this);
    }

    @Override
    public String label() {
        return "MPS Xpress Average Solver.";
    }

    public void solveLP() throws SolverException {
         discSolver.solveLP();
    }

    /**
     * @see jmarkov.jmdp.solvers.MPSXpressDiscounted#buildSolution()
     */
    public Solution<S, A> buildSolution() throws SolverException {
        return discSolver.buildSolution();
    }
    
    

}
