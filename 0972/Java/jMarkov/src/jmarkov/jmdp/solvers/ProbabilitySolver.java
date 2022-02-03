/*
 * Created on 10/08/2005
 */
package jmarkov.jmdp.solvers;

import java.util.Iterator;
import java.util.Map;

import jmarkov.basic.Action;
import jmarkov.basic.DecisionRule;
import jmarkov.basic.State;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.CT2DTConverter;
import jmarkov.jmdp.CTMDP;
import jmarkov.jmdp.DTMDP;

/**
 * This class is designed to calculate the long run probabilities of infinite
 * horizon problem. It uses Jacobian and Power methods for sparse matrixes. Of
 * course this class is not needed To do this jmdp should pass the problem
 * to JMarkov
 * 
 * @author Andres Sarmiento and Germán Riaño - Universidad de Los Andes
 * 
 * @param <S>
 *            State classs
 * @param <A>
 *            Actions class.
 */
// TODO: Replace this class and use jMarkov directly for this purpose
// @Deprecated
public class ProbabilitySolver<S extends State, A extends Action> {

    DecisionRule<S, A> dr;
    private ValueFunction<S> probability = null;
    DTMDP<S, A> problem;
    double epsilon = 0.000001;
    int size = 0;
    boolean gaussSeidel = false;
    boolean jacobi = true;
    long processTime = 0;
    int iterations;
    boolean solved = false;
    int maxIterations = 100;

    /**
     * Initializes a new solver for discrete chains
     * 
     * @param problem
     *            discrete time , infinite horizon problem
     * @param dr
     *            decision rule to be evaluated
     */
    public ProbabilitySolver(DTMDP<S, A> problem, DecisionRule<S, A> dr) {
        this.dr = dr;
        this.problem = problem;
    }

    /**
     * Initializes a new solver for discrete chains and solves the probabilities
     * for the optimal policy.
     * 
     * @param problem
     *            discrete time , infinite horizon problem
     * @throws SolverException 
     */
    public ProbabilitySolver(DTMDP<S, A> problem) throws SolverException {
        this.problem = problem;
        this.dr = problem.getOptimalPolicy().getDecisionRule();
    }

    /**
     * Initializes a new solver for continuous chains and solves the
     * probabilities for a particular decision rule.
     * 
     * @param problem
     *            continuous time , infinite horizon problem
     * @param dr
     */
    public ProbabilitySolver(CTMDP<S, A> problem, DecisionRule<S, A> dr) {
        this.dr = dr;
        this.problem = new CT2DTConverter<S, A>(problem);
    }

    /**
     * Initializes a new solver for continuous chains and solves the
     * probabilities for the optimal policy.
     * 
     * @param problem
     *            continuous time , infinite horizon problem
     * @throws SolverException 
     */
    public ProbabilitySolver(CTMDP<S, A> problem) throws SolverException {
        this.problem = new CT2DTConverter<S, A>(problem);
        this.dr = problem.getOptimalPolicy().getDecisionRule();
    }

    /**
     * 
     * @return true if the probabilities were calculated.
     */
    public boolean isSolved() {
        return solved;
    }

    /**
     * Solves the probabilities
     * 
     */
    public void solve() {
        if (jacobi)
            solveJacobi();
        else
            solvePower();
        solved = true;
    }

    /**
     * Solves discrete chains
     * 
     * @return probabilities
     */
    private ValueFunction<S> solveJacobi() {
        ValueFunction<S> oldVals = initializeProbs(); // equiprobable
        ValueFunction<S> newVals = null;
        double maxDifference = Double.MAX_VALUE;
        iterations = 0;
        long initialTime = System.currentTimeMillis();
        while ((maxDifference > epsilon) && (iterations < maxIterations)) { // WISHLIST:identificar
            // periodicidad
            newVals = new ValueFunction<S>(oldVals);
            ValueFunction<S> piTimesP = piTimesP(oldVals);
            Iterator<Map.Entry<S, Double>> oldIt = oldVals.iterator();
            Iterator<Map.Entry<S, Double>> newIt = newVals.iterator();
            Iterator<Map.Entry<S, Double>> piXpIt = piTimesP.iterator();
            Iterator<Map.Entry<S, A>> drIt = dr.iterator();

            maxDifference = 0;
            while (oldIt.hasNext()) {
                Map.Entry<S, Double> oldE = oldIt.next();
                Map.Entry<S, Double> piXpE = piXpIt.next();
                Map.Entry<S, A> drE = drIt.next();
                Map.Entry<S, Double> newE = newIt.next();

                S i = oldE.getKey();
                double pii = problem.prob(i, i, drE.getValue());
                if (pii == 1)
                    System.out.println("State " + i
                            + " is an absorbing state under action"
                            + drE.getValue());
                //WISHLIST:throw NotUnichainException
                double oldVal = oldE.getValue();
                double newVal = (-piXpE.getValue() + oldVal * pii) / (pii - 1);
                if (Math.abs(-piXpE.getValue() + oldVal * pii) < 1E-10)
                    newVal = 0;
                // double absDiff = Math.abs(newVal-oldVal);
                double diff = Math.abs(newVal - oldVal) / oldVal;
                if (maxDifference < diff)
                    maxDifference = diff;
                if (gaussSeidel)
                    oldE.setValue(newVal);
                else
                    newE.setValue(newVal);
            }
            if (!gaussSeidel)
                oldVals = newVals;
            iterations++;
        }
        processTime = System.currentTimeMillis() - initialTime;
        probability = newVals;
        problem.debug(1, "Probability convergence in " + iterations
                + " iterations");
        problem.debug(1, "Probability convergence in " + processTime
                + " milliseconds");
        return newVals;
    }

    // ValueFunction<S> solveMTJ(){
    // PolicyIterationSolver<S,A> Psolver = new
    // PolicyIterationSolver<S,A>(problem, 0.9);
    // SparseRowColumnMatrix M = Psolver.buildMatrix(dr);
    // SequentialBLAS theSeq = new SequentialBLAS();
    // theSeq.addDiagonal(-1, M);
    // // theSeq.add
    // DenseVector vecValueFunction = new DenseVector(size);
    // LinearSolver solver = new BiCG();
    // // solver.solve(M, costs, vecValueFunction);
    //
    // // return Psolver.buildValueFunction(vecValueFunction);
    // }

    /**
     * Solves discrete chains
     */
    private ValueFunction<S> solvePower() {
        ValueFunction<S> readPi = initializeProbs(); // equiprobable
        ValueFunction<S> writePi = null;
        double maxDifference = Double.MAX_VALUE;
        int iterations = 0;
        long initialTime = System.currentTimeMillis();
        while (maxDifference > epsilon) {
            writePi = piTimesP(readPi);
            // check convergence
            maxDifference = difference(readPi, writePi);
            readPi = writePi;
            iterations++;
        }
        probability = writePi;
        processTime = System.currentTimeMillis() - initialTime;
        problem.debug(1, "Probability convergence in " + iterations
                + " iterations");
        problem.debug(1, "Probability convergence in " + processTime
                + " milliseconds");
        return writePi;
    }

    /**
     * Used for discrete chains operations
     * 
     * @param readPi
     * @return probabilities
     */
    private ValueFunction<S> piTimesP(ValueFunction<S> readPi) {
        ValueFunction<S> writePi = initializeVF(); // ceros
        Iterator<Map.Entry<S, Double>> readIt = readPi.iterator();
        while (readIt.hasNext()) {
            Iterator<Map.Entry<S, Double>> itW = writePi.iterator();
            Map.Entry<S, Double> readE = readIt.next();
            S i = readE.getKey();
            for (S j : problem.reachable(i, dr.getAction(i))) {
                Map.Entry<S, Double> entryW = itW.next();
                while (!entryW.getKey().equals(j))
                    entryW = itW.next();
                entryW.setValue(entryW.getValue() + readE.getValue()
                        * (problem.prob(i, j, dr.getAction(i))));
            }
        }
        return writePi;
    }

    /**
     * Evaluates convergence
     * 
     * @param oldPi
     * @param newPi
     * @return difference
     */
    private double difference(ValueFunction<S> oldPi, ValueFunction<S> newPi) {
        Iterator<Map.Entry<S, Double>> oldIt = oldPi.iterator();
        Iterator<Map.Entry<S, Double>> newIt = newPi.iterator();
        double maxDifference = 0;
        while (oldIt.hasNext()) {
            double oldVal = oldIt.next().getValue();
            double diff = Math.abs(newIt.next().getValue() - oldVal) / oldVal;
            if (maxDifference < diff)
                maxDifference = diff;
        }
        return maxDifference;
    }

    /**
     * The GaussSeidel modification of the ValueIteration method is a change
     * that is garanteed to have a performance at least as good as the methods
     * without the modifications. In many problems, specially the ones with many
     * states, the modification can imply a significant improvement. By default
     * it set to true. It provides no significant improvement if used jointly
     * with the ErrorBounds modification.
     * 
     * @param val
     *            sets whether or not the GaussSeidel modification will be used.
     */
    public void setGaussSeidel(boolean val) {
        this.gaussSeidel = val;
    }

    /**
     * 
     * @param val
     *            true to use jacobi methods
     */
    public void setJacobi(boolean val) {
        this.jacobi = val;
    }

    /**
     * initialization in ceros
     * 
     * @return probabilities in ceros
     */
    ValueFunction<S> initializeVF() {
        Iterator<Map.Entry<S, A>> it = dr.iterator();
        ValueFunction<S> probability = new ValueFunction<S>("Steady state probabilities");
        size = 0;
        while (it.hasNext()) {
            Map.Entry<S, A> entry = it.next();
            probability.set(entry.getKey(), 0.0);
            size++;
        }
        return probability;
    }

    /**
     * equiprobable initialization
     * 
     * @return equiprobable probabilities
     */
    ValueFunction<S> initializeProbs() {
        ValueFunction<S> probability = initializeVF();
        Iterator<Map.Entry<S, Double>> it = probability.iterator();
        while (it.hasNext()) {
            Map.Entry<S, Double> entry = it.next();
            entry.setValue((double) 1 / size);
        }
        return probability;
    }

    /**
     * @return Returns the probability.
     */
    public ValueFunction<S> getProbability() {
        return probability;
    }

}
