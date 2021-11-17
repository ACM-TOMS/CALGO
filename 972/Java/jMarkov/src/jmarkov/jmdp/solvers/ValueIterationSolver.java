package jmarkov.jmdp.solvers;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.DecisionRule;
import jmarkov.basic.Policy;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.CTMDP;
import jmarkov.jmdp.DTMDP;

/*
 * Created on 18/09/2004
 * Last Modified 04/12/2013 by Daniel F. Silva
 */

/**
 * This class belongs to the set of default solvers included in the
 * jmdp package. It extends Solver and should only be used on INFINITE
 * horizon problems. The objective is to be able to return an optimal
 * policy given a problem structure.
 * @author Andres Sarmiento, Germán Riaño and Daniel F. Silva
 * @param <S> States Class.
 * @param <A> Actions class.
 */
public class ValueIterationSolver<S extends State, A extends Action> extends
        AbstractDiscountedSolver<S, A> {

    private double epsilon = 0.0001; // Precision
    private double initVal = 0.0; // initial value function
    private boolean useGaussSeidel = true;
    private boolean useErrorBounds = false;
    /** indicator for average cost objective*/
    private boolean isAverage = false;
    /** stores the differences between value functions in consecutive iterations
    * can be used to track speed of convergence */
    private List<Double> difValues = new ArrayList<Double>();
    /** stores the process time */
    protected long processTime = 0;//
    /** Used to store the number of iterations */
    protected long iterations;
    /**
     * Variable used by bestAction(S). That method returns the value function
     * and stores the action that achieves said value in this field. When
     * solving the average cost problem fixedRelativeCost stores the reference value for each iteration 
     * which is eventually the optimal value for every state
     */
    // TODO: make private
    A bestAction;
    private double fixedRelativeCost = 0.0;
    private ValueFunction<S> relativeValueFunction = new ValueFunction<S>(); 
    private double gain = 0;
    /**
     * Indicates whether we are using modified policy iteration in the average cost case.
     */
    private boolean useModifiedAverage = false;

    boolean printBias = false;
    boolean printGain = false;
    
    /**
     * General constructor for DTMDP. This constructor is not public
     * and is called by other solvers and other constructors.
     * @param problem The problem to solve
     * @param interestRate The interest rate.
     * @param initValue Initial value for all states.
     * @param epsilon Solution tolerance.
     * @param useGaussSeidel whether Gauss Seidel is used.
     * @param useErrorBounds whether Errror bounds are used
     * @param isAverage whether an average problem is solved
     * @param useModifiedAverage whether the Average modification to
     *        prevent periodic problems is used.
     * @see #setInitVal(double)
     * @see #setEpsilon(double)
     * @see #useGaussSeidel
     * @see #useErrorBounds(boolean)
     */
    ValueIterationSolver(DTMDP<S, A> problem, double interestRate,
            double initValue, double epsilon, boolean useGaussSeidel,
            boolean useErrorBounds, boolean isAverage,
            boolean useModifiedAverage) {
        super(problem, interestRate);
        this.initVal = initValue;
        this.epsilon = epsilon;
        this.useGaussSeidel = useGaussSeidel;
        this.useErrorBounds = useErrorBounds;
        this.isAverage = isAverage;
        this.useModifiedAverage = useModifiedAverage;
        problem.setSolver(this);
    }

    /**
     * General constructor for CTMDP. This constructor is not public
     * and is called by other solvers and other constructors.
     * @param problem The problem to solve
     * @param interestRate The interest rate.
     * @param initValue Initial value for all states.
     * @param epsilon Solution tolerance.
     * @param useGaussSeidel whether Gauss Seidel is used.
     * @param useErrorBounds whether Errror bounds are used
     * @param isAverage whether an average problem is solved
     * @param useModifiedAverage whether the Average modification to
     *        prevent periodic problems is used.
     * @see #setInitVal(double)
     * @see #setEpsilon(double)
     * @see #useGaussSeidel
     * @see #useErrorBounds(boolean)
     */
    ValueIterationSolver(CTMDP<S, A> problem, double interestRate,
            double initValue, double epsilon, boolean useGaussSeidel,
            boolean useErrorBounds, boolean isAverage,
            boolean useModifiedAverage) {
        super(problem, interestRate);
        this.initVal = initValue;
        this.epsilon = epsilon;
        this.useGaussSeidel = useGaussSeidel;
        this.useErrorBounds = useErrorBounds;
        this.isAverage = isAverage;
        this.useModifiedAverage = useModifiedAverage;
        problem.setSolver(this);
    }

    /**
     * Default Constructor for Discrte time problems.
     * @param problem the structure of the problem of type DTMDP
     * @param interestRate represents how much less is the reward
     *        received in the next period instead of receiving it in
     *        the present period.
     */
    public ValueIterationSolver(DTMDP<S, A> problem, double interestRate) {
        this(problem, interestRate, 0.0, 0.0001, true, false, false, false);
    }

    /**
     * Default Constructor for continuous time problems.
     * @param problem the structure of the problem of type CTMDP
     * @param interestRate represents how much less is the reward
     *        received in the next period instead of receiving it in
     *        the present period.
     */
    public ValueIterationSolver(CTMDP<S, A> problem, double interestRate) {
        this(problem, interestRate, 0.0, 0.001, true, false, false, false);
    }

    /**
     * Default Constructor for Discrte time average cost problems.
     * This constructor is not public anc can only be called by
     * Relative value iteration solver.
     * @param problem the structure of the problem of type DTMDP
     */
    ValueIterationSolver(DTMDP<S, A> problem, boolean useModifiedAverage) {
        this(problem, 0.0, 0.0, 0.001, true, false, true, useModifiedAverage);
    }

    /**
     * Default Constructor for Discrte time average cost problems.
     * This constructor is not public anc can only be called by
     * Relative value iteration solver.
     * @param problem the structure of the problem of type DTMDP
     */
    ValueIterationSolver(CTMDP<S, A> problem, boolean useModifiedAverage) {
        this(problem, 0.0, 0.0, 0.001, true, false, true, useModifiedAverage);
    }

    /**
     * Value Iteration is a solver method this is theoretically
     * convergent only after infinite iterations. Because of the
     * practical impossibility to do this, the solver is designed to
     * stop when the difference between iterations is as much as
     * epsilon. The smaller epsilon is, the closer the result will be
     * to the actual optimum but it will take a longer time to solve
     * the problem. The default value of epsilon is 0.0001.
     * @param epsilon maximum difference between iterations.
     */
    public synchronized void setEpsilon(double epsilon) {
        this.epsilon = epsilon;
    }

    /**
     * All the states have an initial valueFunction that by default is
     * 1. The solver will converge faster if the initial valueFunction
     * is closer to he optimum. It is a healthy practice to set the
     * initial value in a value that is a bad estimate for the
     * optimum.
     * @param val inital valueFunction for all states.
     */
    private void setInitVal(double val) {
        this.initVal = val;
    }

    /**
     * The GaussSeidel modification of the ValueIteration method is a
     * change that is garanteed to have a performance at least as good
     * as the methos without the modifications. In many problems,
     * specially the ones with many states, the modification can imply
     * a significant improvement. By default it set to true. It
     * provides no significant improvement if used jointly with the
     * ErrorBounds modification.
     * @param val sets whether or not the GaussSeidel modification
     *        will be used.
     * @see #useErrorBounds(boolean)
     */
    public synchronized void useGaussSeidel(boolean val) {
        this.useGaussSeidel = val;
    }

    /**
     * @return Returns the epsilon.
     */
    public final double getEpsilon() {
        return epsilon;
    }

    /**
     * @return Returns the isAverage.
     */
    public final boolean isAverage() {
        return isAverage;
    }

    /**
     * @return Returns true if it uses Modified Average.
     */
    final boolean usesModifiedAverage() {
        return useModifiedAverage;
    }

    /**
     * @return Returns true if uses Error Bounds.
     */
    public final boolean usesErrorBounds() {
        return useErrorBounds;
    }

    /**
     * @return Returns true if Gauss Seidel is active.
     */
    public final boolean usesGaussSeidel() {
        return useGaussSeidel;
    }

    /**
     * The ErrorBounds modification to the ValueIteration method is a
     * change that is garanteed to have a performance at least as good
     * as the methos without the modifications. In many problems,
     * specially the ones with many states, the modification can imply
     * a significant improvement. This method modifies the iteratios
     * and the stopping criterion. It builds upper and lower bounds
     * for the optimal in each iteration and stops when the bounds are
     * only delta apart or less ignoring where the actual
     * valueFunction is. The bounds converge faster than the actual
     * valueFunction. By default it set to false.
     * @param val sets whether or not to use the ErrorBounds
     *        modification.
     */
    public synchronized void useErrorBounds(boolean val) {
        this.useErrorBounds = val;
    }

    /**
     * Solves the problem up to a maximum number of iterations. Main
     * cycle of the algorithm.
     * @param maxIterations maximum number of iterations.
     * @throws SolverException 
     * @see jmarkov.jmdp.solvers.Solver#solve()
     */
    private synchronized Solution<S, A> solve(int maxIterations) throws Exception {
        init();
        double actualDifference = Double.MAX_VALUE;
        long initialTime = 0;
        iterations = 0;
        initialTime = System.currentTimeMillis();
        double bound = (1 - discountFactor) * epsilon / (2 * discountFactor);
        if (isAverage && !useModifiedAverage)
            bound = epsilon;
        while ((actualDifference > bound) && (iterations < maxIterations)) {
            if (useErrorBounds)// si se usan cotas de error
                actualDifference = computeWithErrorBounds();
            else
                // si NO se usan cotas de error
                actualDifference = computeNoErrorBounds();
            difValues.add(new Double(actualDifference));
            getProblem().debug(3,
                    "Max difference from previous value = " + actualDifference);
            iterations++;
        }
        processTime = System.currentTimeMillis() - initialTime;
        if (isAverage){
        	relativeValueFunction = valueFunction;
        	gain = fixedRelativeCost;
        	valueFunction = buildValueFunction(valueFunction, fixedRelativeCost);            
        }
        solved = true;
        return new Solution<S, A>(valueFunction, policy);
    }

    /**
     * Solves the problem.
     * @return returns a Solution with the optimal policy and value
     *         funtion.
     * @throws SolverException
     */
    @Override
    public Solution<S, A> solve() {
        try {
			return solve(Integer.MAX_VALUE);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
    }

    /**
     * Initializes the valueFunction for all the states.
     */
    protected void init() {
        States<S> st = getProblem().getAllStates();
        for (S s : st) {
            valueFunction.set(s, initVal);
        }
    }

    /**
     * Computes an iteration of the Value Iteration Algorithm without
     * the use of error bounds.
     * @return maximum change in value function due to this iteration.
     */
    protected double computeNoErrorBounds() {
        States<S> st = getProblem().getAllStates();
        DecisionRule<S, A> decisionRuleCompute;
        if (iterations > 0)
            decisionRuleCompute = new DecisionRule<S, A>(policy
                    .getDecisionRule());// local decision rule
        else
            decisionRuleCompute = new DecisionRule<S, A>();// empty local decision rule
        ValueFunction<S> vf = new ValueFunction<S>(valueFunction);
        Iterator<Map.Entry<S, Double>> it = vf.iterator();
        Iterator<Map.Entry<S, Double>> itGlobal = valueFunction.iterator();
        Iterator<Map.Entry<S, A>> itDR = decisionRuleCompute.iterator();
        double maxDifference = 0.0;
        int n = 0;
        fixedRelativeCost = 0.0;

        for (S i : st) {
            bestAction = null;
            double newValFunction = bestAction(i);
            if (n == 0 && isAverage)
                fixedRelativeCost = newValFunction;
            double diff;
            if (useGaussSeidel && !isAverage) {
                Map.Entry<S, Double> entryGlobal = itGlobal.next();// over the global Val Function
                diff = Math.abs(newValFunction - entryGlobal.getValue());
                entryGlobal.setValue(newValFunction);// update global value function
            } else {
                Map.Entry<S, Double> entry = it.next();// over local value function
                if (isAverage) {
                    diff = Math.abs(newValFunction - entry.getValue()
                            - fixedRelativeCost);
                    if (useModifiedAverage)
                        // update value function
                        entry.setValue(newValFunction - fixedRelativeCost
                                + (1 - discountFactor) * entry.getValue());
                    else
                        // update value function
                        entry.setValue(newValFunction - fixedRelativeCost);
                } else {
                    diff = Math.abs(newValFunction - entry.getValue());
                    entry.setValue(newValFunction);// update value function
                }
            }
            if (maxDifference < diff)
                maxDifference = diff;
            if (iterations > 0) {
                Map.Entry<S, A> decisionRuleEntry = null;
                decisionRuleEntry = itDR.next();
                // update optimal action
                decisionRuleEntry.setValue(bestAction);
            } else
                // first iteration, no decision rule to update, must be built from scratch
                decisionRuleCompute.set(i, bestAction);
            n++;
        }
        if (!useGaussSeidel || isAverage)
            valueFunction = vf;
        policy = new Policy<S, A>(decisionRuleCompute);
        return maxDifference;
    }

    /**
     * Computes an iteration of the Value Iteration Algorithm with the
     * use of error bounds.
     * @return maximum change in value function due to this iteration.
     */
    protected double computeWithErrorBounds() {
        States<S> st = getProblem().getAllStates();
        DecisionRule<S, A> decisionRuleCompute;
        if (iterations > 0)
            decisionRuleCompute = new DecisionRule<S, A>(policy
                    .getDecisionRule());// update previous decision rule
        else
            decisionRuleCompute = new DecisionRule<S, A>();// first iteration constructs a decision rule from scratch
        ValueFunction<S> vf = new ValueFunction<S>(valueFunction);
        Iterator<Map.Entry<S, Double>> it = vf.iterator();
        Iterator<Map.Entry<S, Double>> itGlobal = valueFunction.iterator();
        Iterator<Map.Entry<S, A>> itDR = decisionRuleCompute.iterator();
        double maxDifference = 0, minDifference = Double.MAX_VALUE;
        int n = 0;
        fixedRelativeCost = 0;

        for (S i : st) {
            bestAction = null;
            double newValueFunction = bestAction(i);
            if (isAverage && n == 0)
                fixedRelativeCost = newValueFunction;
            double diff = 0;
            Map.Entry<S, A> decisionRuleEntry = null;
            if (useGaussSeidel && !isAverage) {
                Map.Entry<S, Double> entryGlobal = itGlobal.next();
                diff = Math.abs(newValueFunction - entryGlobal.getValue());
                entryGlobal.setValue(newValueFunction);
            } else {
                Map.Entry<S, Double> entry = it.next();
                if (isAverage) {
                    diff = Math.abs(newValueFunction - entry.getValue()
                            - fixedRelativeCost);
                    entry.setValue(newValueFunction - fixedRelativeCost);
                } else {
                    diff = Math.abs(newValueFunction - entry.getValue());
                    entry.setValue(newValueFunction);
                }
            }
            if (maxDifference < diff)
                maxDifference = diff;
            if (minDifference > diff)
                minDifference = diff;
            if (iterations > 0) {
                decisionRuleEntry = itDR.next();
                decisionRuleEntry.setValue(bestAction);
            } else
                decisionRuleCompute.set(i, bestAction);
            n++;
        }
        if (!useGaussSeidel || isAverage)
            valueFunction = vf;
        policy = new Policy<S, A>(decisionRuleCompute);
        if (!isAverage)
            return discountFactor / (1 - discountFactor)
                    * (maxDifference - minDifference);
        return (maxDifference - minDifference);
    }

    /**
     * Find the minimal value function for this state and sets the
     * best action to take in state i, in the variable bestAction.
     * @param i state for which the best action is being determined
     * @return the new ValueFunction for this state.
     */
    protected double bestAction(S i) {
        Actions<A> act = getProblem().feasibleActions(i);
        double val = 0.0;
        double minSoFar = Double.MAX_VALUE;
        for (A a : act) {
            val = getProblem().operation(future(i, a, discountFactor),
                    getDiscreteProblem().immediateCost(i, a));
            if (val < minSoFar) {// minimization
                minSoFar = val;
                bestAction = a;
            }
        }
        return minSoFar;
    }

    private ValueFunction<S> buildValueFunction(ValueFunction<S> vf,
            double optimum) {
        ValueFunction<S> result = new ValueFunction<S>(vf);
        Iterator<Map.Entry<S, Double>> it = vf.iterator();
        while (it.hasNext()) {
            Map.Entry<S, Double> m = it.next();
            result.set(m.getKey(), optimum);// + m.getValue());
        }
        return result;
    }

    /**
     * Calculates the difference in value function between the last
     * iteration and the present
     * @param vf value function used
     * @param i current state
     * @param newValFunction new value
     * @return the change in the valueFunction due to this iteration
     *         in the present state i.
     */
    private double difference(ValueFunction<S> vf, S i, double newValFunction) {
        double diff = Math.abs(newValFunction - vf.get(i));
        if (useGaussSeidel)
            valueFunction.set(i, newValFunction);
        else
            vf.set(i, newValFunction);
        return diff;
    }

    // String difsToString() {
    // return difValues.toString();
    // }

    // public double[] getDifs() {
    // int n = difValues.size();
    // double[] vals = new double[n];
    // int i = 0;
    // for (Double d : difValues) {
    // vals[i] = d;
    // i++;
    // }
    // return vals;
    // }

    /**
     * @return Returns the processTime.
     */
    @Override
    public final long getProcessTime() {
        return processTime;
    }

    /**
     * @return Returns the iterations.
     */
    @Override
    public final long getIterations() {
        return iterations;
    }

    @Override
    public String label() {
        StringBuffer buf = new StringBuffer();
        buf.append("Value Iter. Solver ");
        if (isAverage)
            buf.append(" (Avg)");
        else
            buf.append(" (Disc)");
        return buf.toString();
    }

    @Override
    public String description() {
        StringBuffer buf = new StringBuffer();
        if (isAverage) {
            buf.append("Value Iteration Solver for Average Cost problem, ");
            if (useModifiedAverage) {
                buf.append("Factor = " + discountFactor);
            }
        } else {
            buf.append("Value Iteration Solver\n");
            buf.append("Discount Factor = " + discountFactor);
        }
        if (useGaussSeidel)
            buf.append(",\nusing Gauss-Seidel modification\n");
        if (useErrorBounds)
            buf.append(",\nusing Error Bounds convergence\n");
        return buf.toString();
    }
    
    public final double getGain(){
    	return gain;
    }
    
    public final ValueFunction<S> getBias(){
    	return relativeValueFunction;
    }
    
    public void setPrintBias(boolean val) {
        this.printBias = val;
    }
    
    public void setPrintGain(boolean val) {
        this.printGain = val;
    }
    
    /**
     * Prints the solution on a given PrintWriter.
     * 
     * @param pw
     * @see java.io.PrintWriter
     */
    
    @Override
    public void printSolution(PrintWriter pw) {
        pw.println(this);
        try {
            getOptimalPolicy().print(pw);

            if (printValueFunction)
                valueFunction.print(pw);
            if (printBias)
            	pw.println("Bias dor each state:");
            	relativeValueFunction.print(pw);
            if (printGain)
            	pw.println("Gain = " + gain);
            if (printProcessTime) {
                pw.println("Process time = " + getProcessTime()
                        + " milliseconds");
            }
        } catch (SolverException e) {
            pw.print(" Error solving the problem :" +e);
        }
    }

    /**
     * Prints the solution in the default PrintWriter (System.out)
     * 
     * @throws Exception
     */
    @Override
    public void printSolution() throws Exception {
        printSolution(new PrintWriter(System.out));
    }

    
}// class end

