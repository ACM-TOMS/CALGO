package jmarkov.jmdp.solvers;

import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Map;

import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.DecisionRule;
import jmarkov.basic.Policy;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.NonStochasticException;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDP;
import jmarkov.jmdp.solvers.AbstractAverageSolver;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.sparse.BiCG;
import no.uib.cipr.matrix.sparse.BiCGstab;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;
import no.uib.cipr.matrix.sparse.IterativeSolver;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
import no.uib.cipr.matrix.sparse.SparseVector;

/*
 * Created on 04/12/2013
 */

/**
 * This class solves infinite horizon non-discounted problems using the
 * policy iteration algorithm. It extends Solver and should only be
 * used on infinite horizon problems. The objective function the
 * solver uses is the long-run-Average cost. The result is a deterministic
 * optimal policy for the given structure. Policy Iteration is a
 * solver method this is always convergent in a finite number of
 * iterations, under the right conditions. The algorithm has to solve a linear system of equations
 * as big as the amount of states. When there are too many states, it
 * is recommendable to use other solvers, or using the modified policy
 * iteration (by using the second constructor).  The advantage of using
 * Policy Iteration is that the result is the true optimal solution
 * and not an aproximation as in other common methods. The method
 * starts with a policy. It solves the system of linear equations for
 * the value functions for that policy. With this values it looks for
 * a better policy. It then solves the value functions again and looks
 * for a better policy. If this policy is equal to the last policy
 * tried, it stops, in any other case it keeps improving the policy
 * and updating the value functions.
 * @author Daniel F. Silva
 * @param <S> States class.
 * @param <A> Actions class.
 */

public class PolicyIterationSolverAvg<S extends State, A extends Action> extends
        AbstractAverageSolver<S, A> {

    // DenseVector vecValueFunction = null;

    /** Used to store local states */
    // protected List<S> localStates = new ArrayList<S>();
    private DenseVector costs;
    private boolean isOptimal = false;

    /** Used to store the number of iterations */
    protected long iterations;
    /** Used to store process time */
    protected long processTime = 0;

    private DenseVector vecValueFunction = null;
    private FlexCompRowMatrix matrix = null;
    private DecisionRule<S, A> currentDecisionRule = null;

    private ValueFunction<S> relativeValueFunction = new ValueFunction<S>(); 
    private double gain = 0;
    boolean printBias = false;
    boolean printGain = false;
    

    /**
     * The constructor method exclusively receives a problem of the
     * type InfiniteMDP because this solver is only designed to work
     * on infinite horizon problems. This solver solves the average
     * objective function problem. There is no discount factor.
     * @param problem the structure of the problem of type InfiniteMDP
     */
    public PolicyIterationSolverAvg(DTMDP<S, A> problem) {
        super(problem);
    }


    @Override
    public Solution<S, A> solve() throws SolverException {
        long initialTime = System.currentTimeMillis();
        currentDecisionRule = initialDecisionRuleFirst();
        policy = new Policy<S, A>(currentDecisionRule);

        vecValueFunction = new DenseVector(getDiscreteProblem().getNumStates());
        matrix = buildMatrix(currentDecisionRule);

        iterations = 0;
        while (!isOptimal) {
            problem.debug(2, "Iteration " + iterations);
            getProblem().debug(3, "Current Rule = " + currentDecisionRule);
            valueFunction = policyEvaluation();
            getProblem().debug(3, "Current Value function = " + valueFunction);
            currentDecisionRule = policyImprovement();
            policy.setDecisionRule(currentDecisionRule);
            iterations++;
            // assert (iterations < 20);
        }
        policy = new Policy<S, A>(currentDecisionRule);
        updateResults(valueFunction);
        solved = true;
        processTime = System.currentTimeMillis() - initialTime;
        return new Solution<S, A>(valueFunction, policy);
    }


	/**
     * Builds the initial policy, setting the first available action
     * to each state.
     */
    private DecisionRule<S, A> initialDecisionRuleFirst() {
        valueFunction = new ValueFunction<S>();
        DecisionRule<S, A> localDecisionRule = new DecisionRule<S, A>();
        States<S> states = getProblem().getAllStates();
        for (S i : states) {
            Actions<A> availableActions = getProblem().feasibleActions(i);
            // build an initial policy
            Iterator<A> it = availableActions.iterator();
            if(it.hasNext()){
            	A a = it.next();
            	localDecisionRule.set(i, a);
            	if (i.getIndex()==9){
                	a = it.next();
                	localDecisionRule.set(i, a);
            	}
            }
            valueFunction.set(i, 0.0);
        }
        return localDecisionRule;
    }


    private ValueFunction<S> policyEvaluation() throws SolverException {
            valueFunction = solveMatrix();
        return valueFunction;
    }

    private DecisionRule<S, A> policyImprovement() throws SolverException {

        States<S> sts = getProblem().getAllStates();
        DecisionRule<S, A> newDecisionRule;
        newDecisionRule = new DecisionRule<S, A>(currentDecisionRule);
        double val;
        Iterator<Map.Entry<S, A>> itCurDR = currentDecisionRule.iterator();
        Iterator<Map.Entry<S, A>> itNewDR = newDecisionRule.iterator();
        for (S i : sts) {
            Actions<A> actions = getProblem().feasibleActions(i);
            A bestAction = null;
            double bestValue = Double.MAX_VALUE;
            for (A a : actions) {
                val = getProblem().operation(
                        getDiscreteProblem().immediateCost(i, a),
                        future(i, a));
                if (val < bestValue) {// minimization
                    bestValue = val;
                    bestAction = a;
                }
            }

            	Map.Entry<S, A> curDRentry = itCurDR.next();
                A curAction = curDRentry.getValue();
                if (!bestAction.equals(curAction)) {
                    matrix.setRow(i.getIndex(), buildRowVector(i, bestAction));
                    costs.set(i.getIndex(), getDiscreteProblem().immediateCost(
                            i, bestAction));
                }

            Map.Entry<S, A> newDRentry = itNewDR.next();
            newDRentry.setValue(bestAction);
        }
        for (S i : sts) {
            matrix.set(i.getIndex(), 0, 1);
        }
        isOptimal = (currentDecisionRule.equals(newDecisionRule));
        return newDecisionRule;
    }

    /**
     * Builds the i-th row vector of the matrix (I-beta P)
     * corresponding to the given action.
     * @param i the state
     * @param a the action
     * @return A SparseVector.
     */
    private SparseVector buildRowVector(S i, A a) {
        int n = getDiscreteProblem().getNumStates();
        States<S> reachableStates = getDiscreteProblem().reachable(i, a);
        SparseVector vec = new SparseVector(n, reachableStates.size());
        double sum = 0.0;
        for (S j : reachableStates) {
            double probability = getDiscreteProblem().prob(i, j, a);
            sum += probability;
            assert (probability >= 0);
            // make sure is the same.
            j = getDiscreteProblem().getAllStates().get(j);
            if (probability > 0) {
                vec.set(j.getIndex(), probability);
            }

        }
        vec.scale(-1);
        vec.add(i.getIndex(), 1.0);
        if (Math.abs(sum - 1.0) > 1e-5) {
            throw new NonStochasticException(
                    "Probabilities do not add up to 1 for state " + i
                            + ", and action " + a + ", sum = " + sum);
        }
        return vec;
    }

    /**
     * This method builds the Probability Transision matrix for a
     * specified policy. The solver then transforms this matrix and
     * uses it to solve the value functions for each state.
     * @param currentDecisionRule the policy under which the
     *        probability matrix is to be built.
     * @return the probability matrix.
     */

    private FlexCompRowMatrix buildMatrix(DecisionRule<S, A> currentDecisionRule) {

    	StatesSet<S> stts = getDiscreteProblem().getAllStates();
        int n = stts.size();
         
        FlexCompRowMatrix matrix = new FlexCompRowMatrix(n, n);
        costs = new DenseVector(n);
        for (S i : stts) {
            A a = currentDecisionRule.getAction(i);
            matrix.setRow(i.getIndex(), buildRowVector(i, a));
            costs.set(i.getIndex(), getDiscreteProblem().immediateCost(i, a));
        }
        for (S i : stts) {
            matrix.set(i.getIndex(), 0, 1);
        }
        return matrix;
    }

    /**
     * This method is used by the PolicyIterationSolver to solve the
     * linear system of equations to determine the value functions of
     * each state for a given policy.
     * @return a DenseVector (type defined in the JMP package
     *         documentation) with the value functions for each state.
     *         The index for each state are the same ones determined
     *         in the localStates ArrayList
     * @throws SolverException
     */
    protected ValueFunction<S> solveMatrix() throws SolverException {
        getProblem().debug(4, "Matrix to solve:\n" + matrix);
        try {
            IterativeSolver solver = new BiCG(vecValueFunction);
            solver.solve(matrix, costs, vecValueFunction);
        } catch (IterativeSolverNotConvergedException e) {
            throw new SolverException(
                    "Policy iteration Solver: error solving linear system.", e);
        }
        return buildValueFunction(vecValueFunction);
    }

    /**
     * Expected value of valueFunction for the current state and a
     * specified action.
         */
    protected final double future(S i, A a, 
            ValueFunction<S> vf) {
        double sum = 0.0;
        States<S> reachableStates = getDiscreteProblem().reachable(i, a);
        for (S j : reachableStates) {
            sum += getDiscreteProblem().prob(i, j, a) * getGain();
        }
        return  sum;
    }
    
    /**
     * Expected value of valueFunction for the current state and a
     * specified action.
     * @param i The State
     * @param a Action taken
     * @return Expected value of valueFunction.
     */
    protected final double future(S i, A a) {
        return future(i, a, getValueFunction());
    }


    private ValueFunction<S> buildValueFunction(DenseVector vec) {
        ValueFunction<S> vf = new ValueFunction<S>();
        States<S> stts = getDiscreteProblem().getAllStates();
        int i = 0;
        for (S s : stts) {
            vf.set(s, vec.get(i));
            i++;
        }
        return vf;
    }


    protected final double future(S i, A a, double discountF,
            ValueFunction<S> vf) {
        double sum = 0.0;
        States<S> reachableStates = getDiscreteProblem().reachable(i, a);
        for (S j : reachableStates) {
            sum += getDiscreteProblem().prob(i, j, a) * vf.get(j);
        }
        return discountF * sum;
    }
    @Override
    public String description() {
        return "Policy Iteration Solver\n";
    }

    @Override
    public String label() {
        return "Policy Iter. Solver(avg)";
    }

    //
    // Overriden methods
    //

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

    private void updateResults(ValueFunction<S> valueFunction) {
		
    	StatesSet<S> stts = getDiscreteProblem().getAllStates();
    	Iterator<S> it = stts.iterator();
    	S s = it.next();
        gain = valueFunction.get(s);
        relativeValueFunction.set(s,0);
		while (it.hasNext()){
			s = it.next();
			relativeValueFunction.set(s,valueFunction.get(s));
		}
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
            	pw.println("Bias for each state:");
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
