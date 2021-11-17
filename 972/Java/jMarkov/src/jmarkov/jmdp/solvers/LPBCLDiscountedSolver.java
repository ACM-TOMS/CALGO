package jmarkov.jmdp.solvers;

import java.util.Iterator;

import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.DecisionRule;
import jmarkov.basic.Policy;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDP;

import com.dashoptimization.XPRB;
import com.dashoptimization.XPRBctr;
import com.dashoptimization.XPRBlinExp;
import com.dashoptimization.XPRBprob;
import com.dashoptimization.XPRBvar;

/**
 * This solver solves a discounted infinite horizon MDP by building
 * and solving a linear problem using as interface Xpress BCL. It
 * requires the professional version of XpressMP and the JAVA build
 * path must include the xprb.jar libray, for running the
 * applications.
 * @author Diego Bello, Germán Riaño. Universidad de los Andes.
 * @param <S> States class.
 * @param <A> Actions class.
 */
public class LPBCLDiscountedSolver<S extends State, A extends Action> extends
        AbstractDiscountedSolver<S, A> implements LPSolver {
    /** Used to store the process time */
    private long processTime = -1;
    /** Used to store the build time */
    private long buildTime = -1;
    /** Used to store the Linear Programming solve time */
    private long lpSolveTime = -1;
    /** used to store the iterations */
    private long solBuildTime = -1;
    private boolean isAvg = false;

    private boolean showXpressOutput = false;
    private XPRBprob linearProblem = null;// The Xpress LP
    private int solveResult;// Xpress result

    /**
     * The constructor method receives a problem of the type infinite
     * DTMDP and an interest rate that is modified for being used as
     * discount factor. The discount factor and the problem gives
     * necesary information for solving a discounted MDP.
     * @param problem the structure of the problem of type infinite
     *        DTMDP.
     * @param interestRate A rate which is paid for the use of a
     *        resource.
     */

    public LPBCLDiscountedSolver(DTMDP<S, A> problem, double interestRate) {
        super(problem, interestRate);
    }

    /**
     * The constructor method receives a problem of the type infinite
     * DTMDP and solves it using the average cost criteria. The
     * constructor is not public and it is used only by the partner
     * solver.
     * @param problem the structure of the problem of type infinite
     *        DTMDP.
     * @param interestRate A rate which is paid for the use of a
     *        resource.
     */
    LPBCLDiscountedSolver(DTMDP<S, A> problem) {
        this(problem, 0.0);
        isAvg = true;
    }

    @Override
    public long getIterations() {
        return 0;
    }

    @Override
    public Solution<S, A> solve() throws SolverException {
        long startTime = System.currentTimeMillis();
        try {
            linearProblem = new XPRBprob();
            if (!showXpressOutput)
                linearProblem.setMsgLevel(0);
            createObjectiveFunctionAndConstraints(linearProblem);
            long endBuild = System.currentTimeMillis();
            buildTime = endBuild - startTime;
            solveLP();
            long endLP = System.currentTimeMillis();
            lpSolveTime = endLP - endBuild;
            if (solveResult == 0) {
                buildSolution();
                solBuildTime = System.currentTimeMillis() - endLP;
            } else
                throw new SolverException("Xpress failed to solve the problem");
            processTime = System.currentTimeMillis() - startTime;
            return new Solution<S, A>(valueFunction, policy);
        } catch (UnsatisfiedLinkError e) {
            throw new SolverException(
                    "This solver requires XPess professional edition.", e);
        } catch (NoClassDefFoundError e) {
            throw new SolverException(
                    "This solver requires XPess professional edition.", e);
        }
    }

    private void createObjectiveFunctionAndConstraints(XPRBprob linearProblem) {
        // TODO: This needs to be commented more specifically!!
        XPRBlinExp objectiveFunction = new XPRBlinExp();
        XPRBctr constraintNorm = null;// normalization constraint.

        StatesSet<S> statesI = getDiscreteProblem().getAllStates();
        // Map<S, Integer> statesMap = new TreeMap<S, Integer>();
        double alpha = ((isAvg) ? 0.0 : 1.0)
                / (getDiscreteProblem().getAllStates().size());
        // int indexStates = 0;
        int indexVariables = 0;
        int n = getDiscreteProblem().getNumStates();
        for (int indexStates = 0; indexStates < n; indexStates++) {
            // statesMap.put(stateI, indexStates);

            XPRBctr constraint = linearProblem.newCtr(String
                    .valueOf(indexStates));
            constraint.add(alpha);
            constraint.setType(XPRB.E);
            // indexStates++;
        }
        if (isAvg) {
            constraintNorm = linearProblem.newCtr("NORM");
            // constraintNorm = linearProblem.getCtrByName("NORM");
            constraintNorm.setType(XPRB.E);
            constraintNorm.add(1.0);
        }
        for (S stateI : statesI) {
            Actions<A> actions = getProblem().feasibleActions(stateI);
            for (A action : actions) {
                double cost = getDiscreteProblem()
                        .immediateCost(stateI, action);
                String variableName = String.valueOf(indexVariables);
                XPRBvar variable = linearProblem.newVar(variableName, XPRB.PL,
                        0, Double.MAX_VALUE);
                objectiveFunction.add(variable.mul(cost));
                States<S> statesJ = getDiscreteProblem().reachable(stateI,
                        action);
                for (S stateJ : statesJ) {
                    double probability = getDiscreteProblem().prob(stateI,
                            stateJ, action);
                    double coefficient = ((stateI.equals(stateJ)) ? 1.0 : 0.0)
                            - (discountFactor * probability);
                    // String constraintLabel =
                    // String.valueOf(statesMap
                    // .get(stateJ));
                    stateJ = statesI.get(stateJ); 
                    String constraintLabel = String.valueOf(stateJ.getIndex());
                    XPRBctr constraint = linearProblem
                            .getCtrByName(constraintLabel);
                    constraint.add(variable.mul(coefficient));
                }
                if (isAvg) {
                    constraintNorm.add(variable.mul(1.0));
                }
                indexVariables++;
            }
        }
        if (isAvg) {
            XPRBctr delConstraint = linearProblem.getCtrByName(String
                    .valueOf(0));
            linearProblem.delCtr(delConstraint);
        }

        linearProblem.setObj(objectiveFunction);
    }

    private ValueFunction<S> buildValueFunction(XPRBprob linearProblem) {
        int indexStates = 0;
        States<S> states = getDiscreteProblem().getAllStates();
        if (discountFactor < 1.0) {
            for (S state : states) {
                double dualValue = 0;
                dualValue = linearProblem.getCtrByName(
                        String.valueOf(indexStates)).getDual();
                valueFunction.set(state, dualValue);
                indexStates++;
            }
        } else {
            double gain = linearProblem.getObjVal();
            for (S state : states) {
                valueFunction.set(state, gain);
            }
        }
        return valueFunction;

    }

    private DecisionRule<S, A> computeDecisionRule(XPRBprob linearProblem) {
        DecisionRule<S, A> decisionRule = new DecisionRule<S, A>();
        States<S> states = getDiscreteProblem().getAllStates();
        int notTransient = 0;
        int indexVariables = 0;
        for (S state : states) {
            Actions<A> actions = getProblem().feasibleActions(state);
            notTransient = 0;
            for (A action : actions) {
                double answer = linearProblem.getVarByName(
                        String.valueOf(indexVariables)).getSol();
                if (answer > 0.0) {
                    notTransient++;
                    decisionRule.set(state, action);
                }
                indexVariables++;

            }
            if (notTransient == 0) {
                Iterator<A> actionsIterator = actions.iterator();
                decisionRule.set(state, actionsIterator.next());
            }
        }
        return decisionRule;

    }

    @Override
    public String label() {
        return "BCL Solver" + ((isAvg) ? " (avg)" : " (disc)");
    }

    @Override
    public String description() {
        StringBuffer buf = new StringBuffer();
        if (isAvg) {
            buf
                    .append("Linear Programming Solver (BCL) for Average Cost Problem");
        }

        else {
            buf
                    .append("Linear Programming Solver (BCL) for Discounted Cost Problem");
            buf.append(", discount Factor = " + discountFactor);
        }
        return buf.toString();
    }

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

    /**
     * @see jmarkov.jmdp.solvers.LPSolver#solveLP()
     */
    public void solveLP() throws SolverException {
        solveResult = linearProblem.minim("");
    }

    /**
     * @see jmarkov.jmdp.solvers.LPSolver#buildSolution()
     */
    public Solution<S, A> buildSolution() throws SolverException {
        solved = true;
        DecisionRule<S, A> decisionRule = null;
        policy = new Policy<S, A>(decisionRule);
        decisionRule = computeDecisionRule(linearProblem);
        policy.setDecisionRule(decisionRule);
        valueFunction = buildValueFunction(linearProblem);
        return new Solution<S, A>(valueFunction, policy);
    }

}
