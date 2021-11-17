package examples.jmdp;

import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Map.Entry;

import jmarkov.MarkovProcess;
import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDP;

/*
 * Created on 28/09/2004
 *
 */

/**
 * A company receives orders at a random rate. At each stage, the company can
 * decide whether to process no orders incurring in a cost of unfilled order per
 * order per stage, or to process all the pending orders incurring in a setup
 * cost.
 * 
 * This class belongs to the examples supplied in the package jmdp. The
 * objective of this file is to show as clear as possible a simple way to use
 * the jmdp package as a tool for solving problems. The complete details of the
 * present problems are explained in the documentation.
 * 
 * See Bersekas, "Dynamic Programming and Control".
 * 
 * @author Andres Sarmiento and Germán Riaño - Universidad de Los Andes
 */
public class OrderProcessing extends DTMDP<PendingOrders, ProcessAction> {

    double setupCost = 0;
    double unfilledOrderCost = 0;
    int maxOrders = 0;
    double theta = 0;
    double[] demandProbability;
    double[] demandCDF;

    ProcessAction all = new ProcessAction("Process all orders"),
            none = new ProcessAction("Process no orders"),
            vec[] = { all, none };
    Actions<ProcessAction> actions = new ActionsSet<ProcessAction>(vec);

    /**
     * @param initialPendingOrders
     *            Number of initial pending orders
     * @param setupCost
     *            Setup cost
     * @param unfilledOrderCost
     *            Cost of unfilled order per order per stage
     * @param MaxOrders
     *            Maximum orders pending
     * @param theta
     *            Demand rate
     */
    public OrderProcessing(int initialPendingOrders, double setupCost,
            double unfilledOrderCost, int MaxOrders, double theta) {
        super(initialize(initialPendingOrders));
        this.setupCost = setupCost;
        this.unfilledOrderCost = unfilledOrderCost;
        this.maxOrders = MaxOrders;
        this.theta = theta;
        initializeProbabilities();
        // initStates();
    }

    private void initializeProbabilities() {
        demandProbability = new double[maxOrders + 1];// P{demand t = i}
        // P{demand t <= i}
        demandCDF = new double[maxOrders + 1];
        double p = Math.exp(-theta);
        demandProbability[0] = p;
        demandCDF[0] = 1; // P{demand t >= 0}
        for (int i = 1; i <= maxOrders; i++) {
            demandCDF[i] = demandCDF[i - 1] - p;
            p = p * theta / i; // P{demand t = i}
            demandProbability[i] = p;
        }
    }

    void initStates() {
        StatesSet<PendingOrders> states = new StatesSet<PendingOrders>();
        for (int n = 0; n <= maxOrders; n++)
            states.add(new PendingOrders(n));
        super.states = states;
    }

    @Override
    public States<PendingOrders> reachable(PendingOrders i, ProcessAction a) {
        StatesSet<PendingOrders> st = new StatesSet<PendingOrders>();
        if (a == all) {
            for (int n = 0; n < i.getOrders(); n++)
                st.add(new PendingOrders(n));
        }
        for (int n = i.getOrders(); n <= maxOrders; n++)
            st.add(new PendingOrders(n));
        return st;
    }

    @Override
    public double immediateCost(PendingOrders i, ProcessAction a) {
        if (a == all)
            return setupCost;
        return i.getOrders() * unfilledOrderCost; // if (a == none)
    }

    @Override
    public Actions<ProcessAction> feasibleActions(PendingOrders i) {
        if (i.getOrders() == 0)
            return new ActionsSet<ProcessAction>(none);
        else if (i.getOrders() == maxOrders)
            return new ActionsSet<ProcessAction>(all);
        else
            return actions;
    }

    @Override
    public double prob(PendingOrders i, PendingOrders j, ProcessAction a) {
        if (a == none) {
            if (j.getOrders() == maxOrders)
                return demandCDF[maxOrders - i.getOrders()];
            return demandProbability[j.getOrders() - i.getOrders()];
        } else if (a == all) {
            if (j.getOrders() == maxOrders)
                return demandCDF[maxOrders];
            return demandProbability[j.getOrders()];
        }
        return 0;
    }

    /**
     * This method defines the index of MOPS
     * 
     * @author German Riano. Universidad de los Andes. (C) 2006
     * 
     */
    protected enum Measure {
        /** Average order level */
        AVERAGE_ORDER_LEVEL;
    }

    /**
     * @param i
     *            Number of pending orders
     * @param m
     *            Required MOP
     * @return Average order level
     */
    public double computeMOPs(PendingOrders i, Measure m) {
        switch (m) {
        case AVERAGE_ORDER_LEVEL:
            return i.getOrders();
        }
        return 0.0;
    }

    /**
     * @return Average of pending orders
     * @throws SolverException
     */
    public double averageOrders() throws SolverException {
        double queue = 0;
        ValueFunction<PendingOrders> probs = getSteadyStateProbabilities();
        Iterator<Entry<PendingOrders, Double>> it = probs.iterator();
        Entry<PendingOrders, Double> e;
        for (; it.hasNext();) {
            e = it.next();
            queue += e.getValue()
                    * computeMOPs(e.getKey(), Measure.AVERAGE_ORDER_LEVEL);
        }
        return queue;
    }

    /**
     * @param initPendingOrders
     * @return A States set representation of number of pending orders
     */
    public static StatesSet<PendingOrders> initialize(int initPendingOrders) {
        PendingOrders initState = new PendingOrders(initPendingOrders);
        return new StatesSet<PendingOrders>(initState);
    }

    /**
     * @param initPendingOrders
     * @return A state representation of this pnumber of pendinf orders.
     */
    public static PendingOrders initialState(int initPendingOrders) {
        return new PendingOrders(initPendingOrders);
    }

    /**
     * @param args
     * @throws SolverException
     */
    public static void main(String[] args) throws SolverException {
        double setupCost = 60;
        double unfilledOrderCost = 4;
        int MaxOrders = 15;
        double theta = 3;
        // double discountFactor = 0.9;
        OrderProcessing problem = new OrderProcessing(0, setupCost,
                unfilledOrderCost, MaxOrders, theta);

        // PolicyIterationSolver<PendingOrders, ProcessAction> theSolver;
        // theSolver = new PolicyIterationSolver<PendingOrders, ProcessAction>(
        // problem, discountFactor);
        // // ValueIterationSolver<PendingOrders,ProcessAction> theSolver3 = new
        // //
        // ValueIterationSolver<PendingOrders,ProcessAction>(problem,discountFactor);
        // theSolver.setPrintValueFunction(true);
        // problem.setSolver(theSolver);
        // problem.solve();
        // problem.printSolution();

        // AverageRewardSolver<PendingOrders,ProcessAction> theSolver2 = new
        // AverageRewardSolver<PendingOrders,ProcessAction>(problem);
        // theSolver2.setPrintValueFunction(false);
        // theSolver2.setGaussSeidel(false);
        // theSolver2.setPrintProcessTime(true);
        // problem.solve();
        // problem.printSolution();

        // theSolver2 = new
        // AverageRewardSolver<PendingOrders,ProcessAction>(problem);
        // theSolver2.setPrintValueFunction(false);
        // theSolver2.useErrorBounds(true);
        // theSolver2.setPrintProcessTime(true);
        // problem.solve();
        // problem.printSolution();

        // theSolver2 = new
        // AverageRewardSolver<PendingOrders,ProcessAction>(problem, 0.9);
        // theSolver2.setPrintValueFunction(false);
        // theSolver2.useErrorBounds(false);
        // theSolver2.setPrintProcessTime(true);
        problem.solve();
        problem.getSolver().setPrintValueFunction(true);
        problem.printSolution();

        PrintWriter pw = new PrintWriter(System.out, true);
        ValueFunction<PendingOrders> probs = problem
                .getSteadyStateProbabilities();
        probs.print(pw, "%-12S", "%10.5f");

        System.out.printf("%-12.3f average orders", problem.averageOrders());
    }
}

class PendingOrders extends PropertiesState {

    // Constructor
    /**
     * @param n
     *            The number of pending orders.
     */
    public PendingOrders(int n) {
        super(new int[] { n });
    }

    @Override
    public String label() {
        return (prop[0] + " Orders");
    }

    /**
     * Return the orders pending
     * 
     * @return orders.
     */
    public int getOrders() {
        return prop[0];

    }

    /**
     * @see jmarkov.basic.State#computeMOPs(MarkovProcess)
     */
    @Override
    public void computeMOPs(MarkovProcess mp) {
    }

    /**
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        return (0 <= prop[0]);
    }
}

// Action
class ProcessAction extends Action {

    private String name;

    /**
     * @param s Name of the action
     */
    public ProcessAction(String s) {
        name = s;
    }

    /**
     * @see jmarkov.basic.Action#label()
     */
    @Override
    public String label() {
        return name;
    }

    /**
     * @see java.lang.Comparable#compareTo(Object)
     */
    public int compareTo(Action a) {
        return name.compareTo(a.label());
    }

}