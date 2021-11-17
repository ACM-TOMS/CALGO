package examples.jmdp;

import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.FiniteDP;
import jmarkov.jmdp.solvers.FiniteSolver;

/**
 * This example solves a deterministic, dynamic lot-sizing problem, also known
 * as a Wagner Whitin problem.
 * 
 * @author Andres Sarmiento, Germán Riaño - Universidad de Los Andes
 */

public class WagnerWhitin extends FiniteDP<InvLevel, Order> {
    int lastStage, maxInventory, maxBackorders, truckSize;

    double K, b, h, price, cost;

    int[] demand;

    // Constructor

    /**
     * Crates a dynamic economic lot sizing problem to be solved by Wagner
     * Whitin algorithm.
     * 
     * @param initialInventory
     *            Inventory at time t=0.
     * @param lastStage
     *            the last stage of the problem
     * @param maxInventory
     *            maximum physical capacity in inventory, warehouse size.
     * @param maxBackorders
     *            maximum backorders allowed
     * @param truckSize
     *            maximum items in each fixed cost order. Orders can be greater
     *            than this value, but will be charged more than one fixed cost.
     * @param K
     *            fixed cost per order
     * @param b
     *            unit cost per backordered item per stage
     * @param price
     *            unit price for all stages
     * @param cost
     *            unit costo for all stages
     * @param h
     *            inventory percentual holding cost as a fraction of cost
     * @param demand
     *            demand of items in each stage
     */

    public WagnerWhitin(int initialInventory, int lastStage, int maxInventory,
            int maxBackorders, int truckSize, double K, double b, double price,
            double cost, double h, int[] demand) {
        super(new StatesSet<InvLevel>(new InvLevel(initialInventory)),
                lastStage);
        this.maxInventory = maxInventory;
        this.maxBackorders = maxBackorders;
        this.truckSize = truckSize;
        this.K = K;
        this.b = b;
        this.h = h;
        this.demand = demand;
        this.price = price;
        this.cost = cost;
        init();
    }

    void init() {// This method builds all the states and the actions.
        Order acts[] = new Order[this.maxInventory + maxBackorders + 1];
        InvLevel ssts[] = new InvLevel[maxInventory + maxBackorders + 1];
        for (int k = 0; k < maxInventory + maxBackorders + 1; k++) {
            acts[k] = new Order(k);
            ssts[k] = new InvLevel(k - maxBackorders);
        }
        // states = new StatesCollection<InvLevel>(ssts);
    }

 

    private double holdingCost(int x) {
        return (x > 0) ? h * cost * x : 0.0;
    } // holding cost

    private double orderCost(int x) {
        return (x > 0) ? Math.ceil((double) x / truckSize) * K : 0.0;
    } // Order cost

    double backorderCost(int x) {
        return (x < 0) ? -b * x : 0.0;
    }

    double lostOrderCost(int x, int t) {
        return (x + maxBackorders < demand[t]) ? (price - cost)
                * (demand[t] - x - maxBackorders) : 0.0;
    }

    /**
     * Returns the optimal cost for this level of starting inventory.
     * 
     * @param inventory
     * @return The optimal cost for this level of starting inventory.
     * @throws SolverException 
     */
    public double getOptimalCost(int inventory) throws SolverException {
        return getOptimalValueFunction().get(new InvLevel(inventory));
    }

    @Override
    public double immediateCost(InvLevel i, Order a, int t) {
        int s = i.getLevel();
        int o = a.getSize();
        return lostOrderCost(o, t) + orderCost(o) + holdingCost(s + o)
                + backorderCost(s + o);
        // return -2000*(Math.max(s+o-demand[t],0))+holding(s + o, t)+
        // orderCost(o, t);
    }

    @Override
    public double finalCost(InvLevel i) {
        return 0.0;
    }

    @Override
    public Actions<Order> feasibleActions(InvLevel i, int t) {
        ActionsSet<Order> actionSet = new ActionsSet<Order>();
        int min_order = Math.max(-maxBackorders - i.getLevel() + demand[t], 0);
        int max_order = maxInventory - i.getLevel() + demand[t];
        for (int n = min_order; n <= max_order; n++) {
            actionSet.add(new Order(n));
        }
        return actionSet;
    }

    @Override
    public InvLevel destination(InvLevel i, Order a, int t) {
        int o = a.getSize();
        int iLevel = i.getLevel();
        return new InvLevel(Math.max(iLevel + o - demand[t], -maxBackorders));
    }

    /**
     * Test Program.
     * 
     * @param a
     * @throws Exception
     */
    public static void main(String a[]) throws Exception {
        int lastStage = 12;
        int maxInventory = 15;
        int maxBackorders = 5;
        int truckSize = 6;
        double K = 500;
        double b = 2000;
        double p = 22000;
        double c = 20000;
        double h = Math.pow(1.3, 1.0 / 52) - 1.0;
        int[] demand = new int[] { 10, 4, 3, 6, 3, 2, 0, 1, 7, 3, 4, 5 };

        WagnerWhitin prob = new WagnerWhitin(0, lastStage, maxInventory,
                maxBackorders, truckSize, K, b, p, c, h, demand);
         

        FiniteSolver<InvLevel, Order> theSolver = new FiniteSolver<InvLevel, Order>(
                prob);
        prob.setSolver(theSolver);
        prob.solve();
        prob.getSolver().setPrintValueFunction(true);
        // System.out.println(theSolver.bestPolicy(initial));
        prob.printSolution();
        prob.getOptimalCost(0);
    }

}