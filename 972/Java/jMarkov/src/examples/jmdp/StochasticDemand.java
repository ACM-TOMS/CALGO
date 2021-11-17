package examples.jmdp;

import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.jmdp.FiniteMDP;

/**
 * This class belongs to the examples supplied in the package jmdp. The
 * objective of this file is to show as clear as possible a simple way to use
 * the jmdp package as a tool for solving real life problems. The complete
 * details of the present problems are explained in the documentation.
 * 
 * @author Andres Sarmiento, German Riano - Universidad de Los Andes
 */
public class StochasticDemand extends FiniteMDP<InvLevel, Order> {
	//TODO: This example needs more documentation
    int lastStage, maxInventory, maxBackorders, truckSize;

    double K, b, h, theta, price, cost;

    double[] demandProbability, demandCumulativeProbability;

    // demandProbabililty[i] = P[Demand = i]
    // demandCDF[i] = P{Demand >= i}

    // Constructor

    /**
     * @param initSet
     *            Initial level of inventory of the system
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
     * @param theta
     *            demand mean
     */

    public StochasticDemand(States<InvLevel> initSet, int lastStage,
            int maxInventory, int maxBackorders, int truckSize, double K,
            double b, double price, double cost, double h, double theta) {
        super(initSet, lastStage);
        this.maxInventory = maxInventory;
        this.maxBackorders = maxBackorders;
        this.truckSize = truckSize;
        this.K = K;
        this.b = b;
        this.price = price;
        this.cost = cost;
        this.h = h;
        this.theta = theta;
        // initStates();
        initializeProbabilities();
    }

    double holdingCost(int x) {
        double temp = (x > 0) ? h * cost * x : 0.0;
        return temp;
    } // holding cost

    double orderCost(int x) {
        double temp = (x > 0) ? Math.ceil((new Integer(x)).doubleValue()
                / truckSize)
                * K /* + x * cost */: 0.0;
        return temp;
    } // Order cost

    double backorderCost(double x) {
        return (x < 0) ? -b * x : 0.0;
    }

    double lostOrderCost(int x) {
        int mB = maxBackorders;
        double expectedBackorders = 0;
        for (int n = Math.max(x + 1, 0); n <= x + mB; n++)
            expectedBackorders += (n - x) * demandProbability[n];
        double expectedLostDemand = demandCumulativeProbability[x + mB]
                * (theta - x - mB) + (x + mB) * demandProbability[x + mB];
        return (price - cost) * expectedLostDemand
                + backorderCost(-expectedBackorders);
    }

    @Override
    public double finalCost(InvLevel i) {
        return 0.0;
    }

    // see documentation for the explanation for this formula

    @Override
    public double prob(InvLevel i, InvLevel j, Order a, int t) {
        int iLevel = i.getLevel();
        int jLevel = j.getLevel();
        int orderSize = a.getSize();

        // with stock & demand is positive & order is feasable
        if ((-maxBackorders < jLevel) && (jLevel <= orderSize + iLevel)
                && (orderSize + iLevel <= maxInventory))
            return demandProbability[orderSize + iLevel - jLevel];
        else if ((orderSize + iLevel <= maxInventory)
                && (jLevel == -maxBackorders)) // End up stockless
            return demandCumulativeProbability[Math.max(orderSize + iLevel, 0)];
        else
            return 0.0;
    }

    @Override
    public double immediateCost(InvLevel i, Order a, int t) {
        int iLevel = i.getLevel();
        int orderSize = a.getSize();
        double toReturn = orderCost(orderSize)
                + holdingCost(iLevel /* + orderSize */)
                + lostOrderCost(iLevel + orderSize);
        return toReturn;
    }

    void initStates() {
        InvLevel ssts[] = new InvLevel[maxInventory + maxBackorders + 1];
        for (int n = 0; n <= maxInventory; n++) {
            ssts[n] = new InvLevel(n);
        }
        for (int n = maxInventory + 1; n <= maxInventory + maxBackorders; n++) {
            ssts[n] = new InvLevel(n - maxInventory - maxBackorders - 1);
        }
        // states = new StatesCollection<InvLevel>(ssts);
    }

    void initializeProbabilities() {
        demandProbability = new double[maxInventory + maxBackorders + 1];
        demandCumulativeProbability = new double[maxInventory + maxBackorders
                + 1];
        demandProbability[0] = Math.exp(-theta);
        demandCumulativeProbability[0] = 1; // P[demand >= 0]
        double q = 1;
        for (int i = 1; i <= maxInventory + maxBackorders; i++) {
            q = demandCumulativeProbability[i - 1];
            // P{demand >= i}
            demandCumulativeProbability[i] = q - demandProbability[i - 1];
            // P{demand = i}
            demandProbability[i] = demandProbability[i - 1] * theta / i;
        }
    }

    @Override
    public Actions<Order> feasibleActions(InvLevel i, int t) {
        int max = maxInventory - i.getLevel();
        Order[] vec = new Order[max + 1];
        for (int n = 0; n <= max; n++) {
            vec[n] = new Order(n);
        }
        return new ActionsSet<Order>(vec);
    }

    @Override
    public States<InvLevel> reachable(InvLevel i, Order a, int t) {
        StatesSet<InvLevel> statesSet = new StatesSet<InvLevel>();
        for (int n = -maxBackorders; n <= i.getLevel() + a.getSize(); n++) {
            statesSet.add(new InvLevel(n));
        }
        return statesSet;
    }

    /**
     * @param a Not used
     * @throws Exception
     */
    public static void main(String a[]) throws Exception {
        int lastStage = 12;
        int maxInventory = 15;
        int maxBackorders = 5;
        int truckSize = 6;
        int K = 500;
        double b = 1000;
        double h = 0.0050582;// Math.pow(1.3, 1 / 52)-1;
        double theta = 4;
        double price = 22000;
        double cost = 20000;
        InvLevel initial = new InvLevel(0);
        States<InvLevel> initSet = new StatesSet<InvLevel>(initial);

        StochasticDemand pro = new StochasticDemand(initSet, lastStage,
                maxInventory, maxBackorders, truckSize, K, b, price, cost, h,
                theta);
        pro.solve();
        pro.getSolver().setPrintValueFunction(true);
        pro.printSolution();
    }

}