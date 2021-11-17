/*
 * Created on 8/09/2005
 */
package examples.jmdp;

import jmarkov.MarkovProcess;
import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.CTMDP;
import jmarkov.jmdp.solvers.RelativeValueIterationSolver;

/**
 * This class describes an inventory system with two items. Item 1 and item 2
 * have different demand rates. Ordering costs depend on trucks with fixed
 * capacity that may transport any combination of the items. No leadtime.
 * 
 * @author Andres Sarmiento and Germán Riaño. Universidad de los Andes.
 */
public class CTInventory extends CTMDP<CTStock, Order> {

    int maxCapacity, truckSize;
    double K, cost, price, holdingCost, demandProbability[], demandCCDF[],
            demandLoss1[], demandMean, leadTime;
    Actions actions;

    /**
     * @param initSet
     *            Initial state
     * @param maxCapacity
     *            maximum physical capacity in inventory, warehouse size
     * @param truckSize
     *            maximum physical capacity of truck to carrying items
     * @param K
     *            Fixed cost per Order
     * @param cost
     *            Variable cost
     * @param price
     *            Sell price of item A
     * @param holdingCost
     *            Inventory holding cost per period
     * @param demandRateA
     *            Demand rate for item A
     * @param leadTime
     *            Time between the moment when the order is placed and the
     *            moment the order arrives
     */
    public CTInventory(States<CTStock> initSet, int maxCapacity, int truckSize,
            double K, double cost, double price, double holdingCost,
            double demandRateA, double leadTime) {
        super(initSet);
        this.maxCapacity = maxCapacity;// Capacity
        this.truckSize = truckSize;// truck capacity
        this.K = K;// Fixed cost per Order
        this.cost = cost; // variable cost
        this.price = price; // variable cost
        this.holdingCost = holdingCost; // holding cost per computer per period.
        this.demandMean = demandRateA; // demand rate for item A
        this.leadTime = leadTime;
        demandProbability = new double[maxCapacity + 1]; // Demand
        // probability
        // function
        demandCCDF = new double[maxCapacity + 1]; // q(i) = P{Demand >= i}
        demandLoss1 = new double[maxCapacity + 1]; // E(Demand-i)^+
        initializeProbabilities();
    }

    //
    // Maintainment functions
    //

    /**
     * Initialize the probabilities
     */
    public void initializeProbabilities() {
        double p = Math.exp(-demandMean * leadTime), q = p;
        demandProbability[0] = p;
        demandCCDF[0] = 1; // P[demand >= 0]
        demandLoss1[0] = demandMean * leadTime;
        for (int i = 1; i <= maxCapacity; i++) {
            p = p * demandMean * leadTime / i; // P[demand = i]
            demandCCDF[i] = 1 - q; // P[demand >= i]
            q += p; // P[demand <= i]
            demandProbability[i] = p;
            demandLoss1[i] = (demandMean * leadTime - i) * (1 - q) + demandMean
                    * leadTime * p;
        }
    }

    //
    // Inherited methods
    //

    @Override
    public Actions<Order> feasibleActions(CTStock i) {
        ActionsSet<Order> set = new ActionsSet<Order>();
        int iLevel = i.getItems();
        if (iLevel + i.getOrders() < maxCapacity)
            set.add(new Order(1));
        if (iLevel > 1 || i.getOrders() > 0)
            set.add(new Order(0));
        return set;
    }

    @Override
    public States<CTStock> reachable(CTStock i, Order a) {
        int items = i.getItems();
        int orders = i.getOrders();
        StatesSet<CTStock> statesSet = new StatesSet<CTStock>();
        // demand arrival
        if (items > 0)
            statesSet.add(new CTStock(items - 1, orders + a.getSize()));
        // order arrival
        if (orders > 0)
            statesSet.add(new CTStock(items + 1, orders - 1 + a.getSize()));
        return statesSet;
    }

    @Override
    public double rate(CTStock i, CTStock j, Order a) {
        if (j.getItems() < i.getItems())// demand arrival
            return demandMean;
        return i.getOrders() / leadTime; // order arrival. min among the a
        // orders placed.
    }

    @Override
    public double lumpCost(CTStock i, Order a) {
        int x = i.getItems() + i.getOrders();
        double lostOrders = demandMean * leadTime
                * (demandCCDF[x] - demandProbability[x]) - (x) * demandCCDF[x];
        double orderCost = (x > 0) ? K : 0;// Math.ceil(((double)x) /
        // truckSize) * K : 0.0;
        return orderCost + (price - cost) * lostOrders;
    }

    @Override
    public double continuousCost(CTStock i, Order a) {
        return holdingCost * i.getItems();
    }

    /**
     * This method just tests the class.
     * 
     * @param args
     *            Not used
     * @throws SolverException
     */
    public static void main(String[] args) throws SolverException {
        int maxCapacity = 15;// Capacity
        int truckSize = 6; // truck capacity
        double K = 0;// Fixed cost per Order
        double cost = 400; // variable cost
        double price = 1000; // variable cost
        double holdingCost = 80; // holding cost per computer per period.
        double demandMean = 4; // mean of the Poisson demand per unit time
        double leadTime = 2; // expected time units for order arrival
        States<CTStock> initSet = new StatesSet<CTStock>(new CTStock(1, 0));

        CTInventory prob = new CTInventory(initSet, maxCapacity, truckSize, K,
                cost, price, holdingCost, demandMean, leadTime);

        RelativeValueIterationSolver<CTStock, Order> solv = new RelativeValueIterationSolver<CTStock, Order>(
                prob);
        solv.setPrintValueFunction(true);
        prob.solve();
        prob.printSolution();
        System.out.println("time:" + solv.getProcessTime());
    }

}

class CTStock extends PropertiesState {

    /**
     * @param a
     *            Physical stock in inventory
     * @param b
     *            Amount of orders placed
     */
    public CTStock(int a, int b) {
        super(new int[] { a, b });
    }

    /**
     * @return Physical stock in inventory
     */
    public int getItems() {// physical stock in inventory
        return prop[0];
    }

    /**
     * @return Amount of orders placed
     */
    public int getOrders() {// amount of orders placed
        return prop[1];
    }

    @Override
    public String label() {
        return "Stock:" + getItems() + "-Orders:" + getOrders();
    }

    /**
     * @see jmarkov.basic.State#computeMOPs(MarkovProcess)
     */
    @Override
    public void computeMOPs(MarkovProcess mp) {
        // TODO Auto-generated method stub
    }

    /**
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        // TODO Auto-generated method stub
        return true;
    }
}

// class CTOrder extends PropertiesAction{
//    
// public CTOrder(int a, int b) {
// super(new int[] {a,b});
// }
//    
// public int getItemsA(){
// return properties[0];
// }
//    
// public int getItemsB(){
// return properties[1];
// }
//    
// @Override
// public String label() {
// return "Order:(" + getItemsA() + "," + getItemsB() + ")";
// }
// }
