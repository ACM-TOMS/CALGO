package examples.jmdp;

import jmarkov.MarkovProcess;
import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.jmdp.DTMDP;
import jmarkov.jmdp.solvers.ValueIterationSolver;

/**
 * 
 * This file models a multi item inventory problem. The demands are random and
 * independent in between them, with Poisson distribution with means "lambda[]".
 * 
 * @author Andres Sarmiento, Germán Riaño (project advisor) - Universidad de Los
 *         Andes
 */
public class LeadTimeStochasticDemand extends DTMDP<LeadTimeState, Order> {

    int maxInventory = 0;// Capacity

    int maxItemsPerOrder = 0;

    int maxBackorders;// Capacity
    int leadTime;// stages for the order to arrive

    double K = 0;// Fixed cost per Order

    double cost; // variable cost

    double price; // variable cost

    double h = 0; // percentual cost of holding inventory

    double holdingCost; // holding cost per period.

    double demandProbability[]; // Demand probability function P[i][t] =

    // P[D(t) = i]

    double demandCumulativeProbability[]; // q(i) Q[i][t] = P[D(t) >= i]

    double theta; // means of the Poisson demands per stage

    double backOrderCost; // cost incurred when there is no

    // stock to satisfy an order.
    // Charged only once.
    ActionsSet<Order> actionsSet;
    StatesSet<LeadTimeState> statesSet;
    ActionsSet<Order> actions;

    // Constructor

    /**
     * @param initSet
     * @param maxInventory
     *            Capacity
     * @param maxItemsPerOrder
     * @param leadTime
     *            stages for the order to arrive
     * @param maxBackorders
     *            Capacity
     * @param K
     *            Fixed cost per Order
     * @param cost
     *            variable cost
     * @param price
     *            variable cost
     * @param h
     *            percentual cost of holding inventory
     * @param holdingCost
     *            holding cost per period.
     * @param theta
     *            means of the Poisson demands per stage
     * @param backOrderCost
     *            cost incurred when there is no
     */
    public LeadTimeStochasticDemand(States<LeadTimeState> initSet,
            int maxInventory, // Capacity
            int maxItemsPerOrder, int leadTime, // number of SKUs
            int maxBackorders,// Capacity
            double K,// Fixed cost per Order
            double cost,// variable cost
            double price,// variable cost
            double h, // percentual cost of holding inventory
            double holdingCost, // holding cost per period.
            double theta,// means of the Poisson demands per stage
            double backOrderCost) // cost incurred when there is no stock to
    // satisfy an order.
    {
        super(initSet);
        this.maxInventory = maxInventory;
        this.maxItemsPerOrder = maxItemsPerOrder;
        this.leadTime = leadTime;
        this.maxBackorders = maxBackorders;
        this.K = K;
        this.cost = cost;
        this.price = price;
        this.h = h;
        this.holdingCost = holdingCost;
        this.theta = theta;
        this.backOrderCost = backOrderCost;
        states = initializeStates();
        System.out.println("States generated");
        initializeActions();
        System.out.println("Actions generated");
        initializeProbabilities();
    }

    double fixedOrderCost(int totalOrderSize) { // Transport cost or some cost
        // associated with the whole order.
        return (totalOrderSize > 0) ? Math.ceil(totalOrderSize
                / maxItemsPerOrder)
                * K : 0.0;
    }

    double variableOrderCost(int itemOrderSize) { // Order
        // cost due to the purchase cost of each item.
        return cost * itemOrderSize;
    }

    double lostOrderCost(int x) { // see documentation for the
        // explanation for this forumula
        int b = maxBackorders;
        double expectedBackorders = 0;
        for (int k = Math.max(x + 1, 0); k <= x + b; k++) {
            expectedBackorders += (k - x) * demandProbability[k];
        }

        // double expectedLostDemand = demandCDF[x + b]
        // * (theta - (x + b)) + (x + b) * demPMF[x + b];

        // lostDemand is such beyond the maxBackorder point
        // double expectedLostDemand = demandCDF[x + b]
        // * (theta - x) + x * demPMF[x + b];

        double toReturn = // (price - cost) * expectedLostDemand
        +backOrderCost * expectedBackorders;
        return toReturn;
    }

    private void initializeProbabilities() {
        demandProbability = new double[maxBackorders + maxInventory + 1];
        demandCumulativeProbability = new double[maxBackorders + maxInventory
                + 1];

        double p = Math.exp(-theta);
        double q = p;
        demandProbability[0] = p;
        demandCumulativeProbability[0] = 1; // P[demand t >= 0]
        for (int i = 1; i <= maxBackorders + maxInventory; i++) {
            p = p * theta / i; // P[demand t = i]
            demandCumulativeProbability[i] = 1 - q; // P[demand t >= i]
            q += p; // P[demand t <= i]
            demandProbability[i] = p;
        }
    }

    /**
     * This function initializes all the possible actions. Starts by creating
     * all the possible orders of size 0, then size 1, etc.
     * 
     */
    void initializeActions() {
        actionsSet = new ActionsSet<Order>();
        for (int size = 0; size <= maxInventory + maxBackorders; size++) {
            actionsSet.add(new Order(size));
        }
        actions = actionsSet;
    }

    /**
     * This function initializes all the possible states. Starts by creating all
     * the possible states with total inventoty in 0 and iterates over the
     * possible total inventories. For each, it iterates growing the number of
     * possible orders and stages pending.
     * 
     */
    StatesSet<LeadTimeState> initializeStates() {
        statesSet = new StatesSet<LeadTimeState>();
        for (int inventory = -maxBackorders; inventory <= maxInventory; inventory++) {
            for (int order = 0; order <= maxInventory - inventory; order++) {
                if (order == 0)
                    statesSet.add(new LeadTimeState(inventory, order, 0));
                else {
                    for (int pendingStages = 0; pendingStages <= leadTime; pendingStages++) {
                        statesSet.add(new LeadTimeState(inventory, order,
                                pendingStages));
                    }
                }
            }
        }
        return new StatesSet<LeadTimeState>(statesSet);
    }

    @Override
    public States<LeadTimeState> reachable(LeadTimeState i, Order a) {
        int orderSize = a.getSize();
        int level = i.getLevel();
        int pendingOrder = i.getOrder();
        int stagesToOrderArrival = i.getStages();
        StatesSet<LeadTimeState> toReturn = new StatesSet<LeadTimeState>();

        if (orderSize == 0) {// no order placed in this period
            if (pendingOrder > 0) { // pending orders
                if (stagesToOrderArrival == 0) {// order arriving
                    for (LeadTimeState s : states) {
                        if (s.getOrder() == 0 && s.getStages() == 0)
                            toReturn.add(s);
                    }
                } else { // no order arriving
                    for (LeadTimeState s : states) {
                        if (s.getLevel() <= level
                                && s.getOrder() == pendingOrder
                                && s.getStages() == stagesToOrderArrival - 1)
                            toReturn.add(s);
                    }
                }
            } else { // no pending orders
                for (LeadTimeState s : states) {
                    if (s.getLevel() <= level && s.getOrder() == 0
                            && s.getStages() == 0)
                        toReturn.add(s);
                }
            }
        } else { // order placed in this period
            for (LeadTimeState s : states) {
                if (s.getLevel() <= level && s.getOrder() == orderSize
                        && s.getStages() == leadTime)
                    toReturn.add(s);
            }
        }
        return toReturn;
    }

    @Override
    public double prob(LeadTimeState i, LeadTimeState j, Order a) {
        int iLevel = i.getLevel();
        int jLevel = j.getLevel();
        int orderSize = a.getSize();
        int inventoryCapacity = maxInventory - iLevel;

        if ((inventoryCapacity >= orderSize) && (jLevel <= orderSize + iLevel)) {
            // feasable order doesnï¿½t exceed the inventory capacity
            int demand = orderSize + iLevel - jLevel;
            if (jLevel == -maxBackorders)
                return demandCumulativeProbability[demand];
            return demandProbability[demand];
        }
        return 0;
    }

    @Override
    public double immediateCost(LeadTimeState i, Order a) {
        int iLevel = i.getLevel();
        int orderSize = a.getSize();
        double totalCost = fixedOrderCost(a.getSize());
        totalCost += lostOrderCost(iLevel + orderSize);
        // totalCost += variableOrderCost(orderSize);
        totalCost += holdingCost * Math.max(iLevel, 0);
        return totalCost;
    }

    @Override
    public Actions<Order> feasibleActions(LeadTimeState i) {
        int iLevel = i.getLevel();
        int pendingOrders = i.getOrder();
        ActionsSet<Order> toReturn = new ActionsSet<Order>();
        for (Order a : actions) {
            if ((pendingOrders > 0) && (a.getSize() == 0)) {
                toReturn.add(a);
                break;
            }// no pending orders
            if (iLevel + a.getSize() <= maxInventory) {
                toReturn.add(a);
            }
        }
        return toReturn;
    }

    /**
     * @param a Not used
     * @throws Exception
     */
    public static void main(String a[]) throws Exception {
        int maxInventory = 15;// Capacity
        int maxItemsPerOrder = 6;
        int leadTime = 2;
        int maxBackorders = 5;// Capacity
        double K = 500;// Fixed cost per Order
        double cost = 20000; // variable cost
        double price = 23000; // variable cost
        double h = Math.pow(1.3, 1 / 52) - 1; // percentual cost of holding
        // inventory
        double holdingCost = h * cost; // holding cost per period.
        double theta = 4; // means of the Poisson demands per stage
        double backOrderCost = 1000; // cost incurred when there is
        // no

        States<LeadTimeState> initSet = new StatesSet<LeadTimeState>(
                new LeadTimeState(0, 0, 0));
        LeadTimeStochasticDemand pro = new LeadTimeStochasticDemand(initSet,
                maxInventory, maxItemsPerOrder, leadTime, maxBackorders, K,
                cost, price, h, holdingCost, theta, backOrderCost);

        ValueIterationSolver<LeadTimeState, Order> theSolver = new ValueIterationSolver<LeadTimeState, Order>(
                pro, 0.9);
        // PolicyIterationSolver<LeadTimeState,Order> theSolver = new
        // PolicyIterationSolver<LeadTimeState,Order>(pro, 0.9);
        // theSolver.setInitVal(-188000);
        // theSolver.setDelta(0.00001);
        // theSolver.useErrorBounds(false);
        theSolver.useGaussSeidel(false);
        theSolver.setPrintValueFunction(true);
        theSolver.solve();
    }

}// class end

class LeadTimeState extends PropertiesState {
    private int level;
    private int pendingOrder;
    private int stagesToOrderArrival;

    LeadTimeState(int level, int pendingOrder, int stagesToOrderArrival) {
        super(3);
        this.prop[0] = level;
        this.prop[1] = pendingOrder;
        this.prop[2] = stagesToOrderArrival;
    }

    /**
     * @return Inventory level
     */
    public int getLevel() {
        return level;
    }

    /**
     * @return Number of pending orders
     */
    public int getOrder() {
        return pendingOrder;
    }

    /**
     * @return Stages to order arrival
     */
    public int getStages() {
        return stagesToOrderArrival;
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