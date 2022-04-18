package examples.jmdp;

import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.Events;
import jmarkov.basic.EventsSet;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDPEv;
import jmarkov.jmdp.solvers.ValueIterationSolver;

/**
 * This problem is periodic review, a single item, stochastic demand inventory
 * problem. It is modeled like a discounted cost, infinite horizon, Markov
 * Decision Problem. Demand is assumed to be random according to a Poisson
 * Process distribution with given rate per period. The system is a periodic
 * review problem in which an entity periodically checks the inventory level and
 * takes decisions according to the states he finds. There is a price of selling
 * each item and a cost for buying it. Besides, there is a holding cost incurred
 * when holding one item in stock from one period to another. There is also a
 * fixedCost ordering cost independent of the size of the order placed. The
 * objective is to minimize the expected discounted long run cost.
 * 
 * @author Germán Riaño and Andres Sarmiento - Universidad de Los Andes
 */

public class ControlProduccion extends DTMDPEv<InvLevel, Order, DemandEvent> {
    private int maxInventory;// Capacity
    private double fixedCost;// Fixed cost per Order
    double cost; // variable cost
    double price; // variable cost
    double holdingCost; // holding cost per computer per period.
    // double discountFactor;
    double demPMF[]; // Demand probability function
    double demCCDF[]; // q(i) = P{Demand >= i}
    double demLoss1[]; // E(Demand-i)^+
    double expDemand; // mean of the Poisson demand per stage
    Actions actions;

    // see documentation for the explanation for this forumula

    // Constructor

    /**
     * @param maxInventory
     *            Capacity
     * @param fixedCost
     *            Fixed cost per order
     * @param cost
     *            aquisition cost
     * @param price
     *            selling price
     * @param holdingCost
     *            cost per item stored (non-monetary)
     * @param interestRate
     * @param demandMean
     *            Expected demand
     */
    public ControlProduccion(int maxInventory, double fixedCost, double cost,
            double price, double holdingCost, double interestRate,
            double demandMean) {
        super(new StatesSet<InvLevel>(new InvLevel(0)));
        this.maxInventory = maxInventory;// Capacity
        this.fixedCost = fixedCost;// Fixed cost per Order
        this.cost = cost; // variable cost
        this.price = price; // variable cost
        this.holdingCost = holdingCost; // holding cost per computer per period.
        this.expDemand = demandMean; // mean of the Poisson demand per stage
        demPMF = new double[maxInventory + 1]; // Demand
        // probability
        // function
        demCCDF = new double[maxInventory + 1]; // q(i) = P{Demand >= i}
        demLoss1 = new double[maxInventory + 1]; // E(Demand-i)^+
        init();
        setSolver(new ValueIterationSolver<InvLevel, Order>(this, interestRate));
    }

    double orderCost(int x) {
        return (x > 0) ? fixedCost + cost * x : 0.0;
    } // Order cost

    double lostOrders(int x) {
        // double loss1A = demLoss1[x];
        // double loss1B = ((expDemand - x) * demCCDF[x] + expDemand
        // * demPMF[x]);
        return expDemand * (demCCDF[x] - demPMF[x]) - (x) * demCCDF[x];
    }

    private void initializeProbabilities() {
        double p = Math.exp(-expDemand), q = p;
        demPMF[0] = p;
        demCCDF[0] = 1; // P[demand >= 0]
        demLoss1[0] = expDemand;
        for (int i = 1; i <= maxInventory; i++) {
            p = p * expDemand / i; // P[demand = i]
            demCCDF[i] = 1 - q; // P[demand >= i]
            q += p; // P[demand <= i]
            demPMF[i] = p;
            demLoss1[i] = (expDemand - i) * (1 - q) + expDemand * p;
        }
    }

    @Override
    public States<InvLevel> reachable(InvLevel i, Order a, DemandEvent e) {
        StatesSet<InvLevel> stSet = new StatesSet<InvLevel>();
        if (e.getGreaterThan())
            stSet.add(new InvLevel(0));// minimum level
        else
            stSet.add(new InvLevel(i.getLevel() + a.getSize() - e.getDemand()));
        return stSet;
    }

    @Override
    public Events<DemandEvent> activeEvents(InvLevel i, Order a) {
        EventsSet<DemandEvent> eventSet = new EventsSet<DemandEvent>();
        eventSet.add(new DemandEvent(i.getLevel() + a.getSize(), true));
        for (int n = 0; n < i.getLevel() + a.getSize(); n++) {
            eventSet.add(new DemandEvent(n, false));
        }
        return eventSet;
    }

    @Override
    public double prob(InvLevel i, InvLevel j, Order a, DemandEvent e) {
        return (i.getLevel() + a.getSize() - e.getDemand() == j.getLevel()) ? 1.0
                : 0.0;
    }

    @Override
    public double prob(InvLevel i, DemandEvent e) {
        if (e.getGreaterThan())
            return demCCDF[e.getDemand()];
        return demPMF[e.getDemand()];
    }

    @Override
    public double immediateCost(InvLevel i, Order a, DemandEvent e) {
        int iLevel = i.getLevel();
        int orderSize = a.getSize();
        int demand = e.getDemand();
        double cost = orderCost(orderSize) + holdingCost * iLevel;
        // if(e.getGreaterThan()){
        // cost += lostOrderCost(iLevel + orderSize);
        cost -= demand * price;
        // }else{
        // cost -= demand * price;
        // }
        return cost;
    }

    void init() {
        Order acts[] = new Order[maxInventory + 1];
        InvLevel ssts[] = new InvLevel[maxInventory + 1];
        for (int k = 0; k <= maxInventory; k++) {
            acts[k] = new Order(k);
            ssts[k] = new InvLevel(k);
        }
        states = new StatesSet<InvLevel>(ssts);
        actions = new ActionsSet<Order>(acts);
        initializeProbabilities();
    }

    @Override
    public Actions<Order> feasibleActions(InvLevel i) {
        int max = maxInventory - i.getLevel();
        Order[] vec = new Order[max + 1];
        for (int k = 0; k <= max; k++) {
            vec[k] = new Order(k);
        }
        return new ActionsSet<Order>(vec);
    }

    /**
     * Test Program
     * 
     * @param a
     * @throws SolverException
     */
    public static void main(String a[]) throws SolverException {
        int M = 9;// Capacity
        double K = 50;// Fixed cost per Order
        double cost = 400; // variable cost
        double price = 1000; // variable cost
        double holdingCost = 80; // holding cost per computer per period.
        double discountFactor = 0.9;
        double demandMean = 4; // mean of the Poisson demand per stage

        ControlProduccion prob = new ControlProduccion(M, K, cost, price,
                holdingCost, discountFactor, demandMean);

        prob.getSolver().setPrintValueFunction(true);
        prob.solve();
        prob.printSolution();
    }

}// class end
