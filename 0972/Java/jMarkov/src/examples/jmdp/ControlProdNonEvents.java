package examples.jmdp;

import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.DecisionRule;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDP;
import jmarkov.jmdp.solvers.PolicyIterationSolverAvg;
import jmarkov.jmdp.solvers.RelativeValueIterationSolver;
import jmarkov.jmdp.solvers.Solver;
import jmarkov.jmdp.solvers.ValueIterationSolver;
import Jama.Matrix;

/**
 * The problem is a periodic review, single item, stochastic demand inventory
 * problem. It is modeled like a discounted cost, infinite horizon, Markov
 * Decision Problem. Demand is assumed to be random according to a Poisson
 * distribution with mean "lambda". The system is a periodic review problem in
 * which an entity periodically checks the inventory level and takes decisions
 * according to the states he finds. There is a price of selling each item and a
 * cost for buying it. Besides, there is a holding cost incurred when holding
 * one item in stock from one period to another. There is also a fixed ordering
 * cost independent of the size of the order placed. The objective is to
 * minimize the expected discounted long run cost.
 * 
 * @author Germán Riaño and Andres Sarmiento - Universidad de Los Andes
 */

public class ControlProdNonEvents extends DTMDP<InvLevel, Order> {
    private int maxInventory;// Capacity
    private double fixedCost,// Fixed cost per Order
            cost, // variable cost
            price, // variable cost
            holdingCost, // holding cost per computer per period.
            demPMF[], // Demand probability function
            demCCDF[], // q(i) = P{Demand >= i}
            demandLoss1[], // E(Demand-i)^+
            expDemand, // mean of the Poisson demand per stage
            interestRate; // interest rate per period.
    private Actions actions;

    // Constructor
    /**
     * @param maxInventory
     *            maximum physical capacity in inventory, warehouse size. than
     *            this value, but will be charged more than one fixed cost.
     * @param fixedCost
     *            fixed cost per order
     * @param price
     *            unit price for all stages
     * @param cost
     *            unit costo for all stages
     * @param holdingCost
     *            holding cost per unit and per unit of time.
     * @param interestRate
     *            interest rate per period.
     * @param expDemand
     *            demand mean
     */

    public ControlProdNonEvents(int maxInventory, double fixedCost,
            double cost, double price, double holdingCost, double interestRate,
            double expDemand) {
        super(new InvLevel(0));
        this.maxInventory = maxInventory;// Capacity
        this.fixedCost = fixedCost;// Fixed cost per Order
        this.cost = cost; // variable cost
        this.price = price; // variable cost
        this.holdingCost = holdingCost; // holding cost per computer per period.
        this.expDemand = expDemand; // mean of the Poisson demand per stage
        this.interestRate = interestRate;
        demPMF = new double[maxInventory + 1]; // Demand
        // probability
        // function
        demCCDF = new double[maxInventory + 1]; // q(i) = P{Demand >= i}
        demandLoss1 = new double[maxInventory + 1]; // E(Demand-i)^+
        init();
    }

    /**
     * Gets the optimal value function for this initial inventory level.
     * 
     * @param invLevel
     *            inventory level.
     * @return The optimal value fo future costs.
     * @throws SolverException
     */
    public double getValueFunction(int invLevel) throws SolverException {
        return getOptimalValueFunction().get(new InvLevel(invLevel));
    }

    /**
     * Return the optimal order size for this inventory level.
     * 
     * @param invLevel
     * @return order size.
     */
    public int getOptimalOrderSize(int invLevel) {
        try {
            return getOptimalPolicy().getDecisionRule().getAction(
                    new InvLevel(invLevel)).getSize();
        } catch (SolverException e) {
            return -1;
        }
    }

    /**
     * Return the optimal order size for this inventory level.
     * 
     * @return order size.
     */
    public int[] getOptimalOrderSize() {
        States<InvLevel> stts = getAllStates();
        DecisionRule<InvLevel, Order> dr;
        int n = stts.size(), acts[] = new int[n];
        try {
            dr = getOptimalPolicy().getDecisionRule();
            for (int i = 0; i < n; i++)
                acts[i] = dr.getAction(new InvLevel(i)).getSize();
        } catch (SolverException e) {
        }
        return acts;
    }

    private double orderCost(int x) {
        return (x > 0) ? fixedCost + cost * x : 0.0;
    } // Order cost

    private void initializeProbabilities() {
        double p = Math.exp(-expDemand), cdf = p;
        demPMF[0] = p;
        demCCDF[0] = 1.0; // P{demand >= 0}
        demandLoss1[0] = expDemand;
        for (int i = 1; i <= maxInventory; i++) {
            demCCDF[i] = 1 - cdf; // P{demand >= i}
            p = p * expDemand / i; // P{demand = i}
            cdf += p; // P{demand <= i}
            demPMF[i] = p;
            demandLoss1[i] = (expDemand - i) * (1 - cdf) + expDemand * p;
        }
    }

    @Override
    public States<InvLevel> reachable(InvLevel i, Order a) {
        StatesSet<InvLevel> statesSet = new StatesSet<InvLevel>();
        int maxInv = i.getLevel() + a.getSize();
        for (int n = 0; n <= maxInv; n++) {
            statesSet.add(new InvLevel(n));
        }
        return statesSet;
    }

    @Override
    public double prob(InvLevel i, InvLevel j, Order a) {
        int iLevel = i.getLevel();
        int jLevel = j.getLevel();
        int orderSize = a.getSize();

        // with stock & demand is positive & order is feasable
        if ((0 < jLevel) && (jLevel <= orderSize + iLevel)
                && (orderSize + iLevel <= maxInventory))
            return demPMF[orderSize + iLevel - jLevel];
        else if ((orderSize + iLevel <= maxInventory) && (jLevel == 0))
            // End up stockless
            return demCCDF[Math.max(orderSize + iLevel, 0)];
        else
            throw new IllegalArgumentException(
                    "'prob' Called on non-reachable state: i=" + iLevel + ",j="
                            + jLevel + " ,a=" + orderSize);
    }

    @Override
    public double immediateCost(InvLevel i, Order a) {
        int iLevel = i.getLevel();
        int orderSize = a.getSize();
        int available = iLevel + orderSize;
        if (available > maxInventory)
            throw new IllegalArgumentException(
                    "cost called on unavailable action!");
        double expectedSales = expDemand - demandLoss1[available];
        double netProfit = price * expectedSales - orderCost(a.getSize())
                - holdingCost * i.getLevel();
        return -netProfit;

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
     * Very stupid method to see what this is doing!!
     */
    public void printMatrices() {
        double[][] cost = new double[maxInventory + 1][maxInventory + 1];
        double[][][] prb = new double[maxInventory + 1][maxInventory + 1][maxInventory + 1];
        for (int i = 0; i <= maxInventory; i++) {
            for (int a = 0; a + i <= maxInventory; a++) {
                cost[i][a] = immediateCost(new InvLevel(i), new Order(a));
                int invMax = i + a;
                for (int j = 0; j <= invMax; j++) {
                    prb[a][i][j] = prob(new InvLevel(i), new InvLevel(j),
                            new Order(a));
                }
            }
        }
        (new Matrix(cost)).print(8, 2);
        for (int a = 0; a < maxInventory; a++) {
            (new Matrix(prb[a])).print(10, 6);
        }
        (new Matrix(new double[][] { demPMF })).print(10, 6);
        (new Matrix(new double[][] { demCCDF })).print(10, 6);
        (new Matrix(new double[][] { demandLoss1 })).print(10, 6);
    }

    /**
     * Small test program
     * 
     * @param a
     * @throws SolverException
     */
    public static void main(String a[]) throws SolverException {
        int M = 4;// Capacity
        double K = 500;// Fixed cost per Order
        double cost = 400; // variable cost
        double price = 1000; // variable cost
        double holdingCost = 80; // holding cost per computer per period.
        double interestRate = 0.10;
        double demandMean = 4; // mean of the Poisson demand per stage
        ControlProdNonEvents prob = new ControlProdNonEvents(M, K, cost, price,
                holdingCost, interestRate, demandMean);

        prob.printMatrices();

        PolicyIterationSolverAvg<InvLevel, Order> solv = new PolicyIterationSolverAvg<InvLevel, Order>(prob);
        //RelativeValueIterationSolver<InvLevel, Order> solv = new RelativeValueIterationSolver<InvLevel, Order>(prob);
        
        prob.setSolver(solv);
         solv.setPrintGain(true);
         prob.setDebugLevel(4);
        prob.getSolver().setPrintValueFunction(true);
        prob.solve();
        prob.printSolution();
        // System.out.println(solv.getProcessTime());
    }

}// class end

