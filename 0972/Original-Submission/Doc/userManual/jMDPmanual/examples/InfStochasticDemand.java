package examples.jmdp;

import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDP;
import jmarkov.jmdp.solvers.PolicyIterationSolver;
import jmarkov.jmdp.solvers.RelativeValueIterationSolver;
import jmarkov.jmdp.solvers.ValueIterationSolver;
import Jama.Matrix;

/**
 * This problem is a single item, periodic review, stochastic demand inventory
 * problem. It is modeled like a discounted cost, infinite horizon, Markov
 * Decision Problem. Demand is assumed to be random according to a Poisson
 * Process distribution with given rate per period. The system is a periodic
 * review problem in which an entity periodically checks the inventory level and
 * takes decisions according to the states he finds. There is a price of selling
 * each item and a cost for buying it. Besides, there is a holding cost incurred
 * when holding one item in stock from one period to another. There is also a
 * truckCost ordering cost independent of the size of the order placed. The
 * objective is to minimize the expected discounted long run cost.
 * 
 * @author Germán Riaño and Andres Sarmiento - Universidad de Los Andes
 */
public class InfStochasticDemand extends DTMDP<InvLevel, Order> {
	//TODO This example needs more specific documentation.
    // Problem parameters:
    private int maxInv, maxBO, truckSize;
    // Cost and demand parameters:
    private double truckCost, backorderCost, holdingCost, intRate, expDemand,
            price, cost;
    private double[] demPMF, demCDF, demandLoss1;
    private boolean isdisc = false;

    // demPMF[i] = P{Demand = i}
    // demCDF[i] = P{Demand <= i}
    // demandLoss1[i] = E[(Demand - i)^+ ]
    // Constructor

    /**
     * @param maxInv
     *            maximum physical capacity in inventory, warehouse size.
     * @param maxBO
     *            maximum backorders allowed
     * @param truckSize
     *            maximum items in each fixed cost order. Orders can be greater
     *            than this value, but will be charged more than one fixed cost.
     * @param truckCost
     *            fixed cost per order
     * @param backorderCost
     *            unit cost per backordered item per stage
     * @param price
     *            unit price
     * @param cost
     *            unit aquistion costo
     * @param holdingCost
     *            non-finantial holding cost (it does NOT include finantial
     *            cost)
     * @param intRate
     *            interest per period
     * @param expDemand
     *            demand mean
     * @param discounted
     *            Whether a discounted model (rather than average) is to be
     *            used.
     */

    @SuppressWarnings("unchecked")
    public InfStochasticDemand(int maxInv, int maxBO, int truckSize,
            double truckCost, double backorderCost, double price, double cost,
            double holdingCost, double intRate, double expDemand,
            boolean discounted) {
        super(new StatesSet<InvLevel>(new InvLevel(0)));
        this.maxInv = maxInv;
        this.maxBO = maxBO;
        this.truckSize = truckSize;
        this.truckCost = truckCost;
        this.backorderCost = backorderCost;
        this.price = price;
        this.cost = cost;
        this.holdingCost = holdingCost;
        this.expDemand = expDemand;
        // initStates();
        initializeProbabilities();
        this.isdisc = discounted;
        this.intRate = intRate;
        if (discounted)
            setSolver(new ValueIterationSolver(this, intRate));
        else
            setSolver(new RelativeValueIterationSolver(this));
    }

    private void initializeProbabilities() {
        demPMF = new double[maxInv + maxBO + 1];
        demCDF = new double[maxInv + maxBO + 1];
        demandLoss1 = new double[maxInv + maxBO + 1];
        double cdf, p = Math.exp(-expDemand);
        cdf = demCDF[0] = demPMF[0] = p;
        demandLoss1[0] = expDemand;
        int maxlevel = maxInv + maxBO;
        for (int i = 1; i <= maxlevel; i++) {
            demPMF[i] = (p *= expDemand / i); // P{demand = i}
            demCDF[i] = (cdf += p); // P{demand <= i}
            demandLoss1[i] = (expDemand - i) * (1 - cdf) + expDemand * p;
            // = E[(D-i)^+]
        }
    }

    @Override
    public States<InvLevel> reachable(InvLevel i, Order a) {
        StatesSet<InvLevel> statesSet = new StatesSet<InvLevel>();
        // Available inventory upon order receival:
        int maxLevel = i.getLevel() + a.getSize();
        for (int n = -maxBO; n <= maxLevel; n++) {
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
        int demand = orderSize + iLevel - jLevel;
        assert (demand >= 0);
        try {
            if (jLevel == -maxBO)
                return 1.0 - ((demand > 0) ? demCDF[demand - 1] : 0.0);
            else
                // End up stockless
                return demPMF[demand];
        } catch (IndexOutOfBoundsException e) {
            throw new IllegalArgumentException(
                    "'prob' called on non-reachable state!!. i=" + iLevel
                            + ", j=" + jLevel + ", a =" + orderSize, e);
        }
    }

    @Override
    public Actions<Order> feasibleActions(InvLevel i) {
        int max = maxInv - i.getLevel();
        Order[] vec = new Order[max + 1];
        for (int n = 0; n <= max; n++) {
            vec[n] = new Order(n);
        }
        return new ActionsSet<Order>(vec);
    }

    double holdingCost(int x) {
        double totHoldCost = holdingCost + ((isdisc) ? intRate * cost : 0.0);
        return (x > 0) ? totHoldCost * x : 0.0;
    } // holding cost

    double orderCost(int x) {
        return truckCost * Math.ceil((double) x / truckSize) + x * cost;
    } // Order cost

    double backorderCost(double x) {
        return (x < 0) ? -backorderCost * x : 0.0;
    }

    // see documentation for the explanation for this formula
    @Override
    public double immediateCost(InvLevel i, Order a) {
        int maxSale = i.getLevel() + a.getSize() + maxBO;
        double expectedSales = expDemand - demandLoss1[maxSale];
        double netProfit = price * expectedSales - orderCost(a.getSize())
                - holdingCost(i.getLevel()) - backorderCost(i.getLevel());
        return -netProfit;
    }

    /**
     * Very stupid method to see what this is doing!!
     */
    public void printMatrices() {
        double[][] cost = new double[maxBO + maxInv + 1][maxBO + maxInv + 1];
        double[][][] prb = new double[maxBO+maxInv + 1][maxBO+maxInv + 1][maxBO+maxInv + 1];
        for (InvLevel s: getAllStates()) {
            int i = s.getLevel();
            for (Order o: feasibleActions(s)) {
                int a = o.getSize();
                cost[i+maxBO][a] = immediateCost(new InvLevel(i), new Order(a));
                for (InvLevel y:reachable(s,o)) {
                    int j = y.getLevel();
                    prb[a][i+maxBO][j+maxBO] = prob(new InvLevel(i), new InvLevel(j),
                            new Order(a));
                }
            }
        }
        (new Matrix(cost)).print(8, 2);
        for (int a = 0; a < maxInv; a++) {
            (new Matrix(prb[a])).print(10, 6);
        }
        (new Matrix(new double[][] { demPMF })).print(10, 6);
        (new Matrix(new double[][] { demCDF })).print(10, 6);
        (new Matrix(new double[][] { demandLoss1 })).print(10, 6);
    }

    /**
     * Simple test Program.
     * 
     * @param a
     * @throws SolverException
     */
    public static void main(String a[]) throws SolverException {
        int maxInventory = 25;
        int maxBackorders = 0;
        int truckSize = 4;
        int truckCost = 1000;
        double b = 0;// 1000;
        double holdCost = 50;
        double intRate = Math.pow(1.3, 1 / 52);
        double theta = 20;
        double price = 1100;// 22000;
        double cost = 500;// 20000;

        InfStochasticDemand prob = new InfStochasticDemand(maxInventory,
                maxBackorders, truckSize, truckCost, b, price, cost, holdCost, intRate, theta,
                false);

        RelativeValueIterationSolver<InvLevel, Order> solv = new RelativeValueIterationSolver<InvLevel, Order>(
                prob);

        prob.setSolver(solv);
        prob.getSolver().setPrintValueFunction(true);
        prob.solve();
        prob.printSolution();


    }

}