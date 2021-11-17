package examples.jmdp;

import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.jmdp.FiniteMDP;

/**
 * 
 * A periodic review inventory system that is exposed to random demand. There is
 * no lead time and the warehouse capacity is limited.
 * 
 * The parameters of the problem are: Warehouse capacity: int M Number of
 * stages: int N Fixed cost per order: int K Demand probability function:
 * double[] pr Demand cummulative probability function: double[] qr (q(i) =
 * P{Demand >= i}) Selling cost: double[] F
 * 
 * @author Germán Riaño and Andres Sarmiento - Universidad de Los Andes
 */
public class Inventory extends FiniteMDP<InvLevel, Order> {

    int M;// Capacity
    int lastStage;// lastStage
    int K;// Fixed cost per Order
    double pr[]; // Demand probability function
    double qr[]; // q(i) = P{Demand >= i}
    double F[];
    ActionsSet<Order> actions;
    StatesSet<InvLevel> states;

    double cost(int x) {
        return 2 * x;
    } // variable cost

    double holding(int x) {
        return x;
    } // holding cost

    double orderCost(int x) {
        return (x > 0) ? K + cost(x) : 0.0;
    } // Order cost

    // Constructor

    /**
     * @param lastStage
     *            lastStage
     * @param initSet
     * @param M
     *            Capacity
     * @param K
     *            Fixed cost per Order
     * @param pr
     *            Demand probability function
     * @param qr
     *            q(i) = P{Demand >= i}
     * @param F
     */
    public Inventory(int lastStage, States<InvLevel> initSet, int M, int K,
            double[] pr, double[] qr, double[] F) {
        super(initSet, lastStage);
        this.M = M;// Capacity
        this.lastStage = lastStage;// lastStage
        this.K = K;// Fixed cost per Order
        this.pr = pr; // Demand probability function
        this.qr = qr; // q(i) = P{Demand >= i}
        this.F = F;
        initActions();
    }

    @Override
    public double prob(InvLevel i, InvLevel j, Order a, int t) {
        int iLevel = i.getLevel();
        int jLevel = j.getLevel();
        int orderSize = a.getSize();

        // with stock & demand is positive & order is feasable
        if ((0 < jLevel) && (jLevel <= orderSize + iLevel)
                && (orderSize + iLevel <= M))
            return pr[orderSize + iLevel - jLevel];
        else if ((orderSize + iLevel <= M) && (jLevel == 0)) // End up
            // stockless
            return qr[orderSize + iLevel];
        else
            return 0.0;
    }

    @Override
    public double immediateCost(InvLevel i, Order a, int t) {
        int s = i.getLevel();
        int o = a.getSize();
        return -F[s + o] + orderCost(o) + holding(s + o);
    }

    @Override
    public double finalCost(InvLevel i) {
        return 0.0;
    }

    void initActions() {
        actions = new ActionsSet<Order>();
        for (int k = 0; k <= M; k++) {
            actions.add(new Order(k));
        }
    }

    @Override
    public Actions<Order> feasibleActions(InvLevel i, int t) {
        ActionsSet<Order> actionSet = new ActionsSet<Order>();
        for (Order a : actions) {
            if (a.getSize() <= M - i.getLevel())
                actionSet.add(a);
        }
        return actionSet;
    }

    @Override
    public States<InvLevel> reachable(InvLevel i, Order a, int t) {
        StatesSet<InvLevel> stSet = new StatesSet<InvLevel>();
        for (int n = 0; n <= i.getLevel() + a.getSize(); n++) {
            stSet.add(new InvLevel(n));
        }
        return stSet;
    }

    /**
     * @param a Not used
     * @throws Exception
     */
    public static void main(String a[]) throws Exception {
        int M = 3;// Capacity
        int N = 3;// lastStage
        int K = 4;// Fixed cost per Order
        double pr[] = { 0.25, 0.5, 0.25, 0 }; // Demand probability function
        double qr[] = { 1.0, 0.75, 0.25, 0 }; // q(i) = P{Demand >= i}
        double F[] = { 0, 6, 8, 8 };
        StatesSet<InvLevel> initSet = new StatesSet<InvLevel>(new InvLevel(0));

        Inventory pro = new Inventory(N, initSet, M, K, pr, qr, F);
        pro.solve();
        pro.printSolution();
    }

}// class end
