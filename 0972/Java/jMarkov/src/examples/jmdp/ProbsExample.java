package examples.jmdp;

import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.DecisionRule;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDP;
import jmarkov.jmdp.solvers.ProbabilitySolver;

/**
 * Basic example to show convergence of probabilities with probability solver.
 * @author German Riano. Universidad de los Andes. (C) 2006
 * 
 */
public class ProbsExample extends DTMDP<InvLevel, Order> {

    double[][] probs = { { 0, 0.4, 0.6 }, { 0.8, 0.2, 0 }, { 0.3, 0.5, 0.2 } };

    /**
     * @param initial
     *            Initial inventory level
     */
    public ProbsExample(States<InvLevel> initial) {
        super(initial);
    }

    @Override
    public double immediateCost(InvLevel i, Order a) {
        return 0;
    }

    @Override
    public double prob(InvLevel i, InvLevel j, Order a) {
        int iLevel = i.getLevel();
        int jLevel = j.getLevel();
        return probs[iLevel - 1][jLevel - 1];
    }

    @Override
    public States<InvLevel> reachable(InvLevel i, Order a) {
        StatesSet<InvLevel> set = new StatesSet<InvLevel>();
        int iLevel = i.getLevel();

        if (iLevel == 1) {
            set.add(new InvLevel(2));
            set.add(new InvLevel(3));
        } else if (iLevel == 2) {
            set.add(new InvLevel(1));
            set.add(new InvLevel(2));
        } else if (iLevel == 3) {
            set.add(new InvLevel(1));
            set.add(new InvLevel(2));
            set.add(new InvLevel(3));
        }
        return set;
    }

    @Override
    public Actions<Order> feasibleActions(InvLevel i) {
        return new ActionsSet<Order>(new Order(0));
    }

    /**
     * @param args
     * 
     */
    public static void main(String[] args) {
        States<InvLevel> set = new StatesSet<InvLevel>(new InvLevel(1));
        ProbsExample prob = new ProbsExample(set);
        prob.printSolution();
        DecisionRule<InvLevel, Order> dr;
        try {
            dr = prob.getOptimalPolicy().getDecisionRule();

            ProbabilitySolver<InvLevel, Order> solver = new ProbabilitySolver<InvLevel, Order>(
                    prob, dr);
            solver.solve();
            solver.setGaussSeidel(true);
            System.out.println("Gauss");
            solver.solve();
        } catch (SolverException e) {

            e.printStackTrace();
        }
    }

}
