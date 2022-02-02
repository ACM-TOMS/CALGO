package jmarkov.basic;

import java.io.PrintWriter;
import java.io.StringWriter;

/**
 * Policy is a set of "Decision Rules". It contains the Decision Rule for every
 * stage t. When the problem has infinite horizon, only one decision rule is
 * stored.
 * 
 * @author Andres Sarmiento, Germán Riaño
 * @param <S>
 *            The states class
 * @param <A>
 *            Tha Action class
 */

public final class Policy<S extends State, A extends Action> {

    DecisionRule<S, A> decisionRulesP[];
    boolean stationary = false;

    /**
     * Creates a set with the given horizon.
     * 
     * @param stages
     *            The number of stages
     */

    @SuppressWarnings("unchecked")
    public Policy(int stages) {
        decisionRulesP = new DecisionRule[stages];
        stationary = false;
    }

    /**
     * Creates a stationary policy with the given decition rule
     * 
     * @param d
     *            The rule
     */
    @SuppressWarnings("unchecked")
    public Policy(DecisionRule d) {
        decisionRulesP = new DecisionRule[1];
        setDecisionRule(d);
        stationary = true;
    }

    /**
     * Return the time horizon for this Ploicy.
     * @return last stage where actions can be taken
     */
    public int getHorizon() {
        if (stationary)
            return Integer.MAX_VALUE;
        return decisionRulesP.length;
    }

    /*
     * Sets this policy.
     */

    /**
     * Sets a decision rule for stage t in the policy
     * 
     * @param dr
     *            decision rule
     * @param t
     *            stage
     */
    public void setDecisionRule(DecisionRule<S, A> dr, int t) {
        decisionRulesP[t] = dr;
    }

    /**
     * Sets a unique decision rule for the policy, for infinite horizon
     * problems.
     * 
     * @param pol
     */
    public void setDecisionRule(DecisionRule<S, A> pol) {
        decisionRulesP[0] = pol;
        stationary = true;
    }

    /**
     * 
     * @return the unique decision rule, for infinite horizon problems.
     */
    public DecisionRule<S, A> getDecisionRule() {
        assert(stationary);
        if (decisionRulesP.length > 1)
            return null;
        return decisionRulesP[0];
    }

    /**
     * Returns the decision rule for statge t
     * 
     * @param t
     *            stage.
     * @return The decision rule for stage t.
     */
    public DecisionRule<S, A> getDecisionRule(int t) {
        return decisionRulesP[t];
    }

   

    // /**
    // * Sets an action to a stage t/state.
    // */
    //
    // public void set(S i, A a, int t) {
    // DecisionRule<S, A> dr = decisionRule(t);
    // dr.set(i, a);
    // }

    /**
     * Gets the action to be taken in state i at this stage t
     * 
     * @param i
     *            The state
     * @param t
     *            The stage (time) at which action is to be taken.
     * @return The action.
     */

    public A getAction(S i, int t) {
        DecisionRule<S, A> decisionR = getDecisionRule(t);
        return decisionR.getAction(i);
    }

    /**
     * Gives the sting representation of this Policy
     */

    @Override
    public String toString() {
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw);
        print(pw);
        return sw.toString();
    }

    /**
     * Prints the policy to the standard output
     */
    public void print() {
        print(new PrintWriter(System.out, true));
    }

    /**
     * Prints the policy to the given PrintWriter.
     * 
     * @param pw
     *            print writer
     */
    public void print(PrintWriter pw) {
        pw.println("********* Best Policy *********\n");
        int n = decisionRulesP.length;
        for (int t = 0; t < n; t++) {
            if (!stationary) {
                pw.println("Stage: " + t);
            } else {
                pw.println("In every stage do: ");
            }
            getDecisionRule(t).print(pw);
        }
    }

}
