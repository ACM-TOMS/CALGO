//package DPS;
package jmarkov.basic;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.Map.Entry;

/**
 * This class represents a deterministic decision rule which assigns an action
 * to every state.
 * 
 * @author Andres Sarmiento, Germán Riaño - Universidad de Los Andes
 * @param <S>
 *            The states class
 * @param <A>
 *            Tha Action class
 */
public final class DecisionRule<S extends State, A extends Action> implements
        JMarkovElement, Iterable<Map.Entry<S, A>> {

    private Map<S, A> decisionRule;
    private boolean isStationary = true;

    /**
     * Creates a new empty decision rule
     * 
     */
    public DecisionRule() {
        decisionRule = new TreeMap<S, A>();
    }

    /**
     * Creates a decision rule from a given one
     * 
     * @param dr
     *            decision rule
     */
    public DecisionRule(DecisionRule<S, A> dr) {
        decisionRule = new TreeMap<S, A>();
        for (Map.Entry<S, A> entry : dr)
            decisionRule.put(entry.getKey(), entry.getValue());
    }

    /**
     * Maps a given action to a given state
     * 
     * @param s
     *            state
     * @param a
     *            action
     */
    public void set(S s, A a) {
        decisionRule.put(s, a);
    }

    /**
     * Gets the prescribed action for the given State.
     * 
     * @param s
     *            state
     * @return the action corresponding to the given state
     */
    public A getAction(S s) {
        if (decisionRule.containsKey(s)) {
            return decisionRule.get(s);
        }
        return null;
    }

    /**
     * Returns the amount of states linked to actions in the decision rule.
     * 
     * @return Amount of states linked to actions in the decision rule.
     */
    public int size() {
        return decisionRule.size();
    }

    /**
     * Return an iterator over the State-Action pairs.
     * 
     * @return iterator over the entries
     */
	public Iterator<Map.Entry<S, A>> iterator() {
        return decisionRule.entrySet().iterator();
    }

    @Override
    public String toString() {
        return label();
    }

	public String label() {
        if (size() <= 30)
            return decisionRule.toString();
        return "Decision Rule with " + size() + " states.";
    }

    /**
     * Determines if the given decision rules are equal.
     * 
     * @param o
     * @return True, if the decision rules are equal.
     */
    @Override
    public boolean equals(Object o) {
        if (!(o instanceof DecisionRule))
            return false;
        DecisionRule dr = (DecisionRule)o;
        return decisionRule.equals(dr.decisionRule);
    }

    /**
     * Gives the sting representation of this Rule
     */

	public String description() {
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw);
        print(pw);
        return sw.toString();
    }

    /**
     * Prints the Rule to the sandard output
     */
    public void print() {
        print(new PrintWriter(System.out, true));
    }

    /**
     * Prints the policiy to the given PrintWriter.
     * 
     * @param pw
     *            PrintWriter to use
     */
    public void print(PrintWriter pw) {
        print(pw, "%-10S", "%-10S");
    }

    /**
     * Prints the policiy to the given PrintWriter.
     * 
     * @param pw
     *            PrintWriter to use
     * @param statesFormat
     *            format for the states , for example "%-10S" to have 10 width
     *            left aligned states.
     * @param actionFormat
     *            format for the actions , for example "%-10S" to have 10 width
     *            left aligned actions.
     */
    public void print(PrintWriter pw, String statesFormat, String actionFormat) {

        // pw.println("State ------> Action");
        pw.printf(statesFormat + " ------> " + actionFormat, "State", "Action");
        pw.println();
        Set<Entry<S, A>> entries = decisionRule.entrySet();
        for (Map.Entry<S, A> e : entries) {
            pw.printf(statesFormat + " ------> " + actionFormat, e.getKey()
                    .toString(), e.getValue().toString());
            pw.println();
        }
    }
}
