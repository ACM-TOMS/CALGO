/*
 * Created on 19/09/2004
 *
 */
package jmarkov.basic;

import java.io.PrintWriter;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * This structure matches each state with a double number representing its value
 * function, or in some cases the steady state probabilities.
 * 
 * @author Andres Sarmiento and Germán Riaño - Universidad de Los Andes
 * @param <S>
 *            States set.
 */

public class ValueFunction<S extends State> implements JMarkovElement {

    private SortedMap<S, Double> valueFunction;
    private String name = "Value Function";

    /**
     * Creates a new empty value function.
     * 
     */
    public ValueFunction() {
        valueFunction = new TreeMap<S, Double>();
    }

    /**
     * Creates a value function from another given value function
     * 
     * @param vf
     *            value function
     */
    public ValueFunction(ValueFunction<S> vf) {
        valueFunction = new TreeMap<S, Double>(vf.valueFunction);
        name = vf.name;
    }

    /**
     * Creates a new empty value function.
     * 
     * @param name
     *            The name for thei value function
     * 
     */
    public ValueFunction(String name) {
        valueFunction = new TreeMap<S, Double>();
        this.name = name;
    }

    /**
     * Creates a value function from another given value function
     * 
     * @param vf
     *            value function
     * @param name
     *            The name for this value function.
     */
    public ValueFunction(ValueFunction<S> vf, String name) {
        valueFunction = new TreeMap<S, Double>(vf.valueFunction);
        this.name = name;
    }

    /**
     * Associates a state and a double value
     * 
     * @param s
     *            state
     * @param val
     *            value
     */
    public void set(S s, double val) {
        valueFunction.put(s, new Double(val));
    }

    /**
     * Return an iterator used to wakl through the Value Function.
     * 
     * @return iterator over the entries of the map
     */
    public Iterator<Map.Entry<S, Double>> iterator() {
        return valueFunction.entrySet().iterator();
    }

    /**
     * Gets the Value associted with this State.
     * 
     * @param s
     *            given state
     * @return the double value corresponding to the state
     */
    public double get(S s) {
        if (valueFunction.containsKey(s))
            return (valueFunction.get(s)).doubleValue();
        valueFunction.put(s, new Double(0));
        return 0.0;
    }

    /**
     * Gets an array with all the values represented in this value function.
     * 
     * @return an array with the values
     */
    public double[] get() {
        int n = valueFunction.size();
        double[] values = new double[n];
        Collection<Double> vals = valueFunction.values();
        int i = 0;
        for (Double d : vals)
            values[i++] = d.doubleValue();
        return values;
    }

    /**
     * Prints the Value function with the given state format , and values format
     * according to the Format String Syntax.
     * 
     * @param pw
     * @param statesFormat
     *            format for the states , for example "%-10S" to have 10 width
     *            left aligned states.
     * @param valuesFormat
     *            format to use for values. For example us "%6.2" to have 6
     *            width and 2 decimals.
     * @see java.util.Formatter
     */

    public void print(PrintWriter pw, String statesFormat, String valuesFormat) {
        pw.println(name + ":");
        for (Map.Entry<S, Double> entry : valueFunction.entrySet()) {
            pw.printf(statesFormat + " : " + valuesFormat, entry.getKey(),
                    entry.getValue());
            pw.println();
        }
    }

    /**
     * Prints the Value Function. It uses default states and values format.
     * 
     * @param pw
     */
    public void print(PrintWriter pw) {
        print(pw, "%-12S", "%10.2f");
    }

    @Override
    public String toString() {
        return label();
    }

    public String label() {
        int size = valueFunction.size();
        if (size < 20)
            return valueFunction.toString();
        return "Value funciton with " + size + " states";
    }

    public String description() {
        return valueFunction.toString();
    }

}
