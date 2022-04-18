package jmarkov.basic;

import java.util.ArrayList;
import java.util.List;

import jmarkov.MarkovProcess;

/**
 * The Class State represent a state in a MarkovProcess or MDP. The user of the
 * class should estiblish her own coding convention AND code the compareTo
 * method. If the State can be represented with a vector of integers describing
 * its properties, then it might be easier to implement PropertiesState rather
 * than State.
 * 
 * @see PropertiesState PropertiesState class
 * @author German Riaño. Universidad de los Andes.
 * @version 1.0a
 */
public abstract class State implements Comparable<State>, JMarkovElement {

    private List<Double> mops = null;

    private int index = -1;


    /**
     * The method compareTo should be implemented in order to establish a total
     * ordering among the States.
     * 
     * @return A positive integer if this is grater then j, negative if this is
     *         less then j and 0 if this == j.
     * @see Comparable#compareTo(Object)
     */

    public abstract int compareTo(State j);

    // This ensure that the ordering is consistent with equals.
    /**
     * If Object is not State it returns false. Otherwise equals :=
     * (compareTo(o)==0)
     * 
     * @see java.lang.Object#equals(java.lang.Object)
     */
    @Override
    public final boolean equals(Object o) {
        if (!(o instanceof State))
            return false;
        State s1 = (State) o;
        return (compareTo(s1) == 0);
    }

    // MOPS stuff //
   
    /**
     * Sets the value of this MOP.
     * 
     * @param index
     * @param value
     * @return the index where it was added.
     */
    public final int setMOP(int index, double value) {
        if (mops == null)
            mops = new ArrayList<Double>(10);
        while (index >= mops.size()) {
            // make sure array is big enough
            mops.add(null);
        }
        mops.set(index, new Double(value));// set the mop
        return index;
    }

    /**
     * Sets the value of the MOP with this name. If no MOP with this name exists
     * a new one is declared.
     * 
     * @param mopName
     * @param model
     *            The model being solved.
     * @param value
     * @return the index where it was added.
     */
    public int setMOP(MarkovProcess<?, ?> model, String mopName, double value) {
        int i = model.getMOPIndex(mopName);
        if (i == -1) {
            model.addMOP(mopName);
            i = model.numMOPs() - 1; // last position
        }
        setMOP(i, value);
        return i;
    }

    /**
     * This method should be implemented in order to compute all the measures of
     * performance MOPs. Inside it you should issue commands like
     * <code>setMop("Utilization server 1", x, model);</code>. * For large
     * models override this method as empty and rather override getMOP(int). Do
     * NOT mix both approaches!!
     * 
     * @param model
     *            The model being solved.
     * @see #getMOP(int)
     */
    public abstract void computeMOPs(MarkovProcess<?, ?> model);

    /**
     * Gets the value of the MOP with this name, by calling
     * <code>getMOP(int)<code>
     * 
     * @param mopName The name of the MOP.
     * @param model Model being solved.
     * @return current MOP value
     * @see #getMOP(int)
     */
    public final double getMOP(String mopName, MarkovProcess<?, ?> model) {
        int i = model.getMOPIndex(mopName);
        return getMOP(i);
    }

    /**
     * Gets the value of this MOP. The value should had been set via the setMOP
     * method. Alternatively, for better performance define the MOP Names when
     * implementing the MarkovProcess class and override this method. To define
     * the names in the constructor call the method setMOPs(String[]). The index
     * is the same as the one used in the array in the aforementioned method.
     * 
     * @param index
     * @return The value of this MOP.
     * @see MarkovProcess#setMOPs(String[])
     */
    public double getMOP(int index) {
        double val = 0;
        try {
            val = (mops.get(index)).doubleValue();
        } catch (IndexOutOfBoundsException e) {
            throw new IndexOutOfBoundsException(
                    "\n  Measure of Performance Number'" + index
                            + "' was never computed for state: " + this);
        } catch (NullPointerException e) {
            throw new NullPointerException(
                    "\n  Error getting Measure of Performance '" + index
                            + "'  for state: " + this + " "
                            + ((mops == null) ? ". 'mops'" + " is null" : ""));
        }
        return val;
    }

  

    /**
     * This method is called when a state is added to a set, if assertions are
     * enabled. You should include code that checks the consistency of the
     * paprameters entered. It is very helpful during depelopment. Once
     * assertions are disabled, this will not reduce the speed of your program.
     * 
     * @return true if the state is consistent.
     */
    public abstract boolean isConsistent();

    /**
     * Returns a (hopefully short) label that descibes the State. It is used by
     * all print methods and in the GUI.
     * 
     * @return A shor String label that identifies the state.
     */
    public abstract String label();

    /**
     * Returns a String that describes the State. By default it is an empty
     * string, but you should implement it in order to get a meaningful
     * description.
     * 
     * @return A String description of the State
     */
    public String description() {
        return "";
    }

    /**
     * Returns the label.
     * 
     * @see #label()
     */
    @Override
    public final String toString() {
        return label();
    }

    /**
     * @return The index in the State set
     */
    public final int getIndex() {
        return index;
    }

    /**
     * Sets the index number in the state set. This method is called by getStates.
     * 
     * @param i
     *            number in states set.
     */
    final void setIndex(int i) {
        index = i;
    }

} // class
