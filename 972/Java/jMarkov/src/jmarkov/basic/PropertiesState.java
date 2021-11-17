package jmarkov.basic;

import jmarkov.MarkovProcess;

/**
 * The states are characterized by an array of integer-valued properties, whose
 * meaning will change from implementation to implementation. The class must be
 * extended.
 * 
 * @author German Riano, Andres Sarmiento. Universidad de los Andes.
 * @version 1.0a
 */

public class PropertiesState extends State implements PropertiesElement {

    /**
     * This array contains the properties that characterize the state.
     */
    protected final int[] prop; // Properties array

    /**
     * Constructs a State charcterized by K properties. The original values of
     * all the properties is 0.
     * 
     * @param K
     *            Number of Properties
     */
    public PropertiesState(int K) {
        prop = new int[K];
    }

    /**
     * This creates a PropertiesState with the given array. WARNNING: the array
     * is NOT internally copied, so it is assumed then NO changes are made to
     * this array after it is given to the constructor.
     * 
     * @param properties
     *            An integer valued array with the properties that characterize
     *            this state.
     */
    public PropertiesState(int[] properties) {
        this(properties, true);
    }

    /**
     * This creates a PropertiesState with the given array.
     * 
     * @param properties
     *            An integer valued array with the properties that characterize
     *            this state.
     * @param deepCopy
     *            true if the constructor should make a deep copy of the array.
     *            This causes some overhead but increased security.
     */
    public PropertiesState(int[] properties, boolean deepCopy) {
        if (deepCopy) {
            int N = properties.length;
            prop = new int[N];
            System.arraycopy(properties, 0, prop, 0, N);
        } else {
            this.prop = properties;
        }
    }


    /**
     * Constructs a new State by cloning the given State. If you are extending
     * PropertiesState you may want to include a code like:
     * 
     * <pre>
     * public YourState clone() {
     *     return new PropertiesState(this);
     * }
     * </pre>
     * 
     * @param s
     *            A given State.
     */
    public PropertiesState(PropertiesState s) {
        int n = s.prop.length;
        prop = new int[n];
        System.arraycopy(s.prop, 0, prop, 0, n);
    }

    /**
     * By default it computes the long run average for each property. the user
     * should override this method in order to compute more meaningful measures
     * of performance.
     * 
     * @see jmarkov.basic.State#computeMOPs(MarkovProcess)
     */
    @Override
    public void computeMOPs(MarkovProcess mp) {
        int n = prop.length;
        for (int i = 0; i < n; i++) {
            setMOP(mp, "Average for Property " + i, prop[i]);
        }
    }

    /*
     * (non-Javadoc)
     * 
     * @see State#compareTo(State)
     */
    @Override
    public final int compareTo(State s) {
        if (!(s instanceof PropertiesState))
            return -1;// maybe it should return an error?
        return this.compareTo((PropertiesState) s);
    }

    /**
     * Compares according to the internal properties in lexicographic order.
     * 
     * @param ps
     *            The state to compare to
     * @return +1,-1,0 depending on relative ordering.
     * @see jmarkov.basic.State#compareTo(State)
     */
    public final int compareTo(PropertiesState ps) {
        int K = prop.length;
        int K2 = ps.prop.length;
        if (K != K2)
            return (K < K2) ? -1 : 1;// not equal if they are not the same
        // size
        for (int i = 0; i < K; i++) {
            // As soon as a discrepancy is found we exit.
            if (prop[i] < ps.prop[i])
                return -1;
            if (prop[i] > ps.prop[i])
                return +1;
        }
        return 0;
    }

    /**
     * Returns a string representation of this state in vector form. The String
     * will be in the form (p1,p2,..,pK). A Class implementing this Class could
     * give a more meaningful description.
     */
    @Override
    public String label() {
        String stg = "(";
        int n = prop.length;
        for (int i = 0; i < (n - 1); i++) {
            stg += prop[i] + ",";
        }
        stg += prop[n - 1] + ")";
        return stg;
    }

    /**
     * @see jmarkov.basic.PropertiesElement#getNumProps()
     */
    public int getNumProps() {
        return prop.length;
    }

    /**
     * Crates a copy of the properties array.
     * @see jmarkov.basic.PropertiesElement#getProperties()
     */
    public final int[] getProperties() {
        int newProp[] = new int[prop.length];
        System.arraycopy(prop, 0, newProp, 0, prop.length);
        return newProp;
    }

    /**
     * @see jmarkov.basic.PropertiesElement#getProperty(int)
     */
    public int getProperty(int index) {
        return prop[index];
    }

    /**
     * Sets the value of the property at the given index
     * 
     * @param index
     * @param value
     */
    protected void setProperty(int index, int value) {
        prop[index] = value;
    }

    /**
     * @see jmarkov.basic.PropertiesElement#clone()
     */
    @Override
    public PropertiesElement clone() {
        return new PropertiesState(prop);
    }

    /**
     * It is strongly recommended that the user implements this method. If left
     * unimplemented this method returns true.
     * 
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        return true;
    }

}
