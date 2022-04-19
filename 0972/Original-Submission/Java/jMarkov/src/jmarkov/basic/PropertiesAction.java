/*
 * Created on 04/07/2005
 *
 */
package jmarkov.basic;

/**
 * This class is an easy way to use a Action that is represented by an integer
 * valued array.
 * 
 * @author Andres Sarmiento, Germán Riaño - Universidad de Los Andes
 * 
 */
public class PropertiesAction extends Action implements PropertiesElement {

    /** The array of properties */
    private final int[] prop;

    /**
     * Builds an object with the given array.
     * 
     * @param properties
     */
    public PropertiesAction(int[] properties) {
        super();
        this.prop = properties;
    }

    @Override
    public String label() {
        String name = "(";
        for (int i = 0; i < prop.length - 1; i++) {
            name += prop[i] + ",";
        }
        name = name + prop[prop.length - 1] + ")";
        return name;
    }

    /**
     * Creates an Action Object wit an array of the given size.
     * 
     * @param size
     */
    public PropertiesAction(int size) {
        super();
        this.prop = new int[size];
    }

    /**
     * @param a
     *            The action array to compare to
     * @return +1, -1 or 0.
     * @see java.lang.Comparable#compareTo
     */
    public final int compareTo(PropertiesAction a) {
        if (prop.length > a.prop.length)
            return +1;
        for (int i = 0; i < prop.length; i++) {
            if (prop[i] < a.prop[i])
                return 1;
            if (prop[i] > a.prop[i])
                return -1;
        }
        return 0;
    }

    /*
     * (non-Javadoc)
     * 
     * @see jmdp.Action#compareTo(jmdp.Action)
     */
    public final int compareTo(Action a) {
        if (!(a instanceof PropertiesAction))
            return -1;// maybe shoud return an error?
        return this.compareTo((PropertiesAction) a);
    }

    // @Override
    // public int hashCode() {
    // return properties.hashCode();
    // }
    /**
     * Gets thae array of properties.
     * 
     * @return Returns the properties array.
     */
    public final int[] getProperties() {
        int newProp[] = new int[prop.length];
        System.arraycopy(prop, 0, newProp, 0, prop.length);
        return newProp;
    }

    /**
     * Gets the value of this property.
     * 
     * @param index
     * @return the property at the given index
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
     * Returns the number of properties in the array that characterize this
     * element.
     * 
     * @return The number of properties.
     */
    public int getNumProps() {
        return prop.length;
    }

    @Override
    public PropertiesAction clone() {
        return new PropertiesAction(prop);
    }

}
