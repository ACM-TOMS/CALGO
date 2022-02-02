/*
 * Created on 07/08/2005
 *
 */
package jmarkov.basic;

/**
 * This class is an easy way to use an event that is represented by an array of
 * int.
 * 
 * @author Andres Sarmiento, Germán Riaño - Universidad de Los Andes
 */
public class PropertiesEvent extends Event implements PropertiesElement {
    /**
     * The array that characterized the Event.
     */
    private final int[] prop;

    /**
     * Builds a new event with characteristic array as a paramenter.
     * 
     * @param status
     *            characteristic array of the event.
     */
    public PropertiesEvent(int[] status) {
        super();
        this.prop = status;
    }

    @Override
    public String label() {
        String localName = "(";
        for (int i = 0; i < prop.length - 1; i++) {
            localName += prop[i] + ",";
        }
        return localName + prop[prop.length - 1] + ")";
    }

    /**
     * Creates a new PropertiesEvent with an array of the size indicated filled
     * with zeros.
     * 
     * @param size
     *            size of the characteristic array.
     */
    public PropertiesEvent(int size) {
        this(new int[size]);
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

    /**
     * Compares the properties in order.
     * 
     * @param e1
     *            The PropertiesEvent to compare to.
     * @return +1, 0 or -1, according to whether this states dominates e1.
     */
    public final int compareTo(PropertiesEvent e1) {
        for (int i = 0; i < prop.length; i++) {
            if (prop[i] < e1.prop[i])
                return 1;
            if (prop[i] > e1.prop[i])
                return -1;
        }
        return 0;
    }

    /*
     * (non-Javadoc)
     * 
     * @see Event#compareTo(Event)
     */
    @Override
    public final int compareTo(Event e) {
        if (!(e instanceof PropertiesEvent))
            return -1;// maybe it should return an error?
        return this.compareTo((PropertiesEvent) e);
    }

    /**
     * @param e
     *            teh PropertiesEvent to compare To
     * @return true if Events are equal
     */
    public final boolean equals(PropertiesEvent e) {
        return (compareTo(e) == 0);
    }

    // @Override
    // public int hashCode() {
    // return status.hashCode();
    // }
    /**
     * Gets the array of properties.
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

    @Override
    public PropertiesEvent clone() {
        return new PropertiesEvent(prop);
    }

}
