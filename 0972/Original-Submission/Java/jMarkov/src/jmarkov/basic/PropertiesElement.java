/*
 * Created on 11/03/2006
 */
package jmarkov.basic;

/**
 * This interface is a wrapper for elements (States, Actions Events) that can be
 * represented by an arry of integers. Known implentations include
 * PropertiesState, Propertiesevent and PropertiesAction. Basic methods are
 * provided to access the internal array.
 * 
 * @author German Riano. Universidad de los Andes. (C) 2006
 * @see jmarkov.basic.PropertiesState
 * @see jmarkov.basic.PropertiesAction
 * @see jmarkov.basic.PropertiesEvent
 */
public interface PropertiesElement extends JMarkovElement{

    /**
     * Returns the number of properties in the array that characterize this
     * element.
     * 
     * @return The number of properties.
     */
    public abstract int getNumProps();

    /**
     * Gets the array of properties.
     * 
     * @return Returns the properties array.
     */
    public abstract int[] getProperties();

    /**
     * Gets the value of this property at this index.
     * 
     * @param index
     * @return the property at the given index
     */
    public abstract int getProperty(int index);


    public abstract PropertiesElement clone();

}