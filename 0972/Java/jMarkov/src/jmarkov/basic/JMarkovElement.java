/**
 * JMarkovElement.java
 * Created: Feb 23, 2006
 */
package jmarkov.basic;

/**
 * All the elements in JMarkov implement this interface, so they can be easily
 * described in the interface. It is recommended that the method
 * <code>toString()</code> is implemented as final, and calling
 * <code>label()</code>.
 * 
 * 
 * @author Germán Riaño. Universidad de los Andes. (C) 2006
 * 
 */
public interface JMarkovElement {

    /**
     * This method returns a short String used in the user interface to describe
     * this element. It is highly recommended that every class calls label(),
     * using the following code:
     * 
     * <pre>
     * public final String toString() {
     *     return label();
     * }
     * </pre>
     * 
     * @return A String label.
     * @see #label()
     */
    public abstract String toString();

    /**
     * This method returns a short String used in the user interface to describe
     * this element.
     * 
     * @return A String label.
     * @see #description()
     */
    public String label();

    /**
     * This method return a complete verbal describtion of this element. This
     * description may contain multiple text rows.
     * 
     * @return A String describing this element.
     * @see #label()
     */
    public String description();

    /**
     * Returns true if these two elements are equal. If this element implementa a
     *         compareTo() method it is recommended that this method returns
     *         compareTo(o)==0.
     * @param e The Object to compare to.
     * 
     * @return True if the elements are equal. 
     * @see java.lang.Object#equals(java.lang.Object)
     */
    public abstract boolean equals(Object e);

}
