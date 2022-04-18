package examples.jmdp;

import jmarkov.basic.Action;

/**
 * This class represents an order in an inventory system.
 * It is used in many of the Examples.
 * 
 * @author Germán Riano, Andres Sarmiento
 */
public class Order extends Action {
    private int size;

    /**
     * Defualt constructor. Recives the size order.
     * @param k
     */
    Order(int k) {
        size = k;
    }

    @Override
    public String label() {
        return "Order " + size + " Units";
    }

    public int compareTo(Action a) {
        if (a instanceof Order)
            return (size - ((Order) a).size);
        else 
            throw new IllegalArgumentException(
                    "Comparing with different type of Action.");
    }

    /**
     * @return Returns the order size.
     */
    public final int getSize() {
        return size;
    }

}
