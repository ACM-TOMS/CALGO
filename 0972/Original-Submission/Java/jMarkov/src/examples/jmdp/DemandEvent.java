package examples.jmdp;

import jmarkov.basic.PropertiesEvent;

/**
 * This class represent a demand event in an inventory system and it is used by
 * many examples.
 * 
 * @author Andrés Sarmiento, Germán Riaño. Universidad de los Andes. (C) 2006
 * 
 */
public class DemandEvent extends PropertiesEvent {

    private boolean greaterThan;

    /**
     * @param d
     * @param greater
     */
    public DemandEvent(int d, boolean greater) {
        super(new int[] { d });
        greaterThan = greater;
    }

    /**
     * Return the demand size
     * 
     * @return The level.
     */
    public int getDemand() {
        return getProperty(0);

    }

    /**
     * 
     * 
     * @return True when the demand is greater than demand and false if the
     *         demand is equal to demand.
     */
    public boolean getGreaterThan() {
        return greaterThan;
    }

    @Override
    public String label() {
        return "Demand " + getDemand();
    }
}
