package examples.jmdp;

import jmarkov.basic.PropertiesState;

/**
* This class represents a state in a tandem queuing system it is used by the Access Control examples
   * @author Daniel Silva 
     * 
     */

public class TandemQueues extends PropertiesState {
    /**
     * @param custs
     *            State of the system
     */
    public TandemQueues(int[] custs) {
        super(custs);
    }

    /**
     * @return Number of customers in queue 1
     */
    public int getQ1() {
        return prop[0];
    }// in system

    /**
     * @return Number of customers in queue 2
     */
    public int getQ2() {
        return prop[1];
    }// in system

}