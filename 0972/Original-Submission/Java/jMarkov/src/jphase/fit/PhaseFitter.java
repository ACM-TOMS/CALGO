package jphase.fit;

import jphase.PhaseVar;

/**
 * This class defines the methods and attributes for any class implementing 
 * fitting methods for Phase-Type distributions 
 * @author Juan F. Pérez
 * @version 1.0  
 */
public interface PhaseFitter {
   

    /**
     * Executes the fitting procedure to find the parameter set
     * @return Phase variable found
     */
    public PhaseVar fit();
    

}
