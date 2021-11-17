package jphase.GUI;

import jphase.PhaseVar;

/**
 * This class keeps a Phase-type variable together with its name.
 * @author Juan F. Perez
 *
 */
public class PhaseVarInfo {
    /**
     * Name of the PH variable
     */
    public String varName;
    
    /**
     * PH variable
     */
    public PhaseVar var;
    
    /**
     * 
     * @param varName
     * @param var
     */
    public PhaseVarInfo(String varName, PhaseVar var) {
        this.varName = varName;
        this.var = var.copy();        
    }
    
    @Override
    public String toString() {
        return varName;
    }
}
