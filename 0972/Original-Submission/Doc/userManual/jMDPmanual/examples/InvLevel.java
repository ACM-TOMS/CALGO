/*
 * Created on 26/06/2004
 *
 */
package examples.jmdp;

import jmarkov.MarkovProcess;
import jmarkov.basic.PropertiesState;

/**
 * This class allows to represent a State with a songle integer.
 * It's used in many of the examples.
 * @author Daniel Silva, German Riano, Andres Sarmiento. Universida de los Andes
 */
public class InvLevel extends PropertiesState {
    /**
     * Default constructor.
     * @param k The level
     */
	public InvLevel(int k) {
		super(new int[] {k});
	}

	/**
	 * Return the inventory level
	 * 
	 * @return The level.
	 */
	public int getLevel() {
		return prop[0];
	}

	@Override
	public String label() {
		return "Level " + getLevel();
	}

    /**
     * @see jmarkov.basic.State#computeMOPs(MarkovProcess)
     */
    @Override
    public void computeMOPs(MarkovProcess mp) {
        // TODO Auto-generated method stub
    }


    /**
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        // TODO Auto-generated method stub
        return true;
    }

}
