package examples.jmdp;

import jmarkov.basic.Action;

/**
 * This class represents an admission or rejection in an admission control model.
 * It is used in the admission control examples
 * 
 * @author Daniel Silva
 */

public class Admit extends Action {

	private int adm;

	public Admit(int a) {
       adm=a;
    }

    @Override
    public String label() {
        return (adm==0)? "Reject" : "Admit";
    }

    public int compareTo(Action a) {
        if (a instanceof Action)
            return (adm - ((Admit) a).adm);
        else throw new IllegalArgumentException("Comparing to different type of Action.");
    }

    /**
     * @return Admit decision
     */
    public final int getAdm() {
        return adm;
    }
}

