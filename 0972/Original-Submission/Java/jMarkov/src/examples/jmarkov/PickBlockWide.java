package examples.jmarkov;

import java.io.PrintWriter;
import java.io.StringWriter;
import jmarkov.MarkovProcess;
import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.NotUnichainException;
import Jama.Matrix;

/**
 * This class represents a Wide Aisle Picking System with 2 pickers
 * n picklocs and a infinity:1 pick/walk speed ratio. Each picker picks at
 * a given location with probability p, or walks by with probability q=1-p
 * The pickers are blocked when one reaches the other and is unable to pick at
 * a given location because the other is picking at the same location . 
 * We assume a circular warehouse without time spent at the I/O point.
 * 
 * @author Daniel F. Silva (2007)
 * 
 */

public class PickBlockWide extends SimpleMarkovProcess<PickBlockWideState, PickBlockWideEvent> {
    private double p;
    private int picklocs;
    
 
   /**
     * General constructor.
     * 
     * @param p
     *            Pick Probability.
     * @param picklocs
     *            Number of Pick Locations
     * 
     */
    public PickBlockWide(double p, int picklocs) {
    	super(new PickBlockWideState(0,picklocs), PickBlockWideEvent.getAllEvents(picklocs));
        this.p = p;
        this.picklocs = picklocs;
    }
    
    /**
     * Default Constructor used by GUI
     */
    public PickBlockWide() {
    	this(0.1, 20);
    }
    

    /**
     * Determine the active events.
     */

    @Override
    public boolean active(PickBlockWideState i, PickBlockWideEvent e) {
        boolean result = false;
        int stat = i.getDistance();
        int dest = stat + e.getShift();
        result = (dest < picklocs && dest >= 0);
        return result;
    }

    @Override
    public States<PickBlockWideState> dests(PickBlockWideState i, PickBlockWideEvent e) {
    	int origin = i.getDistance();
        int dest = e.getShift();
        return new StatesSet<PickBlockWideState>(i.Trans(origin, dest));
    }

    /**
     * Returns the transition probability from State i to State j.
     */

    @Override
    public double rate(PickBlockWideState i, PickBlockWideState j, PickBlockWideEvent e) {
    	int origin = i.getDistance();
        int shift = e.getShift();
        int dest = origin+shift;
        return getProbs(origin,dest);
    }

    public double getProbs(int i, int j) {
    	int r=0, loc=getPicklocs();
    	double q=0, denom=0, prob=0;
    	double []h;
    	h= new double [loc];
    	q=1.0-p;
    	if (i==0){
    		for (int k=1; k<=loc-1; k++){
    			r=k+i;
    			h[r]=(p*Math.pow(q, (Math.abs(k)-1))+p*Math.pow(q, (loc-Math.abs(k)-1)))/2;
        		denom+=h[r];
    		}
    		h[0]=0;
    		denom+=h[0];
    		prob=h[j]/denom;
    	}
    	else {
    		for (int k=-i; k<=loc-1-i; k++){
    			r=k+i;
    			if (r==i){
    				continue;    				
    			}
    			else{
        		h[r]=(p*Math.pow(q, Math.abs(k))+p*Math.pow(q, (loc-Math.abs(k))))/(1+q);
        		denom+=h[r];
    			}
    		}
    		h[i]=p/(1+q);
    		denom+=h[i];
    		prob=h[j]/denom;
    	}
        return prob;
    }
    
    @Override
    public String description() {
        String stg = "Wide-Aisle-OPS \n\nThere are " + picklocs + " pick locations in the system\n\n";
        stg += "There are 2 Pickers in the system.\n";
        stg += "The pick/walk speed ratio is 1:infinity.\n";
        stg += "Picking Probability is " + p + ".";
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw, true);
        stg += sw.toString();
        return stg;
    }


    /**
     * Replaces the method to printo MOPs. This method is called in
     * <code>toString()</code>, which is also used in the MOP's tab in the
     * GUI.
     */
    
    @Override
    public int printMOPs(PrintWriter out, int width, int decimals) {
        
        try {
         double [] busy = getSteadyState();
         out.printf("    Number of Pick Locations   = %6d ", this.picklocs);
         out.printf("\n    Number of Pickers          = %6d ", 2);
         out.printf("\n    Picking Probability        = %6.4f ", this.p);
         out.printf("\n    Percentage of Time Picker 1 is blocked    = %3.5f ", busy[0]*100/2);
         out.printf("\n    Percentage of Time Picker 2 is blocked    = %3.5f ", busy[0]*100/2);
         out.println();
        }    
         catch (NotUnichainException e) {
            out.println(e);
        }

        return 0;
    } // Main


    
    /**
     * Main method.
     * 
     * @param a
     *            Not used.
     */
    public static void main(String[]a) {
        // as in handout:
        PickBlockWide OPS = new PickBlockWide();
        OPS.setDebugLevel(0);

        OPS.showGUI();
        OPS.printAll();
        OPS.printMOPs();
    }

    
    /**
     * @return Returns the probability.
     */
    public double getProbability() {
        return p;
    }

    /**
     * @return Returns the number of pick locations.
     */
    public int getPicklocs() {
        return picklocs;
    }

} // End of PickBlockWide class

/**
 * A State in a Picking System. This class represents a State in a
 * Wide-Aisle-OPS with 2 pickers in the system. Every picker picks
 * at each location with probability p. The State consist of a scalar
 * denoting the distance between pickers, note that 0 and n mean blocking
 * 
 */

class PickBlockWideState extends PropertiesState {

	int picklocs; 
    /**
     * Constructor for the PickBlockWideState class
     * 
     * @param distance
     *            distance between pickers
     */
    public PickBlockWideState(int distance, int picklocs) {	
        super(1);
        this.picklocs=picklocs;
        setProperty(0, distance);
    }

    /**
     * Returns the distance between the 2 pickers.
     * 
     */

    public int getDistance() {
        return getProperty(0);
    }

    /**
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        // TODO Complete
        return true;
    }


    /**
     * Creates a new state obtained by moving Pickers along the aisle.
     * 
     * @param i
     *            Original distance
     * @param j
     *            New distances
     * @return The new state
     */
    public PickBlockWideState Trans(int i, int j) {
        int newdist=i+j;
        return new PickBlockWideState(newdist, picklocs);
    }

    /**
     * Describes the States
     * 
     * @return a String description of the State.
     */
    @Override
    public String description() {
        int M = getDistance();
        String stg;
        if (M==0){stg="One of the pickers is blocked";}
        else {stg="The distance between pickers is ";stg+=M;}
        return stg;
    }

} // End of PickBlockWideState class.

/**
 * Each Event characterizes the shift in distance from i to j.
 * 
 * @author Daniel Silva
 */
class PickBlockWideEvent extends Event {

    private int shift;

    /**
     * Constructor of PickBlockWideEvent
     * 
     * @param origin
     *            Origin state index
     * @param dest
     *            Destination State index
     */
    public PickBlockWideEvent(int shift) {
        this.shift=shift;
    }

    /**
     * Returns the position shift of this event
     * 
     * @return Position shift index of this event
     */
    public int getShift() {
        return shift;
    }
    /**
     * Returns all the possible events in a given state
     * 
     * @param p
     *            Pick probability
     * @param picklocs
     *            Number of Pick Locations
	 *
     * @return All the possible events in a given state
     */
    public static EventsSet<PickBlockWideEvent> getAllEvents(int picklocs) {
        EventsSet<PickBlockWideEvent> theSet = new EventsSet<PickBlockWideEvent>();
        for (int i = -picklocs+1; i < picklocs; i++) {
           		theSet.add(new PickBlockWideEvent(i));
        }
        return theSet;
    }

    /**
     * @see java.lang.Object#toString()
     */
    @Override
    public String label() {
        return "Distance shifts " + (shift);
    }

}