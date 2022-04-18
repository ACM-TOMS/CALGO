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
 * This class represents a Narrow Aisle Picking System with 2 pickers,
 * n picklocs and a infinity:1 pick/walk speed ratio. Each picker picks at
 * a given location with probability p, or walks by with probability q=1-p
 * The pickers are blocked when one reaches the other and is unable to pass
 * because he is picking. We assume a circular warehouse without time 
 * spent at the I/O point.
 * 
 * @author Daniel F. Silva (2007)
 * 
 */
public class PickBlockNarrow extends SimpleMarkovProcess<PickBlockNarrowState, PickBlockNarrowEvent> {
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
    
    public PickBlockNarrow(double p, int picklocs) {
    	super(new PickBlockNarrowState(0,picklocs), PickBlockNarrowEvent.getAllEvents(picklocs));
        this.p = p;
        this.picklocs = picklocs;
    }
    
    /**
     * Default Constructor used by GUI
     */
    
    public PickBlockNarrow() {
    	this(0.5, 100);
    }
    
    /**
     * Determine the active events.
     */

    @Override
    public boolean active(PickBlockNarrowState i, PickBlockNarrowEvent e) {
        boolean result = false;
        int stat = i.getDistance();
        int dest = stat + e.getShift();
        result = (dest <= picklocs && dest >= 0);
        return result;
    }

    @Override
    public States<PickBlockNarrowState> dests(PickBlockNarrowState i, PickBlockNarrowEvent e) {
        int origin = i.getDistance();
        int dest = e.getShift();
        return new StatesSet<PickBlockNarrowState>(i.Trans(origin, dest));
    }

    /**
     * Returns the transition probability from State i to State j.
     */

    @Override
    public double rate(PickBlockNarrowState i, PickBlockNarrowState j, PickBlockNarrowEvent e) {
    	int origin = i.getDistance();
        int shift = e.getShift();
        int dest = origin+shift;
        return getProbs(origin,dest);
    }

    public double getProbs(int i, int j) {
    	double q, prob=0;
    	q=1.0-p;
    	if (i==0)i=1;
    	else if (i==picklocs)i=picklocs-1;
    	prob = (Math.pow(q,Math.abs(i-j))) /(1+q);  
    	if (j!=0 && j!=picklocs)prob*=p;
        return prob;
    }
    
    @Override
    public String description() {
        String stg = "Narrow-Aisle-OPS \n\nThere are " + picklocs + " pick locations in the system\n";
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
         out.printf("\n    Percentage of Time Picker 1 is Blocked   = %3.6f ", busy[0]*100);
         out.printf("\n    Percentage of Time Picker 2 is Blocked    = %3.6f ", busy[this.picklocs]*100);
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
        PickBlockNarrow OPS = new PickBlockNarrow();
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

} // End of PickBlockNarrow class

/**
 * A State in a Picking System. This class represents a State in a
 * Narrow-Aisle-OPS with 2 pickers in the system. Every picker picks
 * at each location with probability p. The State consist of a scalar
 * denoting the distance between pickers, note that 0 and n mean blocking
 * 
 */

class PickBlockNarrowState extends PropertiesState {

	int picklocs; 
    /**
     * Constructor for the PickBlockNarrowState class
     * 
     * @param distance
     *            distance between pickers
     */
    public PickBlockNarrowState(int distance, int picklocs) {	
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
    public PickBlockNarrowState Trans(int i, int j) {
        int newdist=i+j;
        return new PickBlockNarrowState(newdist, picklocs);
    }

    /**
     * Describes the States
     * 
     * @return a String description of the State.
     */
    @Override
    public String description() {
        int M = getDistance();
        int locs=this.picklocs;
        String stg;
        if (M==0){stg="Picker 2 blocks picker 1.";}
        else if (M==locs){stg="Picker 1 blocks picker 2.";}
        else {stg="The distance between picker 1 and 2 is ";stg+=M;}
        return stg;
    }

} // End of PickBlockNarrowState class.

/**
 * Each Event characterizes the shift in distance from i to j.
 * 
 * @author Daniel Silva
 */
class PickBlockNarrowEvent extends Event {

    private int shift;

    /**
     * Constructor of PickBlockNarrowEvent
     * 
     * @param origin
     *            Origin state index
     * @param dest
     *            Destination State index
     */
    public PickBlockNarrowEvent(int shift) {
        this.shift = shift;
    }

    /**
     * Returns the distance shift of this event
     * 
     * @return Distance shift of this event
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
    public static EventsSet<PickBlockNarrowEvent> getAllEvents(int picklocs) {
        EventsSet<PickBlockNarrowEvent> theSet = new EventsSet<PickBlockNarrowEvent>();
        for (int i = -picklocs; i <= picklocs; i++) {
            		theSet.add(new PickBlockNarrowEvent(i));
            }
        return theSet;
    }

    /**
     * @see java.lang.Object#toString()
     */
    @Override
    public String label() {
    	if (shift>0)
    		return "The distance is increased by " + (shift) + " positions.";
    	else if(shift<0)
    		return "The distance is decreased by " + (-shift) + " positions.";
    	else
    		return "The distance is unchanged.";
    }

}