package examples.jmarkov;

import java.io.PrintWriter;
import java.io.StringWriter;


//import java.net.URL;

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
 * This class represents a Narrow Aisle Picking System with 2 pickers
 * n picklocs and a 1:1 pick/walk speed ratio. Each picker picks at
 * a given location with probability p, or walks by with probability q=1-p
 * The pickers are blocked when one reaches the other and is unable to pass
 * because he is using the next location either picking or walking . 
 * We assume a circular warehouse without time spent at the I/O point.
 * 
 * @author Daniel F. Silva (2007)
 * 
 */
public class PickBlockNarrowSlow extends SimpleMarkovProcess<PickBlockNarrowSlowState, PickBlockNarrowSlowEvent> {
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
    public PickBlockNarrowSlow(double p, int picklocs, int[]a ) {
    	super(new PickBlockNarrowSlowState(a,picklocs), PickBlockNarrowSlowEvent.getAllEvents());
        this.p = p;
        this.picklocs = picklocs;
    }
    
    /**
     * Default Constructor used by GUI
     */
    public PickBlockNarrowSlow() {
    	this(0.1, 250, new int[] { 1, 0, 0 });
    }
    

    /**
     * Determine the active events.
     */

    @Override
    public boolean active(PickBlockNarrowSlowState i, PickBlockNarrowSlowEvent e) {
        boolean result = false;
        if (i.getFirst()==0 && i.getSecond()==0){
	        switch (e.getType()){
	        case WALKWALK:
	        	result = true; break;
	        case WALKPICK:
	        	result = true; break;
	        case PICKPICK:
	        	result = true; break;
	        case PICKWALK:
	        	result = true; break;
	        default:
	        	result = false; break;
	     	}
	    }
        else if (i.getFirst()==1 && i.getSecond()==0){
	        switch (e.getType()){
	        case WALKWALK:
	        	result = true; break;
	        case WALKPICK:
	        	result = true; break;
	        default:
	        	result = false; break;
	        }
        }
        else if (i.getFirst()==0 && i.getSecond()==1){
	        switch (e.getType()){
	        case WALKWALK:
	        	result = true; break;
	        case PICKWALK:
	        	result = true; break;
	        default:
	        	result = false; break;
	        }
        }
        else if (i.getFirst()==1 && i.getSecond()==1){
	        switch (e.getType()){
	        case WALKWALK:
	        	result = true; break;
	        default:
	        	result = false; break;
	        }
        }
        if (i.getDistance()== 1) {
        	if( i.getFirst()==0 && i.getSecond()==1){
		        switch (e.getType()){
		        case WALKWALK:
		        	result = true; break;
		        default:
		        	result = false; break;
		        }
        	}
        	if ( i.getFirst()==0 && i.getSecond()==0){
		        switch (e.getType()){
		        case WALKWALK:
		        	result = true; break;
		        case WALKPICK:
		        	result = true; break;
		        case PICKPICK:
		        	result = true; break;
		        case TOBLOCK1:
		        	result = true; break;
		        default:
		        	result = false; break;
		        }
        	} 
        }
        if (i.getDistance()== picklocs-1) {
        	if( i.getFirst()==1 && i.getSecond()==0){
		        switch (e.getType()){
		        case WALKWALK:
		        	result = true; break;
		        default:
		        	result = false; break;
		        }
        	}
        	if ( i.getFirst()==0 && i.getSecond()==0){
		        switch (e.getType()){
		        case WALKWALK:
		        	result = true; break;
		        case PICKWALK:
		        	result = true; break;
		        case PICKPICK:
		        	result = true; break;
		        case TOBLOCK2:
		        	result = true; break;
		        default:
		        	result = false; break;
		        }
        	} 
        }
        return result;
    }

    @Override
    public States<PickBlockNarrowSlowState> dests(PickBlockNarrowSlowState i, PickBlockNarrowSlowEvent e) {
        int [] origin = i.getProperties();
        return new StatesSet<PickBlockNarrowSlowState>(i.Trans(origin, e));
    }

    /**
     * Returns the transition probability from State i to State j.
     */

    @Override
    public double rate(PickBlockNarrowSlowState i, PickBlockNarrowSlowState j, PickBlockNarrowSlowEvent e) {
    	int [] origin = i.getProperties();
        int [] dest = j.getProperties();
        if (this.active(i,e)==false || Math.abs(origin[0]-dest[0])>1) 
    		return 0;
        return getProbs(origin,dest);
    }

    public double getProbs(int []i, int []j) {
        double pr=0.0, q=1-p;
        if (i[1]==1 && i[2]==1){
	        if (j[1]==0 && j[2]==0 && i[0]==j[0]) pr=1.0;
	        else pr= 0.0;
	    }
        else if (i[1]==1 && i[2]==0){
        	if (j[1]==0 && j[2]==0 && i[0]==j[0])pr= q;
        	else if (j[1]==0 && j[2]==1 && i[0]+1==j[0])pr= p;
	        else pr= 0.0;    
        }
        else if (i[1]==0 && i[2]==1){
        	if (j[1]==0 && j[2]==0 && i[0]==j[0])pr= q;
        	else if (j[1]==1 && j[2]==0 && i[0]-1==j[0])pr= p;
	        else pr= 0.0;    
        }
        else if (i[1]==0 && i[2]==0){
        	if (j[1]==0 && j[2]==0 && i[0]==j[0])pr= q*q;
        	else if (j[1]==0 && j[2]==1 && i[0]+1==j[0])pr= p*q;
        	else if (j[1]==1 && j[2]==0 && i[0]-1==j[0])pr= p*q;
        	else if (j[1]==1 && j[2]==1 && i[0]==j[0])pr= p*p;
	        else pr= 0.0;    
        }
        if (i[0]== 1 && i[1]==0 && i[2]==1) {
        	if( j[0]==1 && i[1]==0 && j[2]==0) pr=1.0;
        	else pr=0;
        	}
         if (i[0]== 1 && i[1]==0 && i[2]==0){
        	if (j[1]==0 && j[2]==0 && i[0]==j[0])pr= q*q;
        	else if (j[1]==0 && j[2]==1 && i[0]+1==j[0])pr= p*q;
        	else if (j[1]==0 && j[2]==1 && i[0]==j[0])pr= p*q;
        	else if (j[1]==1 && j[2]==1 && i[0]==j[0])pr= p*p;
	        else pr= 0.0;
        } 
        if (i[0]== picklocs-1 && i[1]==1 && i[2]==0) {
        	if( j[0]==picklocs-1 && j[1]==0 && j[2]==0) pr=1.0;
        	else pr=0.0;
        	}
        if (i[0]== picklocs-1 && i[1]==0 && i[2]==0){
        	if (j[1]==0 && j[2]==0 && i[0]==j[0])pr= q*q;
        	else if (j[1]==1 && j[2]==0 && i[0]==j[0])pr= p*q;
        	else if (j[1]==1 && j[2]==0 && i[0]-1==j[0])pr= p*q;
        	else if (j[1]==1 && j[2]==1 && i[0]==j[0])pr= p*p;
	        else pr= 0.0;
        } 
        
        return pr;
    }
    
    @Override
    public String description() {
        String stg = "NarrowSlow-Aisle-OPS \n\nThere are " + picklocs + " pick locations in the system\n";
        stg += "There are 2 Pickers in the system.\n";
        stg += "The pick/walk speed ratio is 1:1.\n";
        stg += "The picking Probability is " + p + ".\n";
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw, true);
        stg += sw.toString();
        return stg;
    }


    /**
     * Replaces the method to print MOPs. This method is called in
     * <code>toString()</code>, which is also used in the MOP's tab in the
     * GUI.
     */
    
    @Override
    public int printMOPs(PrintWriter out, int width, int decimals) {
        
        try {
         double [] busy = getSteadyState();
         out.printf("    Number of Pick Locations   = %6d ", this.picklocs);
         out.printf("\n    Number of Pickers          = %6d ", 2);
         out.printf("\n    Picking Probability        = %2.6f ", this.p);
         out.printf("\n    Percentage of Time Picker 1 is blocked    = %2.6f ", busy[1]*100);
         out.printf("\n    Percentage of Time Picker 2 is blocked    = %2.6f ", busy[4*picklocs-6]*100);
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

        PickBlockNarrowSlow OPS = new PickBlockNarrowSlow();
        OPS.setDebugLevel(0);

        OPS.showGUI();
        OPS.printAll();
        OPS.printMOPs();
    }

    
    /**
     * @return Returns the picking probability.
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

} // End of PickBlockNarrowSlow class

/**
 * A State in a Picking System. This class represents a State in a
 * Narrow-Aisle-OPS with 2 pickers in the system. Every picker picks
 * at each location with probability p. The State consist of a scalar
 * denoting the distance between pickers, and a scalar with value either
 * 0 or 1 for each picker depending on whether they just walked(0) or picked(1). 
 * 
 */

class PickBlockNarrowSlowState extends PropertiesState {

	int picklocs; 
    /**
     * Constructor for the PickBlockNarrowSlowState class
     * 
     * @param Property 0
     *            distance between pickers
     * @param Property 1 and 2
     *            Current state of each picker
     *            
     */
    public PickBlockNarrowSlowState(int [] position, int picklocs) {	
        super(position);
        this.picklocs=picklocs;
        setProperty(0, position [0]);
        setProperty(1, position [1]);
        setProperty(2, position [2]);
    }

    /**
     * @see State#computeMOPs(MarkovProcess)
     */

    /**
     * Returns the distance between 2 pickers.
     * 
     */
    public int getDistance() {
        return getProperty(0);
    }
    
    /**
     * Returns the First picker's last activity 1: picked, 0: walked.
     * 
     */
    public int getFirst() {
        return getProperty(1);
    }
    
    /**
	 *	Returns the Second picker's last activity 1: picked, 0: walked.
     * 
     */
    public int getSecond() {
        return getProperty(2);
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
     *            Original position
     * @param e
     *            Event
     * @return The new state
     */
    public PickBlockNarrowSlowState Trans(int [] i, PickBlockNarrowSlowEvent e) {
        int [] newpos;
        newpos = new int [3];
        newpos [0]= i[0] + e.getDistshift();
        newpos [1]= e.getFirst();
        newpos [2]= e.getSecond();
        return new PickBlockNarrowSlowState(newpos, picklocs);
    }

    /**
     * Describes the States
     * 
     * @return a String description of the State.
     */
    @Override
    public String description() {
        int a = getDistance();
        int b = getFirst();
        int c = getSecond();
        int locs=this.picklocs;
        String stg;
        if (a==1 && b==0 && c==1){stg="Picker 2 blocks picker 1.";}
        else if (a==locs-1 && b==1 && c==0){stg="Picker 1 blocks picker 2";}
        else {
        	stg="The distance between pickers 1 and 2 is ";stg+=a;
        	if (b==1)
        		stg+=", Picker 1 just picked";
        	else
        		stg+=", Picker 1 just walked";
        	if (c==1)
        		stg+=", Picker 2 just picked.";
        	else
        		stg+=", Picker 2 just walked.";
        	}
        return stg;
    }

} // End of PickBlockNarrowSlowState class.

/**
 * Each Event characterizes the shift in distance from i to j. and a posiible change
 * in each picker last activity.
 * 
 * @author Daniel Silva
 */
class PickBlockNarrowSlowEvent extends jmarkov.basic.Event {
    
	/** Event types. */
    public static enum Type {
        /** Both pickers just picked */
        PICKPICK,
    	/** Both pickers just walked */
        WALKWALK,
        /** Picker 1 picked, 2 walked. */
        PICKWALK,
        /** Picker 1 walked, 2 picked. */
        WALKPICK,
        /** Picker 1 just got blocked. */
        TOBLOCK1,
        /** Picker 2 just got blocked. */
        TOBLOCK2,;
    }

    private Type type; // event type
    public int distshift; // Shift in distance 
    public int first; // Activity of first picker
    public int second; // Activity of second picker

    /**
     * Creates an Event.
     * 
     * @param type
     */
    public PickBlockNarrowSlowEvent(Type type) {
    	this.type = type;
    	switch (type){
    	case WALKWALK: 
    		this.distshift = 0;
    		this.first=0;
    		this.second=0;
    		break;
    	case WALKPICK: 
    		this.distshift = 1;
    		this.first=0;
    		this.second=1;
    		break;
    	case PICKWALK: 
    		this.distshift = -1;
    		this.first=1;
    		this.second=0;
    		break;
    	case PICKPICK: 
    		this.distshift = 0;
    		this.first=1;
    		this.second=1;
    		break;
    	case TOBLOCK1: 
    		this.distshift = 0;
    		this.first=0;
    		this.second=1;
    		break;
    	case TOBLOCK2: 
    		this.distshift = 0;
    		this.first=1;
    		this.second=0;
    		break;	
    	}
    }


    /**
     * @return event type
     */
    public Type getType() {
        return type;
    }
    
    /**
     * @return distance change
     */
    
    public int getDistshift() {
		return distshift;
	}

    /**
     * @return first picker's last activity
     */
    
	public int getFirst() {
		return first;
	}

    /**
     * @return second picker's last activity
     */
	
	public int getSecond() {
		return second;
	}

	/**
     * 
     * @return A set with all the events in the system.
     */
    public static EventsSet<PickBlockNarrowSlowEvent> getAllEvents() {
        EventsSet<PickBlockNarrowSlowEvent> eSet = new EventsSet<PickBlockNarrowSlowEvent>();
        eSet.add(new PickBlockNarrowSlowEvent(Type.PICKPICK));
        eSet.add(new PickBlockNarrowSlowEvent(Type.PICKWALK));
        eSet.add(new PickBlockNarrowSlowEvent(Type.WALKPICK));
        eSet.add(new PickBlockNarrowSlowEvent(Type.WALKWALK));
        eSet.add(new PickBlockNarrowSlowEvent(Type.TOBLOCK1));
        eSet.add(new PickBlockNarrowSlowEvent(Type.TOBLOCK2));
        return eSet;
    }
    /**
     * 
     * String Description of this event
     */
    
    @Override
    public String label() {
        String stg = "";
        switch (type) {
        case PICKPICK:
        	stg = "Both Pick";
            break;
        case WALKWALK:
            stg = "Both walk";
            break;
        case PICKWALK:
            stg = "The First Picker picks, the second walks";
            break;
        case WALKPICK:
            stg = "The first picker walks the second picks";
            break;
        case TOBLOCK1:
            stg = "Picker 1 gets blocked by picker 2";
            break;
        case TOBLOCK2:
            stg = "Picker 2 gets blocked by picker 1";
            break;
        }
        return stg;
    }

}