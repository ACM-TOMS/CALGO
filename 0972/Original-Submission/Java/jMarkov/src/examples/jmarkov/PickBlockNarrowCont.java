package examples.jmarkov;

import java.io.PrintWriter;
import java.io.StringWriter;
import examples.jmarkov.PickBlockNarrowContEvent.Type;
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
 * n picklocs and exponential picking and walking times, each time can be different 
 * for each picker. Each picker picks at a given location with probability p, 
 * or walks by with probability q=1-p.
 * The pickers are blocked when one reaches the other and is unable to pass
 * because he is picking or walking in the next location. We assume a circular 
 * warehouse without time spent at the I/O point.
 * 
 * @author Daniel F. Silva (2007)
 * 
 */

public class PickBlockNarrowCont extends SimpleMarkovProcess<PickBlockNarrowContState, PickBlockNarrowContEvent> {
    private double p;
    private int picklocs;
    private double pickspeed1, walkspeed1, pickspeed2, walkspeed2;
    
 
   /**
     * General constructor.
     * 
     * @param p
     *            Pick Probability.
     * @param picklocs
     *            Number of Pick Locations
     * @param pickspeed1
     *            Rate of picking picks per minute for picker 1.
     * @param walkspeed1
     *            Rate of walking past locations locations per minute for picker 1.
     * @param pickspeed2
     *            Rate of picking picks per minute  for picker 2.
     * @param walkspeed2
     *            Rate of walking past locations locations per minute for picker 2.
     * 
     */
    
    public PickBlockNarrowCont(double p, int picklocs, double pickspeed1, double walkspeed1, double pickspeed2, double walkspeed2,int[]a ) {
    	super(new PickBlockNarrowContState(a,picklocs), PickBlockNarrowContEvent.getAllEvents());
        this.p = p;
        this.picklocs = picklocs;
        this.pickspeed1 = pickspeed1;
        this.walkspeed1 = walkspeed1;
        this.pickspeed2 = pickspeed2;
        this.walkspeed2 = walkspeed2;
    }
    
    /**
     * Default Constructor used by GUI
     */
    public PickBlockNarrowCont() {
    	this(0.1, 6, 1, 5, 1, 5, new int[] { 1, 0, 0 });
    }
    

    /**
     * Determine the active events.
     */

    @Override
    public boolean active(PickBlockNarrowContState i, PickBlockNarrowContEvent e) {
        boolean result = false;
        if (i.getDistance()==0 && i.getFirst()==1){
	        switch (e.getType()){
	        case FIRSTPICKTOPICK:
	        	result = true; break;
	        case FIRSTPICKTOWALK:
	        	result = true; break;
	        default:
	        	result = false; break;
	     	}
	    }
        else if (i.getDistance()==0 && i.getFirst()==0){
	        switch (e.getType()){
	        case FIRSTWALKTOPICK:
	        	result = true; break;
	        case FIRSTWALKTOWALK:
	        	result = true; break;
	        default:
	        	result = false; break;
	     	}
	    }
        else if (i.getDistance()==picklocs && i.getSecond()==0){
	        switch (e.getType()){
	        case SECWALKTOPICK:
	        	result = true; break;
	        case SECWALKTOWALK:
	        	result = true; break;
	        default:
	        	result = false; break;
	     	}
	    }
        else if (i.getDistance()==picklocs && i.getSecond()==1){
	        switch (e.getType()){
	        case SECPICKTOPICK:
	        	result = true; break;
	        case SECPICKTOWALK:
	        	result = true; break;
	        default:
	        	result = false; break;
	     	}
	    }
        else if (i.getFirst()==0 && i.getSecond()==0){
	        switch (e.getType()){
	        case FIRSTWALKTOPICK:
	        	result = true; break;
	        case FIRSTWALKTOWALK:
	        	result = true; break;
	        case SECWALKTOPICK:
	        	result = true; break;
	        case SECWALKTOWALK:
	        	result = true; break;
	        default:
	        	result = false; break;
	     	}
	    }
        else if (i.getFirst()==0 && i.getSecond()==1){
	        switch (e.getType()){
	        case FIRSTWALKTOPICK:
	        	result = true; break;
	        case FIRSTWALKTOWALK:
	        	result = true; break;
	        case SECPICKTOPICK:
	        	result = true; break;
	        case SECPICKTOWALK:
	        	result = true; break;
	        default:
	        	result = false; break;
	     	}
	    }
        else if (i.getFirst()==1 && i.getSecond()==0){
	        switch (e.getType()){
	        case FIRSTPICKTOPICK:
	        	result = true; break;
	        case FIRSTPICKTOWALK:
	        	result = true; break;
	        case SECWALKTOPICK:
	        	result = true; break;
	        case SECWALKTOWALK:
	        	result = true; break;
	        default:
	        	result = false; break;
	     	}
	    }
        else if (i.getFirst()==1 && i.getSecond()==1){
	        switch (e.getType()){
	        case FIRSTPICKTOPICK:
	        	result = true; break;
	        case FIRSTPICKTOWALK:
	        	result = true; break;
	        case SECPICKTOPICK:
	        	result = true; break;
	        case SECPICKTOWALK:
	        	result = true; break;
	        default:
	        	result = false; break;
	     	}
	    } 
        return result;
    }

    @Override
    public States<PickBlockNarrowContState> dests(PickBlockNarrowContState i, PickBlockNarrowContEvent e) {
        int [] origin = i.getProperties();
        return new StatesSet<PickBlockNarrowContState>(i.Trans(origin, e));
    }

    /**
     * Returns the transition probability from State i to State j, when event e occurs.
     */

    @Override
    public double rate(PickBlockNarrowContState i, PickBlockNarrowContState j, PickBlockNarrowContEvent e) {
    	int [] origin = i.getProperties();
        int [] dest = j.getProperties();
        Type a = e.getType();
        if (this.active(i,e)==false) 
    		return 0;
        return getProbs(origin, dest, a);
    }

    public double getProbs(int []i, int []j, Type e) {
        double pr=0.0, q=1-p;
        if (i[1]==0 && i[2]==0){
	        switch (e){
	        	case FIRSTWALKTOPICK:
	        		pr = p*walkspeed1;
	        		break;
	        	case FIRSTWALKTOWALK:
	        		pr = q*walkspeed1;
	        		break;
	        	case SECWALKTOPICK:
	        		pr = p*walkspeed2;
	        		break;
	        	case SECWALKTOWALK:
	        		pr = q*walkspeed2;
	        		break;
	        	default:
	        		break;
	        }
        }
        else if (i[1]==1 && i[2]==0){
        	switch (e){
	        	case FIRSTPICKTOPICK:
	        		pr = p*pickspeed1;
	        		break;
	        	case FIRSTPICKTOWALK:
	        		pr = q*pickspeed1;
	        		break;
	        	case SECWALKTOPICK:
	        		pr = p*walkspeed2;
	        		break;
	        	case SECWALKTOWALK:
	        		pr = q*walkspeed2;
	        		break;
	        	default:
	        		break;
        	}
        }
        else if (i[1]==0 && i[2]==1){
        	switch (e){
	        	case FIRSTWALKTOPICK:
	        		pr = p*walkspeed1;
	        		break;
	        	case FIRSTWALKTOWALK:
	        		pr = q*walkspeed1;
	        		break;
	        	case SECPICKTOPICK:
	        		pr = p*pickspeed2;
	        		break;
	        	case SECPICKTOWALK:
	        		pr = q*pickspeed2;
	        		break;
	        	default:
	        		break;
        	}
        }
        else if (i[1]==1 && i[2]==1){
        	switch (e){
	        	case FIRSTPICKTOPICK:
	        		pr = p*pickspeed1;
	        		break;
	        	case FIRSTPICKTOWALK:
	        		pr = q*pickspeed1;
	        		break;
	        	case SECPICKTOPICK:
	        		pr = p*pickspeed2;
	        		break;
	        	case SECPICKTOWALK:
	        		pr = q*pickspeed2;
	        		break;
	        	default:
	        		break;
	    	}
        } 
        return pr;
    }
    
    @Override
    public String description() {
        String stg = "NarrowCont-Aisle-OPS \n\nThere are " + picklocs + " pick locations in the system\n\n";
        stg += "There are 2 Pickers in the system.\n";
        stg += "The pick and walk times are exponential RVs.\n";
        stg += "The average pick speed of picker 1 is " + pickspeed1 + " locations per minute \n";
        stg += "The average walk speed of picker 1 is " + walkspeed1 + " locations per minute\n";
        stg += "The average pick speed of picker 2 is " + pickspeed2 + " locations per minute \n";
        stg += "The average walk speed of picker 2 is " + walkspeed2 + " locations per minute\n";
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
         double bloq1, bloq2;
         bloq2=busy[0]+busy[1]+busy[2]+busy[3];
         bloq1=busy[4*picklocs]+busy[4*picklocs+1]+busy[4*picklocs+2]+busy[4*picklocs+3];
         out.printf("    Number of Pick Locations   = %6d ", this.picklocs);
         out.printf("\n    Number of Pickers          = %6d ", 2);
         out.printf("\n    Picking Probability        = %2.6f ", this.p);
         out.printf("\n    Percentage of Time Picker 1 is blocked    = %2.6f ", bloq1*100);
         out.printf("\n    Percentage of Time Picker 2 is blocked    = %2.6f ", bloq2*100);
         out.println();
        }    
         catch (NotUnichainException e) {
            out.println(e);
        }

        return 0;
    } 
    
    /**
     * Main method.
     * 
     * @param a
     *            Not used.
     */
    public static void main(String[]a) {
        PickBlockNarrowCont OPS = new PickBlockNarrowCont();
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

} // End of PickBlockNarrowCont class

/**
 * A State in a Picking System. This class represents a State in a
 * Narrow-Aisle-OPS with 2 pickers in the system. Every picker picks
 * at each location with probability p. The State consist of a scalar
 * denoting the distance between pickers, and a scalar with value either
 * 0 or 1 for each picker depending on whether they just walked or picked. 
 * 
 */

class PickBlockNarrowContState extends PropertiesState {

	int picklocs; 
    /**
     * Constructor for the PickBlockNarrowContState class
     * 
     * @param position
     *            0: distance between pickers
     *            1: last activity of picker 1
     *            2: last activity of picker 2
     */
    public PickBlockNarrowContState(int [] position, int picklocs) {	
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
    public PickBlockNarrowContState Trans(int [] i, PickBlockNarrowContEvent e) {
        int [] newpos;
        newpos = new int [3];
        newpos [0]= i[0] + e.getDistshift();
        newpos [1]= i[1] + e.getFirst();
        newpos [2]= i[2] + e.getSecond();
        return new PickBlockNarrowContState(newpos, picklocs);
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
        if (a==0 && c==1){stg="Picker 1 blocked picker 2, picker 2 is waiting to pick.";}
        else if (a==0 && c==0){stg="Picker 1 blocked picker 2, picker 2 is waiting to walk.";}
        else if (a==locs && b==1){stg="Picker 2 blocked picker 1, picker 1 is waiting to pick.";}
        else if (a==locs && b==0){stg="Picker 2 blocked picker 1, picker 1 is waiting to walk.";}
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

} // End of PickBlockNarrowContState class.

/**
 * Each Event characterizes the shift in distance from i to j and a posible change
 * in each picker's last activity.
 * 
 * @author Daniel Silva
 */
class PickBlockNarrowContEvent extends jmarkov.basic.Event {
    
	/** Event types. */
    public static enum Type {
        /** The first picker finishes a pick and picks*/
        FIRSTPICKTOPICK,
    	/** The first picker finishes a walk and picks*/
        FIRSTWALKTOPICK,
        /** The first picker finishes a pick and walks*/
        FIRSTPICKTOWALK,
    	/** The first picker finishes a walk and walks*/
        FIRSTWALKTOWALK,
        /** The second picker finishes a pick and picks*/
        SECPICKTOPICK,
    	/** The second picker finishes a walk and picks*/
        SECWALKTOPICK,
        /** The second picker finishes a pick and walks*/
        SECPICKTOWALK,
    	/** The second picker finishes a walk and walks*/
        SECWALKTOWALK,
        ;
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
    public PickBlockNarrowContEvent(Type type) {
    	this.type = type;
    	switch (type){
    	case FIRSTPICKTOPICK: 
    		this.distshift = 1;
    		this.first=0;
    		this.second=0;
    		break;
    	case FIRSTWALKTOPICK: 
    		this.distshift = 1;
    		this.first=1;
    		this.second=0;
    		break;
    	case FIRSTPICKTOWALK: 
    		this.distshift = 1;
    		this.first=-1;
    		this.second=0;
    		break;
    	case FIRSTWALKTOWALK: 
    		this.distshift = 1;
    		this.first=0;
    		this.second=0;
    		break;
    	case SECPICKTOPICK: 
    		this.distshift = -1;
    		this.first=0;
    		this.second=0;
    		break;
    	case SECWALKTOPICK: 
    		this.distshift = -1;
    		this.first=0;
    		this.second=1;
    		break;
    	case SECPICKTOWALK: 
    		this.distshift = -1;
    		this.first=0;
    		this.second=-1;
    		break;
    	case SECWALKTOWALK: 
    		this.distshift = -1;
    		this.first=0;
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
     * @return first picker's current activity
     */
	public int getFirst() {
		return first;
	}

    /**
     * @return second picker's current activity
     */
	public int getSecond() {
		return second;
	}

	/**
     * @return A set with all the events in the system.
     */
    public static EventsSet<PickBlockNarrowContEvent> getAllEvents() {
        EventsSet<PickBlockNarrowContEvent> eSet = new EventsSet<PickBlockNarrowContEvent>();
        eSet.add(new PickBlockNarrowContEvent(Type.FIRSTPICKTOPICK));
        eSet.add(new PickBlockNarrowContEvent(Type.FIRSTPICKTOWALK));
        eSet.add(new PickBlockNarrowContEvent(Type.FIRSTWALKTOPICK));
        eSet.add(new PickBlockNarrowContEvent(Type.FIRSTWALKTOWALK));
        eSet.add(new PickBlockNarrowContEvent(Type.SECPICKTOPICK));
        eSet.add(new PickBlockNarrowContEvent(Type.SECPICKTOWALK));
        eSet.add(new PickBlockNarrowContEvent(Type.SECWALKTOPICK));
        eSet.add(new PickBlockNarrowContEvent(Type.SECWALKTOWALK));
        return eSet;
    }

    @Override
    public String label() {
        String stg = "";
        switch (type) {
        case FIRSTPICKTOPICK:
        	stg = "The first picker finishes a pick and picks";
            break;
        case FIRSTWALKTOPICK:
            stg = "The first picker finishes a walk and picks";
            break;
        case FIRSTPICKTOWALK:
            stg = "The first picker finishes a pick and walks";
            break;
        case FIRSTWALKTOWALK:
            stg = "The first picker finishes a walk and walks";
            break;
        case SECPICKTOPICK:
        	stg = "The second picker finishes a pick and picks";
            break;
        case SECWALKTOPICK:
            stg = "The second picker finishes a walk and picks";
            break;
        case SECPICKTOWALK:
            stg = "The second picker finishes a pick and walks";
            break;
        case SECWALKTOWALK:
            stg = "The second picker finishes a walk and walks";
            break;    
        }
        return stg;
    }
}
