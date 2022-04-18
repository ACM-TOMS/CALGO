package examples.jmarkov;

import java.io.PrintWriter;
import java.io.StringWriter;
import examples.jmarkov.PickBlockWideContKPickers_V3Event.Type;
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
 * n picklocs and exponential picking and walking times, each time is different 
 * for each picker. Each picker picks at a given location with probability p, 
 * or walks by with probability q=1-p.
 * The pickers are blocked when one reaches the other and both need to pick at the 
 * same location. We assume a circular warehouse without time spent at the I/O point.
 * 
 * @author Daniel F. Silva (2007)
 * 
 */

public class PickBlockWideContKPickers_V3 extends SimpleMarkovProcess<PickBlockWideContKPickers_V3State, PickBlockWideContKPickers_V3Event> {
    private double p;
    private int picklocs,k;
    private double pickspeed, walkspeed;
    
 
   /**
     * General constructor.
     * 
     * @param p
     *            Pick Probability.
     * @param picklocs
     *            Number of Pick Locations
     * @param pickspeed
     *            Rate of picking picks per minute for picker 1.
     * @param walkspeed
     *            Rate of walking past locations locations per minute for picker 1.
     *@param k
     *            Number of pickers.
     * 
     */
    
    public PickBlockWideContKPickers_V3(double p, int picklocs, int k, double pickspeed, double walkspeed ,int[]a ) {
    	super(new PickBlockWideContKPickers_V3State(a,picklocs,k), PickBlockWideContKPickers_V3Event.getAllEvents(k));
        this.p = p;
        this.picklocs = picklocs;
        this.pickspeed = pickspeed;
        this.walkspeed = walkspeed;
        this.k=k;

    }
    
    /**
     * Default Constructor used by GUI
     */
    public PickBlockWideContKPickers_V3() {
    	this(1, 20, 3, 1, 10, new int[] {0,1,-1,0,0,0});
    }
    

    /**
     * Determine the active events.
     */

    @Override
    public boolean active(PickBlockWideContKPickers_V3State i, PickBlockWideContKPickers_V3Event e) {
        boolean result = false; 
        boolean blockable = false;
        int act = i.getState(e.getPicker()+k), pick=e.getPicker(),distmin=k;
        int []dude=new int[k];
        for (int j = 0; j<k; j++){
        	if ((Math.abs(i.getState(j)-(i.getState(pick)+1))==0||Math.abs(i.getState(j)-(i.getState(pick)+1))==picklocs) && j!=pick){distmin=0;dude[j]=j;}
        	else dude[j]=k;
        }
        for (int j = 0; j<k; j++){
        	if (dude[j]<k){if (distmin==0 && i.getState(dude[j]+k)>=1)blockable=true;}
        }
        if (act==2){
	      	result = false;
	    }
        else if (act==0){
	        switch (e.getType()){
	        case WALKTOPICK:
	        	if (blockable==true) result =false;
	        	else result = true; break;
	        case WALKTOWALK:
	        	result = true; break;
	        case WALKTOBLOCK:
	        	if (blockable==true)
	        		result=true;break;
	        default:
	        	result = false; break;
	     	}
	    }
        else if (act==1){
	        switch (e.getType()){
	        case PICKTOPICK:
	        	if (blockable==true) result =false;
	        	else result = true; break;
	        case PICKTOWALK:
	        	result = true; break;
	        case PICKTOBLOCK:
	        	if (blockable==true)
	        		result=true;break;
	        default:
	        	result = false; break;
	     	}
	    } 
        return result;
    }

    @Override
    public States<PickBlockWideContKPickers_V3State> dests(PickBlockWideContKPickers_V3State i, PickBlockWideContKPickers_V3Event e) {
        int [] origin = i.getProperties();
        return new StatesSet<PickBlockWideContKPickers_V3State>(i.Trans(origin,k, e));
    }

    /**
     * Returns the transition probability from State i to State j, when event e occurs.
     */

    @Override
    public double rate(PickBlockWideContKPickers_V3State i, PickBlockWideContKPickers_V3State j, PickBlockWideContKPickers_V3Event e) {
    	if (this.active(i,e)==false) 
    		return 0;
        double pr=0.0, q=1-p;
        boolean blockable = false;
        int act = i.getState(e.getPicker()+k), pick=e.getPicker(),distmin=k;
        int []dude=new int[k];
        for (int l = 0; l<k; l++){
        	if ((Math.abs(i.getState(l)-(i.getState(pick)+1))==0||Math.abs(i.getState(l)-(i.getState(pick)+1))==picklocs) && l!=pick){distmin=0;dude[l]=l;}
        	else dude[l]=k;
        }
        for (int l = 0; l<k; l++){
        	if (dude[l]<k){if (distmin==0 && i.getState(dude[l]+k)>=1)blockable=true;}
        }  
        if (act==2){
        	pr=0;
	    }
        else if (act==0){
	        switch (e.getType()){
	        case WALKTOPICK:
	        	if (blockable==true) pr=0;
	        	else pr=p*walkspeed; break;
	        case WALKTOWALK:
	        	pr=q*walkspeed; break;
	        case WALKTOBLOCK:
	        	if (blockable==true)
	        		pr=p*walkspeed;break;
	        default:
	        	pr=0; break;
	     	}
	    }
        else if (act==1){
	        switch (e.getType()){
	        case PICKTOPICK:
	        	if (blockable==true) pr=0;
	        	else pr=p*pickspeed; break;
	        case PICKTOWALK:
	        	pr=q*pickspeed; break;
	        case PICKTOBLOCK:
	        	if (blockable==true)
	        		pr=p*pickspeed;break;
	        default:
	        	pr=0; break;
	     	}
	    }         return pr;
    }
    
    @Override
    public String description() {
        String stg = "WideCont-Aisle-OPS \n\nThere are " + picklocs + " pick locations in the system\n\n";
        stg += "There are "+k+" Pickers in the system.\n";
        stg += "The pick and walk times are exponential RVs.\n";
        stg += "The average pick speed of a picker is " + pickspeed + " locations per minute \n";
        stg += "The average walk speed of a picker is " + walkspeed  + " locations per minute\n";
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
            double block=0;
            State [] c;
            StatesSet a=getStates();

            c=a.toStateArray();
            
            for (int i = 0 ; i<k ; i++){
   	         for (int j=0; j<getNumStates(); j++ ){
   	        	 if(c[j].getMOP(k+i)==2){
   	        	 	block+=busy[j];
   	        	 }
   	         }
            }
         out.printf("    Numero de Pick Locations   = %6d ", this.picklocs);
         out.printf("\n    Numero de Pickers          = %6d ", k);
         out.printf("\n    Picking Probability        = %2.6f ", this.p);
         out.printf("\n    Porcentaje de Tiempo de Bloqueo de cada Picker     = %2.6f ", block*100/k);
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

        PickBlockWideContKPickers_V3 OPS = new PickBlockWideContKPickers_V3();
        OPS.setDebugLevel(0);

        OPS.setMaxStates(8000);
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

} // End of PickBlockWideContKPickers_V3 class


class PickBlockWideContKPickers_V3State extends PropertiesState {
/**
 * A State in a Picking System. This class represents a State in a
 * WideCont-Aisle-OPS with 2 pickers in the system. Every picker picks
 * at each location with probability p. The State consist of a scalar
 * denoting the distance between pickers, and a scalar with value either
 * 0 or 1 for each picker depending on whether they just walked or picked. 
 * 
 */

int picklocs, k; 
/**
 * Constructor for the PickBlockWideContKPickers_V3State class
 * 
 * @param position
 *            0 to k-1: distance between each picker and the next
 *            k to 2k-1: last activity of each picker
 */
public PickBlockWideContKPickers_V3State(int [] position, int picklocs,int k) {	
    super(position);
    for (int i = 0 ; i <= (2*k-1) ; i++){
    	setProperty(i, position [i]);	
    }
    this.picklocs=picklocs;
    this.k=k;

}

/**
 * Returns the the current state of a given position.
 * 
 */
public int getState(int i) {
    return getProperty(i);
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
 *            Original state
 * @param e
 *            Event
 * @return The new state
 */
public PickBlockWideContKPickers_V3State Trans(int [] i, int k, PickBlockWideContKPickers_V3Event e) {
    int [] newpos;
    newpos = new int[2*k];
    if (e.getPicker()!=0){
	    for (int j=1; j<k ; j++){
	    	if (j==e.getPicker()){
	        	newpos[j]=i[j]+1;
	    	}
	    	else newpos[j]=i[j];
	    }
    }
    else{
	    for (int j=1; j<k ; j++){
	       	newpos[j]=i[j]-1;
	    }
    }
    for (int j=1; j<k ; j++){
    	if (picklocs%2==0){
	    	if (newpos[j]==(-picklocs/2)) newpos[j]=picklocs/2;
	    	else if (newpos[j]==(picklocs/2+1)) newpos[j]=(-picklocs/2+1);
	    }
	    else if (picklocs%2==1){
	    	if (newpos[j]==(-picklocs/2-1)) newpos[j]=picklocs/2;
	    	else if (newpos[j]==(picklocs/2+1)) newpos[j]=(-picklocs/2);
	    }
    }
    for (int j=k; j<=(2*k-1) ; j++){
    	if (j==k+e.getPicker() && i[j]!=2) newpos[j] = i[j] + e.getShift();
    	else newpos[j]=i[j];
    }
    int pick=e.getPicker();
    int []dude=new int[k],distmin=new int[k];
    for (int j=0; j<k ; j++){dude[j]=k;distmin[j]=k;}
    for (int j=0; j<k ; j++){ 
        for (int r =0; r<k; r++){
        	if (Math.abs(i[r]-i[j])==0 && r!=j){
        		distmin[j]=0;
        		dude[j]=r;
            	if (i[j+k]==2 && distmin[j]==0 && dude[j]==pick && i[dude[j]+k]==1){newpos[j+k]=1;break;}
        	}
        }
        if (i[j+k]==2 && distmin[j]==0 && dude[j]==pick && i[dude[j]+k]==1){break;}
    }
    return new PickBlockWideContKPickers_V3State(newpos, picklocs, k);
}

/**
 * Describes the States
 * 
 * @return a String description of the State.
 */
@Override
public String description() {
    String stg = "";
    for (int i=0; i<k ; i++){
    	stg+="Picker "+(i+1)+ " is ";
    	if (getState(i+k)==0)stg+="walking, ";
    	else stg+= "picking, ";
    	if (i==k-1)stg+="the Distance between Picker " + (i+1) + " and picker " + (1) + " is " + getState(i) +". ";
    	else stg+="the Distance between Picker " + (i+1) + " and picker " + (i+2) + " is " + getState(i) +". ";
    }
    return stg;
}

} // End of PickBlockWideContKPickers_V3State class.

/**
 * Esta clase implementa eventos en un Drive Thru. Extendiende la clase
 * SimpleMarkovProcess.
 */
class PickBlockWideContKPickers_V3Event extends jmarkov.basic.Event {
    
	/** Event types. */
    public static enum Type {
        /** A picker finishes a pick and picks*/
        PICKTOPICK,
    	/** A picker finishes a walk and picks*/
        WALKTOPICK,
        /** A picker finishes a pick and walks*/
        PICKTOWALK,
    	/** A picker finishes a walk and walks*/
        WALKTOWALK,
        /** A picker moves to a blocked state*/
        PICKTOBLOCK,
        /** A picker moves to a blocked state*/
        WALKTOBLOCK,

    }

    public Type type;
	public int picker, shift;
	
	/** Create an event
	 * 
	 * @param picker
	 * 			Picker that finishes an action
	 * 
	 * @param Type
	 * 			The action that is finiched and which is begun
	 */
	
    public PickBlockWideContKPickers_V3Event(int picker, Type type) {
        this.picker = picker;
        this.type = type;
        switch (type){
	    	case PICKTOPICK: 
	    		this.shift = 0;
	    		break;
	    	case WALKTOPICK: 
	    		this.shift = 1;
	    		break;
	    	case PICKTOWALK: 
	    		this.shift = -1;
	    		break;
	    	case WALKTOWALK: 
	    		this.shift = 0;
	    		break;
	    	case PICKTOBLOCK: 
	    		this.shift = 1;
	    		break;
	    	case WALKTOBLOCK: 
	    		this.shift = 2;
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
     * @return picker that moves 
     */
    public int getPicker() {
		return picker;
	}

    /**
     * @return shift
     */
    public int getShift() {
		return shift;
	}

	/**
     * @return A set with all the events in the system.
     */
    public static EventsSet<PickBlockWideContKPickers_V3Event> getAllEvents(int k) {
        EventsSet<PickBlockWideContKPickers_V3Event> eSet = new EventsSet<PickBlockWideContKPickers_V3Event>();
        for (int i=0;i<k; i++){
	        eSet.add(new PickBlockWideContKPickers_V3Event(i, Type.PICKTOPICK));
	        eSet.add(new PickBlockWideContKPickers_V3Event(i, Type.PICKTOWALK));
	        eSet.add(new PickBlockWideContKPickers_V3Event(i, Type.WALKTOPICK));
	        eSet.add(new PickBlockWideContKPickers_V3Event(i, Type.WALKTOWALK));
	        eSet.add(new PickBlockWideContKPickers_V3Event(i, Type.PICKTOBLOCK));
	        eSet.add(new PickBlockWideContKPickers_V3Event(i, Type.WALKTOBLOCK));
	        
        }
        return eSet;
    }

    @Override
    public String label() {
        String stg = "";
        switch (type) {
        case PICKTOPICK:
        	stg = "Picker "+(picker+1)+" finishes a pick and picks";
            break;
        case WALKTOPICK:
            stg = "Picker "+(picker+1)+" finishes a walk and picks";
            break;
        case PICKTOWALK:
            stg = "Picker "+(picker+1)+" finishes a pick and walks";
            break;
        case WALKTOWALK:
            stg = "Picker "+(picker+1)+" finishes a walk and walks";
            break;
        case PICKTOBLOCK:
            stg = "Picker "+(picker+1)+" finishes a pick and gets blocked";
            break;
        case WALKTOBLOCK:
            stg = "Picker "+(picker+1)+" finishes a walk and gets blocked";
            break;
        }
        return stg;
    }
}
