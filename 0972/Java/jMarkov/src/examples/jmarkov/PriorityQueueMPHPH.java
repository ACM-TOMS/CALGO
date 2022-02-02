package examples.jmarkov;

import jmarkov.GeomProcess;
import jmarkov.GeomRelState;
import jmarkov.MarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.StatesSet;
import jphase.DenseContPhaseVar;
import jphase.PhaseVar;

import static examples.jmarkov.PriorityQueueMPHPHEvent.Type.ARRIVAL_HI;
import static examples.jmarkov.PriorityQueueMPHPHEvent.Type.SERVICE_END_HI;
import static examples.jmarkov.PriorityQueueMPHPHEvent.Type.SERVICE_PHASECHG_HI;
import static examples.jmarkov.PriorityQueueMPHPHEvent.Type.ARRIVAL_LOW;
import static examples.jmarkov.PriorityQueueMPHPHEvent.Type.SERVICE_END_LOW; 
import static examples.jmarkov.PriorityQueueMPHPHEvent.Type.SERVICE_PHASECHG_LOW; 


public class PriorityQueueMPHPH extends GeomProcess<PriorityQueueMPHPHState, PriorityQueueMPHPHEvent>{

	/**
	 * High-priority arrival rate
	 */
	double lambda_hi;
	
	/**
	 * Low-priority arrival rate
	 */
    double lambda_low;
    
    /**
     * High-priority service time PH variable
     */
    PhaseVar servTime_hi;
    
    /**
     * Low-priority service time PH variable
     */
    PhaseVar servTime_low;
    
    /**
     * 
     */
    int bufferCapacity;  
    
    
    /**
     * 
     * @param lambda_hi
     * @param lambda_low
     * @param servTime_hi
     * @param servTime_low
     */
    public  PriorityQueueMPHPH(double lambda_hi, double lambda_low, PhaseVar servTime_hi, PhaseVar servTime_low, int bufferCapacity) {
        super(new PriorityQueueMPHPHState(0,0,0), 
        		PriorityQueueMPHPHEvent.getAllEvents(servTime_hi, servTime_low));
        this.lambda_hi = lambda_hi;
        this.lambda_low = lambda_low;
        this.servTime_hi = servTime_hi;
        this.servTime_low = servTime_low;
        this.bufferCapacity = bufferCapacity; 
    }

    
    /**
     * Used by GUI.
     */
    public  PriorityQueueMPHPH() {
        this(1.0, 0.5, 
        		DenseContPhaseVar.HyperExpo(new double[] { 5.0, 8.0 }, new double[] { 0.5, 0.5 }),
        		DenseContPhaseVar.HyperExpo(new double[] { 3.0, 5.0 }, new double[] { 0.5, 0.5 }), 
        		10 );
    }
	
    
    @Override
    public boolean active(PriorityQueueMPHPHState state, int absLevel, PriorityQueueMPHPHEvent event) {

        boolean result = false;
        switch (event.eventType) {
	        case ARRIVAL_HI:
	        	if ( state.getNumberHiJobs() < bufferCapacity )
	        		result = true;
            	break;
	        case SERVICE_END_HI:
	            result =  (state.getServiceType()==1 && state.getServicePhase() == event.eventPhase);
	            //check if service completion rate in this phase is non-zero
	            result = result && servTime_hi.getMat0().get(state.getServicePhase()-1) > 0;
	            break;
	        case SERVICE_PHASECHG_HI:
	            result =  (state.getServiceType()==1 && state.getServicePhase() == event.eventPhase);
	            // check if phase change probability is non-zero
	            result = result && servTime_hi.getMat0().get(state.getServicePhase()-1) < servTime_hi.getMatrix().get(state.getServicePhase()-1, event.eventPhase-1);
	            break;
			case ARRIVAL_LOW:
				result = true; 
				break;
			case SERVICE_END_LOW:
				result =  (state.getServiceType()==2 && state.getServicePhase() == event.eventPhase);
				//check if service completion rate in this phase is non-zero
				result = result && servTime_low.getMat0().get(state.getServicePhase()-1) > 0;
				break;
			case SERVICE_PHASECHG_LOW:
				result =  (state.getServiceType()==2 && state.getServicePhase() == event.eventPhase);
				//check if phase change probability is non-zero
				result = result && servTime_low.getMat0().get(state.getServicePhase()-1) < servTime_low.getMatrix().get(state.getServicePhase()-1, event.eventPhase-1);
				break;
        }
        return result;
    }

    
    @Override
    public GeomRelState<PriorityQueueMPHPHState>[] dests(PriorityQueueMPHPHState state,
            											int absLevel, PriorityQueueMPHPHEvent event) {
    	
        StatesSet<GeomRelState<PriorityQueueMPHPHState>> destStates = new StatesSet<GeomRelState<PriorityQueueMPHPHState>>();

        int curPhase = state.getServicePhase(); 
        int curNumHiJobs = state.getNumberHiJobs(); 
        int curServiceType = state.getServiceType(); 
        int rLevel = 0;

        switch (event.eventType) {
	        case ARRIVAL_LOW:
	            rLevel = +1;
	            if (absLevel == 0 && curNumHiJobs == 0) 
            		// low-priority arrival in an empty system: level increases and starts service 
	                addDestsServiceEnd(rLevel, destStates, servTime_low, curNumHiJobs, 2);
	            else // low-priority arrival in non-empty system: only modifies the level by 1
	                destStates.add(
	                			new GeomRelState<PriorityQueueMPHPHState>(
	                				new PriorityQueueMPHPHState(curNumHiJobs, curPhase, curServiceType), rLevel));
	            break;
	            
	        case SERVICE_END_LOW:
	            rLevel = -1;
	            if (absLevel == 1){
	            	if( curNumHiJobs == 0){ // low-priority service completion that leaves the system empty
	            		destStates.add(new GeomRelState<PriorityQueueMPHPHState>(
                						new PriorityQueueMPHPHState(0, 0, 0))
	            						);
	            	}else{	// low-priority service completion that allows a high-priority job to start service 
	            		addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs, 1);
	            	}
	            }
	            else{
	            	if( curNumHiJobs == 0){ // low-priority service completion, a new low-priority job starts service 
	            		addDestsServiceEnd(rLevel, destStates, servTime_low, curNumHiJobs, 2);
	            	}else{ // low-priority service completion, a new high-priority job starts service
	            		addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs, 1);
	            	}
	            }
	            break;
	            
	        case SERVICE_PHASECHG_LOW:
	            addDestsServiceChg(rLevel, destStates, servTime_low, curNumHiJobs, curPhase, curServiceType); 
	            break;
	            
			case ARRIVAL_HI:
				if (absLevel == 0 && curNumHiJobs == 0) 
            		// hi-priority arrival in an empty system: number of high jobs increases and starts service 
	                addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs+1, 1);
	            else // hi-priority arrival in non-empty system: only increases the number of high-priority jobs
	                destStates.add(
	                			new GeomRelState<PriorityQueueMPHPHState>(
	                				new PriorityQueueMPHPHState(curNumHiJobs+1, curPhase, curServiceType), rLevel));
				break;
			case SERVICE_END_HI:
				if (absLevel == 0){
	            	if( curNumHiJobs == 1){ // hi-priority service completion that leaves the system empty
	            		destStates.add(new GeomRelState<PriorityQueueMPHPHState>(
                						new PriorityQueueMPHPHState(0, 0, 0)) );
	            	}else{	// hi-priority service completion that allows a high-priority job to start service 
	            		addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs-1, 1);
	            	}
	            }
	            else{
	            	if( curNumHiJobs == 1){ // hi-priority service completion, a new low-priority job starts service 
	            		addDestsServiceEnd(rLevel, destStates, servTime_low, curNumHiJobs-1, 2);
	            	}else{ // hi-priority service completion, a new high-priority job starts service
	            		addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs-1, 1);
	            	}
	            }
				
				break;
			case SERVICE_PHASECHG_HI:
	            addDestsServiceChg(rLevel, destStates, servTime_hi, curNumHiJobs, curPhase, curServiceType); 
	            break;
				
        }
        return destStates.toStateArray();
    }

    
    
    /**
     * Adds to the state set all destinations generated by a service initiation. 
     * @param relLevel: change in the level
     * @param destStates
     * @param serVar: service-time PH variable of the new job 
     */
    void addDestsServiceEnd(int relLevel, StatesSet<GeomRelState<PriorityQueueMPHPHState>> destStates, PhaseVar serVar, int numHiJobs, int serType) {
        double alpha[] = serVar.getVectorArray(); // initial probability distribution 
        int m = serVar.getNumPhases();
        for (int n = 0; n < m; n++) {
            if (alpha[n] > 0) {
                destStates.add( 
                		new GeomRelState<PriorityQueueMPHPHState>(new PriorityQueueMPHPHState(numHiJobs, n+1, serType), relLevel)
                		);
            }
        }
    }
    
    /**
     * Adds to the state set all destinations generated by a service phase transition. 
     * @param relLevel: change in the level 
     * @param destStates
     * @param serVar: service-time PH variable of the new job 
     */
    void addDestsServiceChg(int relLevel, StatesSet<GeomRelState<PriorityQueueMPHPHState>> destStates, PhaseVar serVar, int numHiJobs, int currPhase, int curServiceType) {
        double transMatrix[][] = serVar.getMatrixArray(); // sub-generator matrix 
        currPhase = currPhase - 1; 
        int m = serVar.getNumPhases();
        for (int n = 0; n < m; n++) {
            if (transMatrix[currPhase][n] > 0) {
                destStates.add( 
                		new GeomRelState<PriorityQueueMPHPHState>(new PriorityQueueMPHPHState(numHiJobs, n+1, curServiceType), relLevel)
                		);
            }
        }
    }
    
    
    
    @Override
    public double rate(PriorityQueueMPHPHState curState, int curLevel, 
    						PriorityQueueMPHPHState destState, int destLevel, PriorityQueueMPHPHEvent event) { 
    	
		double rate = -1; 
		int curPhase = curState.getServicePhase();
		int curNumHiJobs = curState.getNumberHiJobs();
		int newPhase = destState.getServicePhase();
        switch (event.eventType) {
	        case ARRIVAL_LOW:
	            if (curLevel == 0 && curNumHiJobs == 0){
            		// low-priority arrival in an empty system: level increases and starts service
	            	double initVector[] = servTime_low.getVectorArray(); 
	                rate = lambda_low*initVector[newPhase-1];
	            }
	            else{ // low-priority arrival in non-empty system: only modifies the level by 1
	            	rate = lambda_low; 
	            }
	            break;
	            
	        case SERVICE_END_LOW:
	            if (curLevel == 1){
	            	if( curNumHiJobs == 0){ // low-priority service completion that leaves the system empty
	            		double exitVector[] = servTime_low.getMat0Array(); 
	            		rate = exitVector[curPhase-1]; 
	            	}else{	// low-priority service completion that allows a high-priority job to start service 
	            		double exitVector[] = servTime_low.getMat0Array();
	            		double initVector[] = servTime_hi.getVectorArray();
	            		rate = exitVector[curPhase-1] * initVector[newPhase-1];
	            	}
	            }
	            else{
	            	if( curNumHiJobs == 0){ // low-priority service completion, a new low-priority job starts service 
	            		double exitVector[] = servTime_low.getMat0Array();
	            		double initVector[] = servTime_low.getVectorArray();
	            		rate = exitVector[curPhase-1] * initVector[newPhase-1];
	            	}else{ // low-priority service completion, a new high-priority job starts service
	            		double exitVector[] = servTime_low.getMat0Array();
	            		double initVector[] = servTime_hi.getVectorArray();
	            		rate = exitVector[curPhase-1] * initVector[newPhase-1];
	            	}
	            }
	            break;
	            
	        case SERVICE_PHASECHG_LOW:
	        	double transMatrixLow[][] = servTime_low.getMatrixArray(); 
	        	rate = transMatrixLow[curPhase-1][newPhase-1]; 
	            break;
	            
			case ARRIVAL_HI:
				if (curLevel == 0 && curNumHiJobs == 0){ 
            		// hi-priority arrival in an empty system: number of high jobs increases and starts service
					double initVector[] = servTime_hi.getVectorArray(); 
					rate = lambda_hi*initVector[newPhase-1];
				}
	            else // hi-priority arrival in non-empty system: only increases the number of high-priority jobs
	            	rate = lambda_hi;
				break;
			case SERVICE_END_HI:
				if (curLevel == 0){
	            	if( curNumHiJobs == 1){ // hi-priority service completion that leaves the system empty
	            		double exitVector[] = servTime_hi.getMat0Array(); 
	            		rate = exitVector[curPhase-1];
	            	}else{	// hi-priority service completion that allows a high-priority job to start service 
	            		double exitVector[] = servTime_hi.getMat0Array();
	            		double initVector[] = servTime_hi.getVectorArray();
	            		rate = exitVector[curPhase-1] * initVector[newPhase-1];
	            	}
	            }
	            else{
	            	if( curNumHiJobs == 1){ // hi-priority service completion, a new low-priority job starts service 
	            		double exitVector[] = servTime_hi.getMat0Array();
	            		double initVector[] = servTime_low.getVectorArray();
	            		rate = exitVector[curPhase-1] * initVector[newPhase-1];
	            	}else{ // hi-priority service completion, a new high-priority job starts service
	            		double exitVector[] = servTime_hi.getMat0Array();
	            		double initVector[] = servTime_hi.getVectorArray();
	            		rate = exitVector[curPhase-1] * initVector[newPhase-1];
	            	}
	            }
				
				break;
			case SERVICE_PHASECHG_HI:
				double transMatrixHi[][] = servTime_hi.getMatrixArray(); 
	        	rate = transMatrixHi[curPhase-1][newPhase-1];  
	            break;
				
        }
        return rate;
    }

    
    
    
    
	public static void main(String[] a) {
        double lambda_hi = 3;
        double lambda_low = 2.0;
        PhaseVar servTime_hi = DenseContPhaseVar.HyperExpo(new double[] { 5.0, 8.0 },
                new double[] { 0.5, 0.5 });
        PhaseVar servTime_low = DenseContPhaseVar.HyperExpo(new double[] { 3.0, 5.0 },
                new double[] { 0.5, 0.5 });
        int bufferCapacity = 10; 
        
        PriorityQueueMPHPH model = new PriorityQueueMPHPH(lambda_hi, lambda_low, servTime_hi, servTime_low, bufferCapacity);
        
        
        model.showGUI();
        model.setDebugLevel(5);
        model.generate();
        model.printAll();
    }


	@Override
	public String description() {
		String desc = "Priority Queue\n" + "lambda high: " + this.lambda_hi + "\nlambda low: " + this.lambda_low + "\n";
		desc = "Service time low:" + this.servTime_low.toString() +"\n";
		desc = "Service time high:" + this.servTime_hi.toString() +"\n";
		return desc;
	}
}






class PriorityQueueMPHPHState extends PropertiesState {
	
	
	/**
	 * 
	 * @param numberHiJobs
	 * @param servicePhase
	 * @param serviceType: type of service in execution: 0 (no service), 1 (high priority), 2 (low priority)
	 */
	public PriorityQueueMPHPHState(int numberHiJobs, int servicePhase, int serviceType) {
		super(3);
        setProperty(0, numberHiJobs);
        setProperty(1, servicePhase);
        setProperty(2, serviceType);
    }
	
	
	public int getNumberHiJobs() {
        return this.prop[0];
    }
	
	public int getServicePhase() {
        return this.prop[1];
    }
	
	public int getServiceType() {
        return this.prop[2];
    }
	
	
	@Override
	public void computeMOPs(MarkovProcess mp) {
		setMOP(mp,"Number High Jobs", getNumberHiJobs());
		setMOP(mp,"High Jobs Blocking Probability", 
				getNumberHiJobs() == ((PriorityQueueMPHPH)mp).bufferCapacity ? 1 : 0);
		setMOP(mp,"Utilization High Jobs ", 
				getServiceType() == 1 ? 1 : 0);
		setMOP(mp,"Utilization Low Jobs ", 
				getServiceType() == 2 ? 1 : 0);
		setMOP(mp,"Utilization ", 
				getServiceType() > 0 ? 1 : 0);
		
	}
	
}

class PriorityQueueMPHPHEvent extends Event {
	
	
	public enum Type {
        /** High-priority arrivals. */
        ARRIVAL_HI,
        /** High-priority service completion. */
        SERVICE_END_HI,
        /** High-priority service phase change. */
        SERVICE_PHASECHG_HI,
        /** Low-priority arrivals. */
        ARRIVAL_LOW,
        /** Low-priority service completion. */
        SERVICE_END_LOW,
        /** Low-priority service phase change. */
        SERVICE_PHASECHG_LOW
        
    }
	
	/**
	 * Type of event
	 */
	Type eventType; 
	
	/**
	 * Phase where event occurs
	 */
	int eventPhase;
	
	
	/** Arrival event */
	PriorityQueueMPHPHEvent(Type eventType) {
		if( eventType == ARRIVAL_HI || eventType == ARRIVAL_LOW){
			this.eventType = eventType;
		}
    }

    /** Service event */
    PriorityQueueMPHPHEvent(Type type, int phase) {
        this.eventType = type;
        this.eventPhase = phase;
    }

	
	static EventsSet<PriorityQueueMPHPHEvent> getAllEvents(PhaseVar servTime_hi, PhaseVar servTime_low) {
		
        EventsSet<PriorityQueueMPHPHEvent> E = new EventsSet<PriorityQueueMPHPHEvent>();
        
        // high-priority arrival event 
        E.add(new PriorityQueueMPHPHEvent(ARRIVAL_HI));
        
        // high-priority service completion and service change in each phase 
        int numPhases_hi = servTime_hi.getNumPhases();
        for (int n = 1; n <= numPhases_hi; n++) {
            // high-priority service completion in phase n
            E.add(new PriorityQueueMPHPHEvent(SERVICE_END_HI, n));
            E.add(new PriorityQueueMPHPHEvent(SERVICE_PHASECHG_HI, n));
        }
        
        // low-priority arrival event 
        E.add(new PriorityQueueMPHPHEvent(ARRIVAL_LOW));
        
        // high-priority service completion and service change in each phase 
        int numPhases_low = servTime_low.getNumPhases();
        for (int n = 1; n <= numPhases_low; n++) {
            // high-priority service completion in phase n
            E.add(new PriorityQueueMPHPHEvent(SERVICE_END_LOW, n));
            E.add(new PriorityQueueMPHPHEvent(SERVICE_PHASECHG_LOW, n));
        }
        
        return E;
    }
	
	@Override
    public String label() {
        String stg = "";
        switch (eventType) {
        case ARRIVAL_HI:
            stg = "Arrival Hi";
            break;
        case SERVICE_END_HI:
            stg = "Service End Hi (" + eventPhase + ")";
            break;
        case SERVICE_PHASECHG_HI:
            stg = "Service Ph-Ch Hi (" + eventPhase + ")";
            break;
        case ARRIVAL_LOW:
            stg = "Arrival Low";
            break;
        case SERVICE_END_LOW:
            stg = "Service End Low (" + eventPhase + ")";
            break;
        case SERVICE_PHASECHG_LOW:
            stg = "Service Ph-Ch Low (" + eventPhase + ")";
            break;
        }
        return stg;
    }
	
}