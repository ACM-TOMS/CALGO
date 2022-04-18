package examples.jmarkov;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.StringTokenizer;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Vector.Norm;
import jmarkov.GeomProcess;
import jmarkov.GeomRelState;
import jmarkov.MarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.StatesSet;
import jphase.ContPhaseVar;
import jphase.DenseContPhaseVar;
import jphase.PhaseVar;
import jphase.fit.EMHyperErlangFit;
import jphase.fit.MomentsACPH2Fit;
import jphase.fit.MomentsACPHFit;

import static examples.jmarkov.PriorityQueueMPHPHPreemptEvent.Type.ARRIVAL_HI;
import static examples.jmarkov.PriorityQueueMPHPHPreemptEvent.Type.SERVICE_END_HI;
import static examples.jmarkov.PriorityQueueMPHPHPreemptEvent.Type.SERVICE_PHASECHG_HI;
import static examples.jmarkov.PriorityQueueMPHPHPreemptEvent.Type.ARRIVAL_LOW;
import static examples.jmarkov.PriorityQueueMPHPHPreemptEvent.Type.SERVICE_END_LOW; 
import static examples.jmarkov.PriorityQueueMPHPHPreemptEvent.Type.SERVICE_PHASECHG_LOW; 


public class PriorityQueueMPHPHPreempt extends GeomProcess<PriorityQueueMPHPHPreemptState, PriorityQueueMPHPHPreemptEvent>{

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
    public  PriorityQueueMPHPHPreempt(double lambda_hi, double lambda_low, PhaseVar servTime_hi, PhaseVar servTime_low, int bufferCapacity) {
        super(new PriorityQueueMPHPHPreemptState(0,0,0), 
        		PriorityQueueMPHPHPreemptEvent.getAllEvents(servTime_hi, servTime_low));
        this.lambda_hi = lambda_hi;
        this.lambda_low = lambda_low;
        this.servTime_hi = servTime_hi;
        this.servTime_low = servTime_low;
        this.bufferCapacity = bufferCapacity; 
    }

    
    /**
     * Used by GUI.
     */
    public  PriorityQueueMPHPHPreempt() {
        this(1.0, 0.5, 
        		DenseContPhaseVar.HyperExpo(new double[] { 5.0, 8.0 }, new double[] { 0.5, 0.5 }),
        		DenseContPhaseVar.HyperExpo(new double[] { 3.0, 5.0 }, new double[] { 0.5, 0.5 }), 
        		10 );
    }
	
    
    @Override
    public boolean active(PriorityQueueMPHPHPreemptState state, int absLevel, PriorityQueueMPHPHPreemptEvent event) {

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
	            result = result && servTime_hi.getMat0().get(state.getServicePhase()-1) < -servTime_hi.getMatrix().get(state.getServicePhase()-1, state.getServicePhase()-1);
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
				result = result && servTime_low.getMat0().get(state.getServicePhase()-1) < -servTime_low.getMatrix().get(state.getServicePhase()-1, state.getServicePhase()-1);
				break;
        }
        return result;
    }

    
    @Override
    public GeomRelState<PriorityQueueMPHPHPreemptState>[] dests(PriorityQueueMPHPHPreemptState state,
            											int absLevel, PriorityQueueMPHPHPreemptEvent event) {
    	
        StatesSet<GeomRelState<PriorityQueueMPHPHPreemptState>> destStates = new StatesSet<GeomRelState<PriorityQueueMPHPHPreemptState>>();

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
	                			new GeomRelState<PriorityQueueMPHPHPreemptState>(
	                				new PriorityQueueMPHPHPreemptState(curNumHiJobs, curPhase, curServiceType), rLevel));
	            break;
	            
	        case SERVICE_END_LOW:
	            rLevel = -1;
	            if (absLevel == 1){
	            	if( curNumHiJobs == 0){ // low-priority service completion that leaves the system empty
	            		destStates.add(new GeomRelState<PriorityQueueMPHPHPreemptState>(
                						new PriorityQueueMPHPHPreemptState(0, 0, 0))
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
	            else if (absLevel > 0 && curNumHiJobs == 0)  // hi-priority arrival in non-empty system with low priority job in service: preemts job in service, starts service and increases the number of high-priority jobs
	            	addDestsServiceEnd(rLevel, destStates, servTime_hi, curNumHiJobs+1, 1);
	            else // hi-priority arrival in non-empty system with hi priority job in service: increases the number of high-priority jobs
	                destStates.add(
	                			new GeomRelState<PriorityQueueMPHPHPreemptState>(
	                				new PriorityQueueMPHPHPreemptState(curNumHiJobs+1, curPhase, curServiceType), rLevel));
				break;
			case SERVICE_END_HI:
				if (absLevel == 0){
	            	if( curNumHiJobs == 1){ // hi-priority service completion that leaves the system empty
	            		destStates.add(new GeomRelState<PriorityQueueMPHPHPreemptState>(
                						new PriorityQueueMPHPHPreemptState(0, 0, 0)) );
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
    void addDestsServiceEnd(int relLevel, StatesSet<GeomRelState<PriorityQueueMPHPHPreemptState>> destStates, PhaseVar serVar, int numHiJobs, int serType) {
        double alpha[] = serVar.getVectorArray(); // initial probability distribution 
        int m = serVar.getNumPhases();
        for (int n = 0; n < m; n++) {
            if (alpha[n] > 0) {
                destStates.add( 
                		new GeomRelState<PriorityQueueMPHPHPreemptState>(new PriorityQueueMPHPHPreemptState(numHiJobs, n+1, serType), relLevel)
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
    void addDestsServiceChg(int relLevel, StatesSet<GeomRelState<PriorityQueueMPHPHPreemptState>> destStates, PhaseVar serVar, int numHiJobs, int currPhase, int curServiceType) {
        double transMatrix[][] = serVar.getMatrixArray(); // sub-generator matrix 
        currPhase = currPhase - 1; 
        int m = serVar.getNumPhases();
        for (int n = 0; n < m; n++) {
            if (transMatrix[currPhase][n] > 0) {
                destStates.add( 
                		new GeomRelState<PriorityQueueMPHPHPreemptState>(new PriorityQueueMPHPHPreemptState(numHiJobs, n+1, curServiceType), relLevel)
                		);
            }
        }
    }
    
    
    
    @Override
    public double rate(PriorityQueueMPHPHPreemptState curState, int curLevel, 
    						PriorityQueueMPHPHPreemptState destState, int destLevel, PriorityQueueMPHPHPreemptEvent event) { 
    	
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
				else if (curLevel > 0 && curNumHiJobs == 0){  // hi-priority arrival in non-empty system with low priority job in service: preemts job in service, starts service and increases the number of high-priority jobs
					double initVector[] = servTime_hi.getVectorArray();
					rate = lambda_hi*initVector[newPhase-1];
				}else // hi-priority arrival in non-empty system with hi priority job in service: increases the number of high-priority jobs
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
		
		int version = 1; //1: paper, 2: timing, 3: scalability
		
		double lambda_hi;
	    double lambda_low;
	    ContPhaseVar servTime_hi;
	    ContPhaseVar servTime_low;
	    int bufferCapacity;
	    double lambda[] = {0.05, 0.25, 0.45};
		double meanExecs[];
	    PriorityQueueMPHPHPreempt model;
		switch (version){ 
			case 1:
				// parameters used in the example in the paper
				lambda_hi = 0.2;
			    lambda_low = 0.2;
			    
			    double[] data = readTextFile("src/examples/jphase/W2.txt");
		        EMHyperErlangFit fitter_hi = new EMHyperErlangFit(data); 
			    servTime_hi = fitter_hi.fit(4);
			    
			    System.out.println("mean hi: "+servTime_hi.expectedValue());
			    System.out.println("var hi: "+servTime_hi.toString());
			    
			    
			    MomentsACPHFit fitter_low = new MomentsACPHFit(2, 6, 25);
			    //MomentsACPH2Fit fitter_low = new MomentsACPH2Fit(2, 6, 22);
				servTime_low = fitter_low.fit();
				System.out.println("mean low: "+servTime_low.expectedValue());
				System.out.println("var hi: "+servTime_low.toString());
				
					
			    bufferCapacity = 10;
			    
			    model = new PriorityQueueMPHPHPreempt(lambda_hi, lambda_low, servTime_hi, servTime_low, bufferCapacity);
			    model.setMaxStates(10000);
			    model.setDebugLevel(0);
			    model.generate();
		    	model.printMOPs();
		    	model.showGUI();
			    break;
			case 2:
				// parameters used for testing, especially  vs matlab
				//double lambda[] = {0.05, 0.25, 0.45};
				meanExecs = new double[lambda.length];
				for (int i = 0; i < lambda.length; i++){
				    
				    
				    lambda_hi = lambda[i];
				    lambda_low = lambda[i];
				    
				    /*
				    PhaseVar servTime_hi = DenseContPhaseVar.HyperExpo(new double[] { 5.0, 8.0 },
				            new double[] { 0.5, 0.5 });
				    PhaseVar servTime_low = DenseContPhaseVar.HyperExpo(new double[] { 3.0, 5.0 },
				            new double[] { 0.5, 0.5 });
				    */
				    
				    servTime_hi = DenseContPhaseVar.HyperExpo(new double[] { 1.577350269, 0.422649730},
				            new double[] { 0.788675134594813,   0.211324865405187 });
				    
				    servTime_low = DenseContPhaseVar.HyperExpo(new double[] { 1.577350269, 0.422649730 },
				            new double[] { 0.788675134594813,   0.211324865405187 });
		            
				    
				    
				    double sumAlpha = 0;
				    for(int j = 0; j < servTime_hi.getNumPhases(); j++)
				    	sumAlpha += servTime_hi.getVectorArray()[j];
				    DenseVector temp = (DenseVector)servTime_hi.getVector();
				    for(int j = 0; j < servTime_hi.getNumPhases(); j++)
				    	temp.set(j, temp.get(j)/sumAlpha);
				    servTime_hi = new DenseContPhaseVar(temp, servTime_hi.getMatrix());
				    
				    bufferCapacity = 10;
				    
				    
				    int reps = 1;
				    long meanExecTime = 0;
				    for(int rep = 0; rep < reps; rep++){
					    model = new PriorityQueueMPHPHPreempt(lambda_hi, lambda_low, servTime_hi, servTime_low, bufferCapacity);
					    model.setMaxStates(10000);
					    
					    boolean useGUI = false;
					    if(useGUI){
					        model.showGUI();
					        model.setDebugLevel(5);
					        model.generate();
					        model.printAll();}
					    else{
					    	model.setDebugLevel(0);
					    	model.generate();
					    	long startTime = System.currentTimeMillis();
					    	model.printMOPs();
					    	long stopTime = System.currentTimeMillis();
					        long elapsedTime = stopTime - startTime;
					        System.out.println("Execution time: "+elapsedTime+" ms");
					        meanExecTime += elapsedTime;  
					    }
				    }
				    meanExecTime /= reps; 
				    System.out.println("Mean execution time (ms): "+meanExecTime);
				    System.out.println("Mean execution time (s): "+(double)meanExecTime/1000);
				    meanExecs[i] = (double)meanExecTime/1000;
				     
				}
				
				for (int i = 0; i < lambda.length; i++){
					System.out.println("Mean execution time (s): "+ meanExecs[i]);
				}
				
				break;
			case 3:
				// parameters used for testing, especially for scalability and vs matlab
				//int Cset[] = {10, 20, 50, 100, 200, 1000};
				int Cset[] = {10, 20, 50, 100, 200, 500, 1000}; 
				double meanExecScala[][] = new double[lambda.length][Cset.length];
				for (int k = 0; k < Cset.length; k++){
					for (int i = 0; i < lambda.length; i++){
					    lambda_hi = lambda[i];
					    lambda_low = lambda[i];
					    servTime_hi = DenseContPhaseVar.HyperExpo(new double[] { 1.577350269, 0.422649730},
					            new double[] { 0.788675134594813,   0.211324865405187 });
					    
					    servTime_low = DenseContPhaseVar.HyperExpo(new double[] { 1.577350269, 0.422649730 },
					            new double[] { 0.788675134594813,   0.211324865405187 });
					    
					    double sumAlpha = 0;
					    for(int j = 0; j < servTime_hi.getNumPhases(); j++)
					    	sumAlpha += servTime_hi.getVectorArray()[j];
					    DenseVector temp = (DenseVector)servTime_hi.getVector();
					    for(int j = 0; j < servTime_hi.getNumPhases(); j++)
					    	temp.set(j, temp.get(j)/sumAlpha);
					    servTime_hi = new DenseContPhaseVar(temp, servTime_hi.getMatrix());
					    
					    //bufferCapacity = 10;
					    bufferCapacity = Cset[k];
					    
					    
					    int reps = 1;
					    long meanExecTime = 0;
					    for(int rep = 0; rep < reps; rep++){
						    model = new PriorityQueueMPHPHPreempt(lambda_hi, lambda_low, servTime_hi, servTime_low, bufferCapacity);
						    model.setMaxStates(10000);
						    
						    boolean useGUI = false;
						    if(useGUI){
						        model.showGUI();
						        model.setDebugLevel(5);
						        model.generate();
						        model.printAll();}
						    else{
						    	model.setDebugLevel(0);
						    	model.generate();
						    	long startTime = System.currentTimeMillis();
						    	model.printMOPs();
						    	long stopTime = System.currentTimeMillis();
						        long elapsedTime = stopTime - startTime;
						        System.out.println("Execution time: "+elapsedTime+" ms");
						        meanExecTime += elapsedTime;  
						    }
					    }
					    meanExecTime /= reps; 
					    //System.out.println("Mean execution time (ms): "+meanExecTime);
					    //System.out.println("Mean execution time (s): "+(double)meanExecTime/1000);
					    meanExecScala[i][k] = (double)meanExecTime/1000;
					     
					}
				}
				
				System.out.println("Mean execution times (s): ");
				for (int k = 0; k < Cset.length; k++){
					for (int i = 0; i < lambda.length; i++){
						//System.out.print("Mean execution time (s): "+ meanExecScala[i][k]+"\t");
						System.out.print(meanExecScala[i][k]+"\t");
					}
					System.out.println();
				}
				
				break;
								 
		}
		
		
    }


	@Override
	public String description() {
		String desc = "Priority Queue\n" + "lambda high: " + this.lambda_hi + "\nlambda low: " + this.lambda_low + "\n";
		desc = "Service time low:" + this.servTime_low.toString() +"\n";
		desc = "Service time high:" + this.servTime_hi.toString() +"\n";
		return desc;
	}
	
	private static double[] readTextFile(String fileName){
		ArrayList<Double> data = new ArrayList<Double>();
		
		try{
			FileReader archivo = new FileReader(fileName);
			BufferedReader entrada = new BufferedReader(archivo);
			String s;
			StringTokenizer str;
			System.out.println("Data file found");
											
			while (entrada.ready()){ 
				s=entrada.readLine();
				str=new StringTokenizer (s);
				while(str.countTokens() == 0){
					s=entrada.readLine();
					str=new StringTokenizer (s);
				}
				if (str.countTokens() != 1)throw new Exception ();
				data.add( new Double(Double.parseDouble(str.nextToken())) );
			}
			entrada.close();
			
		}catch(Exception e){
			System.out.println("File could not be read");
			return null;
		}
		double[] datos = new double[data.size()];
		for(int i = 0; i < data.size();i++)datos[i] = data.get(i).doubleValue();
		return datos;
	}
}






class PriorityQueueMPHPHPreemptState extends PropertiesState {
	
	
	/**
	 * 
	 * @param numberHiJobs
	 * @param servicePhase
	 * @param serviceType: type of service in execution: 0 (no service), 1 (high priority), 2 (low priority)
	 */
	public PriorityQueueMPHPHPreemptState(int numberHiJobs, int servicePhase, int serviceType) {
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
				getNumberHiJobs() == ((PriorityQueueMPHPHPreempt)mp).bufferCapacity ? 1 : 0);
		/*setMOP(mp,"Utilization High Jobs ", 
				getServiceType() == 1 ? 1 : 0);
		setMOP(mp,"Utilization Low Jobs ", 
				getServiceType() == 2 ? 1 : 0);*/
		setMOP(mp,"Utilization ", 
				getServiceType() > 0 ? 1 : 0);
		
	}
	
}

class PriorityQueueMPHPHPreemptEvent extends Event {
	
	
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
	PriorityQueueMPHPHPreemptEvent(Type eventType) {
		if( eventType == ARRIVAL_HI || eventType == ARRIVAL_LOW){
			this.eventType = eventType;
		}
    }

    /** Service event */
    PriorityQueueMPHPHPreemptEvent(Type type, int phase) {
        this.eventType = type;
        this.eventPhase = phase;
    }

	
	static EventsSet<PriorityQueueMPHPHPreemptEvent> getAllEvents(PhaseVar servTime_hi, PhaseVar servTime_low) {
		
        EventsSet<PriorityQueueMPHPHPreemptEvent> E = new EventsSet<PriorityQueueMPHPHPreemptEvent>();
        
        // high-priority arrival event 
        E.add(new PriorityQueueMPHPHPreemptEvent(ARRIVAL_HI));
        
        // high-priority service completion and service change in each phase 
        int numPhases_hi = servTime_hi.getNumPhases();
        for (int n = 1; n <= numPhases_hi; n++) {
            // high-priority service completion in phase n
            E.add(new PriorityQueueMPHPHPreemptEvent(SERVICE_END_HI, n));
            E.add(new PriorityQueueMPHPHPreemptEvent(SERVICE_PHASECHG_HI, n));
        }
        
        // low-priority arrival event 
        E.add(new PriorityQueueMPHPHPreemptEvent(ARRIVAL_LOW));
        
        // high-priority service completion and service change in each phase 
        int numPhases_low = servTime_low.getNumPhases();
        for (int n = 1; n <= numPhases_low; n++) {
            // high-priority service completion in phase n
            E.add(new PriorityQueueMPHPHPreemptEvent(SERVICE_END_LOW, n));
            E.add(new PriorityQueueMPHPHPreemptEvent(SERVICE_PHASECHG_LOW, n));
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