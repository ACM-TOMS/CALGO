package examples.jphase;

import jmarkov.MarkovProcess;
import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jphase.ContPhaseVar;
import jphase.fit.MomentsECCompleteFit;
import no.uib.cipr.matrix.Matrix;

/**
 * Example of a system with two servers (one for the cash register, and another one delivering the products) 
 * each one with a separate queue.
 */
class RealSystemState extends PropertiesState {

    /**
      *We define each state with a vector of 4 components: the state of the first server, the state of the second, the first queue
      *length and the second queue length
     */

    /**
     * Constructor for an empty system 
     */
	RealSystemState() {
		this(new int[2],0,0);
    }
   
	/**
	 * Constructor that receives status and queue lengths
	 */
	RealSystemState(int[] status, int Q1size,int Q2size) {
		super(4);
		for (int j = 0; j < 2; j++) {
			prop[j] = status[j];
		}
		prop[2] = Q1size;
		prop[3] = Q2size;
    }
 
    @Override
    public void computeMOPs(MarkovProcess mp) {
    	setMOP(mp,"Queue1 Length", getQ1Size());
    	setMOP(mp,"Queue2 Length", getQ2Size());
    }

    /**
     * @param i is position in the vector
     * @return an integer that shows the current state
     */
    public int getProp(int i) {
    	return prop[i];
    }
   
    /**
     * 
     * @return the size of the first queue 
     * 
     */
    public int getQ1Size() {
    	return prop[2];
    }
   
    /**
     * @return the size of the Second queue 
     * 
     */
    public int getQ2Size() {
        return prop[3];
   }  
   
    @Override
    public boolean isConsistent() {
    	return true;
    }

    @Override
    public String label() {
    	 String stg = "";
     	for(int k=0; k<2; k++){
     		stg+=" Server: ";
     		stg+=(k+1)+" Phase: "+prop[k]+" ";
     		}
       return stg + " Q1: " + getQ1Size()+" Q2: "+getQ2Size();
    }
    


}
/**
* 
* This class defines the events. 
* An event has two components: Type which can have 3 values 
* depending whether it represents an arrival, a departure or a Phase change, and
* the server number (1 or 2) where it is occurring.
* 
* @author Laura and Leonardo
*
*/
class RealSystemEvent extends Event{
	static final int ARRIVAL = 1;
	static final int DEPARTURE = 2;
	static final int PHASECHANGE = 3;
	int type; // ARRIVAL, DEPARTURE, PHASECHANGE
	int server;//server number

	RealSystemEvent(int type, int server) {
        this.type = type;
        this.server=server;}

	static EventsSet<RealSystemEvent> getAllEvents() {
        EventsSet<RealSystemEvent> eSet = new EventsSet<RealSystemEvent>();
        eSet.add(new RealSystemEvent(ARRIVAL,0)); // arrivals are independent of the server number, so we let the arrivals happen only for the first server 
        for(int i=0; i<2;i++){     
        	eSet.add(new RealSystemEvent(PHASECHANGE,i));// phase change for each server
        	eSet.add(new RealSystemEvent(DEPARTURE,i));// Departure for each server
        }
        return eSet;
	}
	@Override
    public String label() {
        String stg = "";
        switch (type) {
        		case (ARRIVAL) :
                        stg += "Arrival ";
                        break;
                case (PHASECHANGE) :
                        stg += "Phase change in server"+(server+1);
                        break;
                case (DEPARTURE) :
                        stg += "Departure from server "+(server+1);
                        break;
        }
		return stg;
	}

}

/**
 * @author Leonardo Lozano and Laura Vielma
 *
 */
public class PhaseSystem extends SimpleMarkovProcess<RealSystemState, RealSystemEvent>{

	double lambda;//arrival rate
     int Max;//Max System Capacity  
     ContPhaseVar[] Distributions; // Service phase distributions
     
    static final int ARRIVAL = 1;
	static final int DEPARTURE = 2;
	static final int PHASECHANGE = 3;
     
     
     /**
     * @param Distributions
     * @param lambda
     * @param Max
     */
    public PhaseSystem (ContPhaseVar[] Distributions, double lambda,int Max) {
    	  super(new RealSystemState(), //
	                RealSystemEvent.getAllEvents());
            
             this.lambda = lambda;
             this.Max=Max;
             this.Distributions=Distributions;
     }
     /**
      * Determines the active events. 
      */
     @Override
     public boolean active(RealSystemState i, RealSystemEvent e) {
             boolean result = false;
             switch (e.type) {
             
             case ARRIVAL:// arrivals can happen if there is room in the system
            	if(i.getQ1Size()+i.getQ2Size()<Max){
            		result =true;
 	            }
            	break;
 	       
 	        case DEPARTURE:  // A departure can occur if the server is active
 	        	if(i.getProp(e.server)!= 0){
 	        		result=true;
 	        	}
 	        	break;
 	        	
 	        case PHASECHANGE:// can occur if a server is busy 
 	        	if(i.getProp(e.server)!=0){
 	        		result = true;
 	        		}
 	            break;
 	         }
             return result;
}
 
//Determines the possible destinations
    @Override
     public States<RealSystemState> dests(RealSystemState i, RealSystemEvent e) {
             
    	     StatesSet<RealSystemState> set = new StatesSet<RealSystemState>();//New empty space set
             int[] props = new int[2];
             double alpha1[]=Distributions[0].getVectorArray();//Innitial probabilities vector for the first server 
             double alpha2[]=Distributions[1].getVectorArray();//Innitial probabilities vector for the second server 
             int Q1 = i.getQ1Size();
             int Q2 = i.getQ2Size();
             for(int k=0;k<2;k++){
             props[k] = i.getProp(k);}        //copy current values
             switch (e.type) {
             case ARRIVAL:
    	      	 if (i.getProp(0)== 0){//if the first server is idle, the client enters service in any possible phase
        		for(int k =1; k<=Distributions[0].getNumPhases();k++){
              		if(alpha1[k-1]>0){
        			props[0]=k;
        			set.add(new RealSystemState(props,Q1,Q2));}}
        		}
        	 	                	
              	 if(i.getProp(0)!= 0){// if the first server is busy the client must stay in queue
            	 Q1++;
            	 set.add(new RealSystemState(props,Q1,Q2));}
             break;
            
        case DEPARTURE:
          if(e.server==0){
        	  if (i.getQ1Size()==0 && i.getProp(1)==0){// if there is no queue 1 and the second server is idle, the first server 
        		  									  //state is set to idle and the client enters service with the second server in any possible phase
        		props[e.server]=0;
               	for(int k =1; k<=Distributions[1].getNumPhases();k++){
               		if(alpha2[k-1]>0){
               		props[1]=k;
        			set.add(new RealSystemState(props,Q1,Q2));}}
        	  }
        	 
        	  if (i.getQ1Size()==0 && i.getProp(1)!=0){// if there is no queue 1 and the second server is busy, the first server state
                Q2++;  								  // is set to idle, but the client must enter the second queue
          		props[e.server]=0;
             	set.add(new RealSystemState(props,Q1,Q2));}
        	
        	  if(i.getQ1Size()!=0 && i.getProp(1)==0   ){//if there is queue 1 and the second server is idle then a new client enters the first
        		  									     //server in any possible phase, and the departing client enters the second server in any possible phase
            	Q1--;
            	for(int k =1; k<=Distributions[0].getNumPhases();k++){
            		for(int u =1; u<=Distributions[1].getNumPhases();u++){
            		 if(alpha1[k-1]>0 && alpha2[u-1]>0){	
            			props[0]=k;
            			props[1]=u;
            			set.add(new RealSystemState(props,Q1,Q2));}}}
            	  }
        	  
        	  if(i.getQ1Size()!=0 && i.getProp(1)!=0 ){//if there is queue 1 and the second server is busy then a new client enters the first
				     									//server in any possible phase, and the departing client enters queue 2
            	Q1--;
            	Q2++;
            	for(int k =1; k<=Distributions[0].getNumPhases();k++){
            		if(alpha1[k-1]>0){
            		props[0]=k;
        			set.add(new RealSystemState(props,Q1,Q2));}}
                }
            
          }
         else{// when departure i from second server, the client leaves the system
        	 if (i.getQ2Size()==0){// if there is no queue 2, the second server state is set to idle
                 
             	props[e.server]=0;
             	set.add(new RealSystemState(props,Q1,Q2));}
        	
        	 else{// if there is queue 2, a waiting client enters second server in any possible phase 
        		 Q2--;
        		for(int k =1; k<=Distributions[1].getNumPhases();k++){
        			if(alpha2[k-1]>0){
        				props[1]=k;
            			set.add(new RealSystemState(props,Q1,Q2));}}
        		 }
         }     	    
            break;
            
        case PHASECHANGE:
        	int aux=props[e.server]; //Saves actual phase 	
        	for(int k =1; k<=Distributions[e.server].getNumPhases();k++){//phase can change to any other phase except the actual phase
          		if(k!=aux){
        		props[e.server]=k;
    			set.add(new RealSystemState(props,Q1,Q2));}}
           break;  }
             
             return set;
     }

    // Determines the events rates
     @Override
     public double rate(RealSystemState i, RealSystemState j, RealSystemEvent e) {
             double rate = 0;
             Matrix A= Distributions[e.server].getMatrix();// Phase change rate matrix
             switch (e.type) {
             case ARRIVAL:
            	double prob[]=Distributions[e.server].getVectorArray();//Innitial probabilities 
            	if(i.getProp(0)==0){
            	rate = lambda*prob[j.getProp(0)-1];} 
            	else{
            		rate=lambda;
            	}
 	            break;
 	            
 	        case DEPARTURE:
 	        	
 	        	double a[] = Distributions[e.server].getMat0Array();//Exit vector
 	        	rate=a[i.getProp(e.server)-1];// Exit vector evaluated in actual phase
 	        	break;
 	       
 	        case PHASECHANGE:
 	        	
	            rate= A.get(i.getProp(e.server)-1, j.getProp(e.server)-1);//Rate matrix evaluated in actual phase and the next phase
 	        	
 	            break;                                 
             }
             return rate;
     }
  
     @Override
     public String description() {
    	 return "Real system with phase distributions";
     }
	
	/**
	 * @param s
	 */
	public static void main(String s[]) {
		double[] data1 ={0.28,	0.19,	0.19,	0.31,	0.26,	0.24,	0.34,	0.18,	0.12,	0.26,	0.13,	0.27,	0.35,	0.2,	0.27,	0.55,	0.19,	0.19,	0.08,	1.03,	0.22,	0.14,	0.38,	0.25,	0.15,	0.17,	0.24,	0.11,	0.11,	0.17,	0.13,	0.3,	0.15,	0.2, 2.09,	0.21,	0.28,	0.15,	0.22,	0.21,	0.18,	0.38,	0.14,	0.21,	0.07,	0.11,	0.35,	0.07,	0.25,	0.3, 0.21,	0.18,	0.16,	0.12,	0.18,	0.28,	0.14,	0.25,	0.33,	0.15,	0.45,	0.13,	0.37,	0.27,	0.26,	0.17,	0.15,	0.43,	0.22,	0.33,	0.41,	0.31,	0.23,	0.15,	0.31};
		double[] data2={0.29,	0.23,	0.23,	0.04,	0.04,	0.19,	0.08,	0.03,	0.03,	0.04,	0.1,	0.08,	0.08,	0.11,	0.21,	0.08,	0.14,	1.03,	0.08,	0.52,	0.59,	0.1,	0.1,	0.2, 0.22,	0.29,	0.11,	0.32,	0.1,	0.08,	0.05,	0.06,	0.19,	0.06,	0.11,	0.07,	0.09,	0.06,	0.21,	0.11,	0.04,	0.18,	0.16,	0.09, 	0.06,	0.11,	0.55,	0.08,	0.12,	0.04,	0.12,	0.1,	0.14,	0.11,	0.09,	0.07,	0.13,	0.11,	0.57,	0.05,	0.07,	0.08,	0.25,	0.03};
        MomentsECCompleteFit fitter1 = new MomentsECCompleteFit(data1);
        MomentsECCompleteFit fitter2 = new MomentsECCompleteFit(data2);
        ContPhaseVar Servicios1 = fitter1.fit();
        ContPhaseVar Servicios2 = fitter2.fit();
        
        ContPhaseVar[] Distributions = {Servicios1, Servicios2};
        double lambda =2.63;
        int Max=11;
        	
        PhaseSystem TheModel = new PhaseSystem(Distributions, lambda,Max);
        TheModel.printAll();
        TheModel.showGUI();
	}
}