package examples.jmarkov;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.DecimalFormat;
import java.util.Random;

import examples.jmdp.InvLevel;
import examples.jmdp.Order;

import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.NotUnichainException;
import jmarkov.solvers.JamaSolver;
import jmarkov.solvers.MtjSolver;
import jmarkov.solvers.Solver;
import jmarkov.solvers.SteadyStateSolver;

/**
 * This class represents a network with 2 queues in tandem, where the second queue has a buffer 
 * of size B and the first has a unit buffer.
 * 
 * @author Daniel F. Silva (2011)
 * 
 */
public class AccessControl2Queues extends SimpleMarkovProcess<AccessControl2QueuesState, AccessControl2QueuesEvent> {
    private int B;
    private double Mu1;
    private double Mu2;
    private double L;
    private double c1;
    private double c2;
    private boolean pol;
    
    /**
     * General constructor.
     * 
     * @param B
     *            Buffer Size
     * @param Mu1
     *            First Server Rate.
     * @param Mu2
     *            Second Server Rate.
     * @param L
     *            Arrival Rate.
     * @param c1
     *            Cost of turning away.
     * @param c2
     *            Cost of lost customer.
     * @param pol
     *            Policy (true is greedy, false is prudent)
     *                 
     */
    
    public AccessControl2Queues(int B, double L, double Mu1, double Mu2, double c1, double c2, boolean pol,int[]length) {
    	super(new AccessControl2QueuesState(length), AccessControl2QueuesEvent.getAllEvents());
        this.B = B;
        this.Mu1 = Mu1;
        this.Mu2 = Mu2;
        this.c1 = c1;
        this.c2 = c2;
        this.L = L;
        this.pol = pol;
    }
    
    /**
     * Default Constructor used by GUI
     */
    
    public AccessControl2Queues() {
    	this(1, 1, 10, 1, 1, 20, true, new int[] {0,0} );
    }
    
    /**
     * Determine the active events.
     */

    @Override
    public boolean active(AccessControl2QueuesState i, AccessControl2QueuesEvent e) {
        boolean result = false;
        int Q1 = i.getFirst()+e.getFirst();
        int Q2 = i.getSecond()+e.getSecond();
        if(e.getType() == AccessControl2QueuesEvent.Type.LOSS)
        	result = (i.getFirst() == 1 && Q1 == 0 && i.getSecond() == B && Q2 == B);
        else
        	result = (0 <= Q1 && Q1 <= 1 && 0 <= Q2 && Q2 <= B);
        if (!pol && Q1 == 1 && Q2 == B)
        	result = false;
        return result;
    }

    @Override
    public States<AccessControl2QueuesState> dests(AccessControl2QueuesState i, AccessControl2QueuesEvent e) {
        int [] origin = i.getProperties();
        return new StatesSet<AccessControl2QueuesState>(i.Trans(origin, e));
    }

    /**
     * Returns the transition probability from State i to State j.
     */

    @Override
    public double rate(AccessControl2QueuesState i, AccessControl2QueuesState j, AccessControl2QueuesEvent e) {
    	int [] origin = i.getProperties();
        int [] dest = j.getProperties();
        switch (e.getType()){
        case ADMIT:
        	if(origin[0] == dest[0] -1) return L; 
        	break;
        case COMP_1:
        	if((origin[0] == dest[0] +1)&&(origin[1] == dest[1] -1)) return Mu1;
        	break;
        case COMP_2:
        	if(origin[1] == dest[1] +1) return Mu2;   
    		break;
    	default:
    		if(origin[0] == 1 && origin [1] == B && dest[0] == 0 && dest[1] == B) return Mu1;   
			break;
        }
        return 0;
    }

    
    @Override
    public String description() {
        String stg = "Access control to a 2 tandem queue system \n\nThere are " + B + " spaces in the buffer\n";
        stg += "The arrival rate is " + L + ".\n";
        stg += "The service rate at station 1 is " + Mu1 + ".\n";
        stg += "The service rate at station 2 is " + Mu2 + ".\n";
        StringWriter sw = new StringWriter();
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
         out.printf("    Buffer Size   = %6d ", this.B);
         out.printf("\n    Total Cost   = %6f ", getCost());
         out.println();
        return 0;
    } 
    
    /**
     * Main method.
     * 
     * @param a
     *            Not used.
     */
    public static void main(String[]a) {
    	boolean runOne = false;
    	boolean runRandom = true;
    	
    	if(runOne){
    		AccessControl2Queues ACQ = new AccessControl2Queues();
    		double cost = -1;
    		int tries=1;
    		while (cost == -1 && tries <= 10){
    			ACQ.reset();
    			cost = ACQ.getCost();
    			tries++;
    		}
    		ACQ.setDebugLevel(4);
    		ACQ.showGUI();
    		ACQ.printAll();
    		ACQ.printMOPs();
    	}
    	else{
	    	double greed=0;
	    	double prude=0;
	    	boolean greedy = true;
	    	String stg = "Buffer \t Lam \t Mu1 \t Mu2 \t C1 \t C2 \t GreedyCost \t PrudentCost \t Prudent_0_B \t Greedy_0_B \t Greedy_1_B \t cStar \n";

	    	if (runRandom){
		    	double pZeroBPrudent=0;
		    	double pZeroBGreedy=0;
		    	double pOneBGreedy=0; 
		    	Random randGen = new Random(123456);
		    	Random randTwo = new Random(654321);
	    		for (int i=1; i<=1000; i++){
		    		double c1= randGen.nextDouble()*999+1;
		    		double c2= randGen.nextDouble()*999+1;
		    		if(c1>=c2) {i--; continue;}
		    		System.out.print("\n i = ");
	    			System.out.print(i);
		    		double l= randGen.nextDouble()*99+1;
		    		double m1= randGen.nextDouble()*99+1;
		    		double m2= randGen.nextDouble()*99+1;
		    		int x = 1;
		    		if(m2 <m1 && m2< l) x=5;
		    		else if(m1 <m2 && m1< l) x=2;
		    		c1= randTwo.nextDouble()*999+1;
		    		c2= c1*(1+randTwo.nextDouble()/x);
		    		for(int b= 1; b <=500; b++){
	    				if(b>50)b+=49;
			    		AccessControl2Queues ACQ1 = new AccessControl2Queues(b, l, m1, m2, c1, c2, greedy, new int[] {0,0} );
			    		ACQ1.setSteadyStateSolver(new JamaSolver(ACQ1));
//			    		final MtjSolver.EnumSolver sol = MtjSolver.EnumSolver.valueOf("QMR");
//			    		MtjSolver Solver1 = new MtjSolver(ACQ1, sol, false);
//			    		ACQ1.setSteadyStateSolver(Solver1);
			    		ACQ1.setMaxStates(10010);
			    		int tries=1;
		        		greed=-1;
		        		while (greed == -1 && tries <= 10){
		        			ACQ1.reset();
		        			greed = ACQ1.getCost();
		        			tries++;
		        		}
		        		DecimalFormat df1 = new DecimalFormat("#.###");
		        		DecimalFormat df2 = new DecimalFormat("#.#######");
	        			pZeroBGreedy = ACQ1.getProb(0, b);
	        			pOneBGreedy = ACQ1.getProb(1, b);
		        		tries = 1;
		                AccessControl2Queues ACQ2 = new AccessControl2Queues(b, l, m1, m2, c1, c2, !greedy, new int[] {0,0} );
		                ACQ2.setSteadyStateSolver(new JamaSolver(ACQ2));
//			    		final MtjSolver.EnumSolver sol2 = MtjSolver.EnumSolver.valueOf("QMR");
//			    		MtjSolver Solver2 = new MtjSolver(ACQ2, sol2, false);
//			    		ACQ2.setSteadyStateSolver(Solver2);
			    		ACQ2.setMaxStates(10010);
		                prude=-1;
		                while (prude == -1 && tries <= 10){
		        			ACQ2.reset();
		        			prude = ACQ2.getCost();
		        			tries++;
		        		}
//		                double sum = 0;
//		                for (int k=0; k <=b-1;k++){
//		                	sum += (l*ACQ2.getProb(1, k)-l*ACQ1.getProb(1, k));
//		                }
//		                sum = sum+l*ACQ2.getProb(0, b)-l*ACQ1.getProb(1, b);
//		                double den = m1*ACQ1.getProb(1, b);
//			    		double cStar = sum/den;
	                    pZeroBPrudent = ACQ2.getProb(0, b);
		                tries = 0;
	                    stg+=b+"\t "+df1.format(l)+"\t "+df1.format(m1)+"\t "+df1.format(m2)+"\t "+df1.format(c1)+"\t "+df1.format(c2)+"\t "+df1.format(greed)+"\t "+df1.format(prude)+" \t "+df2.format(pZeroBPrudent)+" \t "+df2.format(pZeroBGreedy)+" \t "+df2.format(pOneBGreedy)+"\n"; //" \t "+df2.format(cStar)+"\n";                        		
		    		}
	    		}
	    	}
	    	else{
		    	double pZeroBPrudent=0;
		    	double pZeroBGreedy=0;
		    	double pOneBGreedy=0;
		    	int[] e = new int [] {1,2};
		    	
		    	for (int i=0; i < e.length; i++){
		    		int b=e[i];
		    		for (int j=0; j < e.length; j++){
		        		float l=e[j];
		        		for (int k=0; k < e.length; k++){
		            		float m1=e[k];
		            		for (int r=0; r < e.length; r++){
		                		float m2=e[r];
		                		for (int t=0; t < e.length; t++){
		                    		float c1=e[t];
		                    		for (int u=0; u < e.length; u++){
		                        		float c2= e[u];
		                            	AccessControl2Queues ACQ1 = new AccessControl2Queues(b, l, m1, m2, c1, c2, greedy, new int[] {0,0} );
		                        		int tries=1;
		                        		greed=-1;
		                        		while (greed == -1 && tries <= 10){
		                        			ACQ1.reset();
		                        			greed = ACQ1.getCost();
		                        			tries++;
		                        		}
	                        			pZeroBGreedy = ACQ1.getProb(0, b);
	                        			pOneBGreedy = ACQ1.getProb(1, b);
		                        		tries = 1;
		                                AccessControl2Queues ACQ2 = new AccessControl2Queues(b, l, m1, m2, c1, c2, !greedy, new int[] {0,0} );
		                                prude=-1;
		                                while (prude == -1 && tries <= 10){
		                        			ACQ2.reset();
		                        			prude = ACQ2.getCost();
		                        			tries++;
		                        		}
		                                pZeroBPrudent = ACQ2.getProb(0, b);
		                                tries = 0;
		                                stg+=b+"\t "+l+"\t "+m1+"\t "+m2+"\t "+c1+"\t "+c2+"\t "+greed+"\t "+prude+" \t "+pZeroBPrudent+" \t "+pZeroBGreedy+" \t "+pOneBGreedy+"\n";                        		
		                    		}
		                		}
		            		}
		        		}
		    		}
		    	}
	    	}    
	    	FileWriter outFile;
			try {
				outFile = new FileWriter("C:/Users/Owner/Documents/HomePC Files/Daniel Universidad/MIIND/Tesis/jMarkov2_0a/examples/jmarkov/AccessControlFiles/result.txt");
				PrintWriter out = new PrintWriter(outFile);
		        StringWriter sw = new StringWriter();
		        stg += sw.toString();
				out.print(stg);
				out.close();
			} catch (IOException EE) {
				EE.printStackTrace();
			}
    	}
    }

    
    /**
     * @return Returns the size of the buffer at the second station.
     */
    public int getB() {
        return B;
    }

    /**
     * @return Returns rate of the server at the first station.
     */
    public double getMu1() {
        return Mu1;
    }


    /**
     * @return Returns rate of the server at the second station.
     */
    public double getMu2() {
        return Mu2;
    }
    

    /**
     * @return Returns the job arrival rate.
     */
    public double getLambda() {
        return L;
    }
    
    public double getCost() {
    	try{
		    double [] probs =  getSteadyState();
		    AccessControl2QueuesState[] states = getStates().toStateArray();
		    double cost=0;
		    double totprob=0;
		    for (double i : probs) {
		    	totprob += i;
		    }
		    if (totprob < 0.99 || totprob > 1.01)
		    	return -1;
		    for(int i=0; i<states.length; i++){
			   	 if (states[i].getFirst()==1){
			   		 cost += c1*L*probs[i];
			   		 if (states[i].getSecond()==B)
			   			 cost += c2*Mu1*probs[i];
			   	 }
			   	 else if(!pol && states[i].getFirst()==0 && states[i].getSecond()==B)
			   		 cost += c1*L*probs[i];
		    }
		    return cost;
    	}
        catch (NotUnichainException e) {
            return 0;
        }
    }

    public double getProb(int k, int l) {
    	try{
		    double [] probs =  getSteadyState();
		    AccessControl2QueuesState[] states = getStates().toStateArray();
		    double totprob=0;
		    double prob=0;
		    for (double j : probs) {
		    	totprob += j;
		    }
		    if (totprob < 0.99 || totprob > 1.01)
		    	return -1;
		    for(int i=0; i<states.length; i++){
			   	 if (states[i].getFirst()==k && states[i].getSecond()==l)
			   			 prob = probs[i];
		    }
		    return prob;
    	}
        catch (NotUnichainException e) {
            return 0;
        }
    }


} // End of class

/**
 * A State in a 2 tandem queue system. Each number represents the number of jobs
 * at the current station at time t.
 * 
 */

class AccessControl2QueuesState extends PropertiesState {

    /**
     * Constructor for the AccessControl2QueuesState class
     * 
     * @param length
     *            Length of each queue
     */
    public AccessControl2QueuesState(int [] length) {	
        super(length);
        setProperty(0, length[0]);
        setProperty(1, length[1]);
    }

    /**
     * Returns the length of the first queue
     * 
     */
    public int getFirst() {
        return getProperty(0);
    }
    
    /**
	 *	Returns the length of the second queue
     * 
     */
    public int getSecond() {
        return getProperty(1);
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
    public AccessControl2QueuesState Trans(int [] i, AccessControl2QueuesEvent e) {
        int [] newlength;
        newlength = new int [2];
        newlength [0]= i[0] + e.getFirst();
        newlength [1]= i[1] + e.getSecond();
        return new AccessControl2QueuesState(newlength);
    }

    /**
     * Describes the States
     * 
     * @return a String description of the State.
     */
    @Override
    public String description() {
        int a = getFirst();
        int b = getSecond();
        String stg;
        stg="There are ";
        stg+=a;
        stg+=" jobs in Q1 and ";
        stg+=b;
        stg+=" jobs in Q2.";
        return stg;
    }
} // End of AccessControl2QueuesState class.

/**
 * Each Event characterizes one of the following three situations:
 * 
 * 		ADMIT:  A new job goes into Station 1.
 * 		Complete 1: A job is completed at station 1.
 * 		Complete 2: A job is completed at station 2.
 * 
 * @author Daniel Silva
 */
class AccessControl2QueuesEvent extends Event {

	/** Event types. */
    public static enum Type {
        /** Arrival is admitted in system */
        ADMIT,
    	/** A service is completed in station 1 */
        COMP_1,
        /** A service is completed in station 2. */
        COMP_2,
        /** A job is lost between stations 1 and 2. */
        LOSS,
    }

    private Type type; // event type
    public int first; // Change in queue 1
    public int second; // Change in queue 2
    
    /**
     * Creates an Event.
     * 
     * @param type
     */
    public AccessControl2QueuesEvent(Type type) {
    	this.type = type;
    	switch (type){
    	case ADMIT: 
    		this.first=+1;
    		this.second=0;
    		break;
    	case COMP_1: 
    		this.first=-1;
    		this.second=+1;
    		break;
    	case COMP_2: 
    		this.first=0;
    		this.second=-1;
    		break;
    	case LOSS: 
    		this.first=-1;
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
     * @return Change in length of first queue
     */
    
	public int getFirst() {
		return first;
	}

    /**
     * @return Change in length of second queue
     */
	
	public int getSecond() {
		return second;
	}

	/**
     * 
     * @return A set with all the events in the system.
     */
    public static EventsSet<AccessControl2QueuesEvent> getAllEvents() {
        EventsSet<AccessControl2QueuesEvent> eSet = new EventsSet<AccessControl2QueuesEvent>();
        eSet.add(new AccessControl2QueuesEvent(Type.ADMIT));
        eSet.add(new AccessControl2QueuesEvent(Type.COMP_1));
        eSet.add(new AccessControl2QueuesEvent(Type.COMP_2));
        eSet.add(new AccessControl2QueuesEvent(Type.LOSS));
        return eSet;
    }

    @Override
    public String label() {
        String stg = "";
        switch (type) {
        case ADMIT:
        	stg = "A new job enters the system.";
            break;
        case COMP_1:
            stg = "A job is completed at station 1.";
            break;
        case COMP_2:
            stg = "A job is completed at station 2.";
            break;
        case LOSS:
            stg = "A job is lost between stations 1 and 2.";
            break;        }
        return stg;
    }

} // End of AccessControl2QueuesEvent class.