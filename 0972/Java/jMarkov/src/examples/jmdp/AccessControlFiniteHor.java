/*
 * Created on 10/20/2013
 */
package examples.jmdp;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.DecimalFormat;
import java.util.Random;

import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.Events;
import jmarkov.basic.EventsSet;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.FiniteMDPEv;
import examples.jmdp.TandemEvent.TEvent;
import jmarkov.jmdp.solvers.FiniteSolver;
 

/**
 * The present example describes a two station tandem queuing system, with capacity b1 in the 
 * first queue and b2, in the second. We must decide whether on an optimal admission control policy
 * based on the cost of rejecting an arrival at the first station vs. admitting it and then 
 * a loss ocurring at the second station. This is a finite horizon model, we solve for discounted cost.
 * 
 * @author Daniel Silva
 * 
 */
public class AccessControlFiniteHor extends
        FiniteMDPEv<TandemQueues, Admit, TandemEvent> {

    double lambda, mu1, mu2, c1, c2;
    int b1,b2, N;

    /**
     * Builds a AccessControlFiniteHor
     * @param init Initial state
     * 
     * @param c1
     *            Cost of refusing to admit an arriving customer
     * @param c2
     *            Cost of losing a customer at the second station
     * @param b1
     *            Capacity of first station
     * @param b2
     *            Capacity of second station
     * @param lambda
     *            Arrival rate to the system
     * @param mu1
     *            Service rate of station 1
     * @param mu2
     *            Service rate at station 2
     * @param N
     *            Time Horizon
     */

    public AccessControlFiniteHor(States<TandemQueues> init, double c1, double c2,double lambda, double mu1, double mu2, 
            int b1, int b2, int N) {
    	super(init, N);
        this.c1 = c1;
        this.c2 = c2;
        this.lambda = lambda;
        this.mu1 = mu1;
        this.mu2 = mu2;
        this.b1 = b1;
        this.b2 = b2;
        this.N = N;
    }

    @Override
    public Events<TandemEvent> activeEvents(TandemQueues i, Admit a, int t) {
    	EventsSet<TandemEvent> eventSet = new EventsSet<TandemEvent>();
    	eventSet.add(new TandemEvent(TEvent.Arrival));
        eventSet.add(new TandemEvent(TEvent.Service1));
        eventSet.add(new TandemEvent(TEvent.Service2));
        return eventSet;
    }


    @Override
    public States<TandemQueues> reachable(TandemQueues i, Admit a,
            TandemEvent e, int t) {
        StatesSet<TandemQueues> set = new StatesSet<TandemQueues>();
        int ac = a.getAdm();
        switch (e.get()) {
        case Arrival: {
            int[] temp = { Math.min(i.getQ1() + ac,b1), i.getQ2()};
            set.add(new TandemQueues(temp));
            break;
        }
        case Service1: {
            int[] temp =  {(i.getQ1()==0) ? 0:i.getQ1()-1,(i.getQ1()==0) ? i.getQ2(): Math.min(i.getQ2()+1,b2)};
            set.add(new TandemQueues(temp));
            break;
        }
        case Service2: {
            int[] temp = {i.getQ1(),(i.getQ2()==0) ? 0: i.getQ2()-1};
            set.add(new TandemQueues(temp));
            break;
        }
        }
        return set;
    }

    @Override
    public Actions<Admit> feasibleActions(TandemQueues i, int t) {
        ActionsSet<Admit> set = new ActionsSet<Admit>();
        if (i.getQ1() < b1) {// is admit possible?
           set.add(new Admit(1));
        }
        set.add(new Admit(0));
        return set;
    }

    @Override 
    public double immediateCost(TandemQueues i, Admit a,
            TandemEvent e, int t) {
        int ac = a.getAdm();
        switch (e.get()) {
        case Arrival: {
            if (ac==0) return c1;
            else return 0;
        }
        case Service1: {
        	if (i.getQ2()==b2 && i.getQ1()>0) return c2;
            else return 0;
        }
        default: {
        	return 0;
        }
        }
    }

    
    @Override
    public double prob(TandemQueues i, TandemQueues j, Admit a,
            TandemEvent e, int t) {
        int ac = a.getAdm();
        switch (e.get()) {
        case Arrival: {
            if (ac==0 && i.getQ1()==j.getQ1() && i.getQ2()==j.getQ2()) return 1;
            else if (ac==1 && i.getQ1()+1==j.getQ1() && i.getQ2()==j.getQ2()) return 1;
            else return 0;
        }
        case Service1: {
            if (i.getQ1()==0 && i.getQ1()==j.getQ1() && i.getQ2()==j.getQ2()) return 1;
            else if (i.getQ1()>0 && i.getQ2()<b2 && i.getQ1()-1==j.getQ1() && i.getQ2()+1==j.getQ2()) return 1;
            else if (i.getQ1()>0 && i.getQ2()==b2 && i.getQ1()-1==j.getQ1() && i.getQ2()==j.getQ2()) return 1;
            else return 0;        }
         default: {
        	if (i.getQ2()==0 &&i.getQ1()==j.getQ1() && i.getQ2()==j.getQ2()) return 1;
        	else if (i.getQ2()>0  && i.getQ1()==j.getQ1() && i.getQ2()-1==j.getQ2()) return 1;
            else return 0;
        }
        }
    }    
    
    @Override
    public double prob(TandemQueues i, TandemEvent e, int t) {
        double d = lambda + mu1 +mu2;
    	switch (e.get()) {
        case Arrival: {
        	return lambda/d;
        }
        case Service1: {
            return mu1/d;        }
         default: {
        	return mu2/d;
        }
        }
    }   

    /**
     * @see jmarkov.jmdp.FiniteMDP#finalCost(jmarkov.basic.State)
     */
    @Override
    public double finalCost(TandemQueues i) {
        return 0;
    }

    
    /**
     * This method just tests the class.
     * 
     * @param args
     *            Not used
     * @throws SolverException
     */
    public static void main(String[] args) throws SolverException {
    	boolean runOne = false;

    	if(runOne){
    		int b1 = 3, b2 = 3, N=10;
	        double c1 = 10, c2 = 11, lambda = 8, mu1 = 5, mu2 = 3;
	        int[] initState = { 1, 1 };
	        States<TandemQueues> init = new StatesSet<TandemQueues>(new TandemQueues(
	                initState));
	
	        AccessControlFiniteHor AC = new AccessControlFiniteHor(init, c1,c2, lambda, mu1, mu2, b1, b2, N);
	
	        FiniteSolver<TandemQueues, Admit> theSolver = 
	        		new FiniteSolver<TandemQueues, Admit>(AC);
	        
	        AC.setDebugLevel(4);
	        AC.setSolver(theSolver);
	        AC.getSolver().setPrintValueFunction(true);
	        AC.solve();
	        AC.printSolution();
    	}
    	else{
			int b1 = 6, b2 = 6, N=10;
//	        double c1 = 183, c2 = 552, l = 90, m1 = 60, m2 = 2, x=.95;
	    	String stg = "Iter \t Lam \t Mu1 \t Mu2 \t C1 \t C2 \t C2+E \t MAX \n";
	    	DecimalFormat df1 = new DecimalFormat("#.###");

	    	Random randGen = new Random(123456);
			for (int i=1; i<=10000; i++){

	    		double c1= randGen.nextDouble()*999+1;
	    		double c2= randGen.nextDouble()*(999-c1)+c1;
	    		double l= randGen.nextDouble()*99+1;
	    		double m1= randGen.nextDouble()*99+1;
	    		double m2= randGen.nextDouble()*99+1;
	    		double x=randGen.nextDouble()*0.1;
	    		double [][] beta = new double [b1+1][b2+1];
	    		double [][] beta2 = new double [b1+1][b2+1];
	    		if(i%1000==0)System.out.println(i+", "); else if(i%100==0)System.out.print(i+", ");
   	    		for (int k=0 ; k <=b2 ; k++){		    	
   	    			for (int r=0 ; r <=b1 ; r++){
		    		int[] initState = { r, k };
		    		States<TandemQueues> init = new StatesSet<TandemQueues>(new TandemQueues(initState));
		    		AccessControlFiniteHor AC = new AccessControlFiniteHor(init, c1,c2, l, m1, m2, b1, b2, N);
		    		FiniteSolver<TandemQueues, Admit> theSolver = new FiniteSolver<TandemQueues, Admit>(AC);
		    		AccessControlFiniteHor AC2 = new AccessControlFiniteHor(init, c1,c2+c2*x, l, m1, m2, b1, b2, N);
		    		FiniteSolver<TandemQueues, Admit> theSolver2 = new FiniteSolver<TandemQueues, Admit>(AC2);
		                    
		            AC.setDebugLevel(0);
		            AC.setSolver(theSolver);
		            AC.solve();
		            AC2.setDebugLevel(0);
		            AC2.setSolver(theSolver2);
		            AC2.solve();
		            double[] g= theSolver.getOptimalValueFunction().get();
		            double[] g2= theSolver2.getOptimalValueFunction().get();
		            beta[r][k] = -g[0];
		            beta2[r][k] = -g2[0];
		            
//		            stg+=df1.format(l)+"\t "+df1.format(m1)+"\t "+df1.format(m2)+"\t "+df1.format(c1)+"\t "+df1.format(c2)+"\t "+df1.format(g[0])+"\t ("+r+","+k+") \n"; 
//		            stg+=df1.format(l)+"\t "+df1.format(m1)+"\t "+df1.format(m2)+"\t "+df1.format(c1)+"\t "+df1.format(c2+c2*x)+"\t "+df1.format(g2[0])+"\t ("+r+","+k+") \n";
		    		}
		    	}// end for
    		double mx=-1000000;
    		for (int k=0 ; k <=b2 ; k++){
    			for (int r=0 ; r <b1 ; r++){
    				if (k+r< b2)continue;
	    			double dif = (beta2[r+1][k]-beta2[r][k])-(beta[r+1][k]-beta[r][k]);
	    			if(dif>=mx) mx=dif;
	    		}
    		}
	        stg+=i+"\t"+df1.format(l)+"\t "+df1.format(m1)+"\t "+df1.format(m2)+"\t "+df1.format(c1)+"\t "+df1.format(c2)+"\t "+df1.format(c2+x*c2)+"\t "+mx+"\n";
			}// end for i
	    	FileWriter outFile;
	    	try {
		    	outFile = new FileWriter("C:/Users/Owner/Documents/HomePC Files/Daniel Universidad/MIIND/Tesis/jMarkov2_0a/examples/jmarkov/AccessControlFiles/resultFinite.txt");
		    	PrintWriter out = new PrintWriter(outFile);
			    StringWriter sw = new StringWriter();
			    stg += sw.toString();
				out.print(stg);
				out.close();
			} catch (IOException EE) {
				EE.printStackTrace();
			}

	}//end else    	
    	
    }//end main
}// end class
