/*
 * Created on 3/15/2013
 * Test commit.
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
import jmarkov.basic.DecisionRule;
import jmarkov.basic.Events;
import jmarkov.basic.EventsSet;
import jmarkov.basic.Policy;
import jmarkov.basic.Solution;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDPEv;
import examples.jmdp.TandemEvent.TEvent;
import jmarkov.jmdp.solvers.PolicyIterationSolver;
import jmarkov.jmdp.solvers.PolicyIterationSolverAvg;
import jmarkov.jmdp.solvers.RelativeValueIterationSolver;
import jmarkov.jmdp.solvers.ValueIterationSolver;
 

/**
 * The present example describes a two station tandem queuing system, with capacity b1 in the 
 * first queue and b2, in the second. We must decide whether on an optimal admission control policy
 * based on the cost of rejecting an arrival at the first station vs. admitting it and then 
 * a loss ocurring at the second station. We solve for the long-run average cost in the infinite horizon.
 * 
 * @author Daniel Silva
 * 
 */
public class AccessControl extends
        DTMDPEv<TandemQueues, Admit, TandemEvent> {

    double lambda, mu1, mu2, c1, c2;
    int b1,b2;

    /**
     * Builds a AccessControl
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
     */

    public AccessControl(States<TandemQueues> init, double c1, double c2,double lambda, double mu1, double mu2, 
            int b1, int b2) {
    	super(init);
        this.c1 = c1;
        this.c2 = c2;
        this.lambda = lambda;
        this.mu1 = mu1;
        this.mu2 = mu2;
        this.b1 = b1;
        this.b2 = b2;
    }

    @Override
    public Events<TandemEvent> activeEvents(TandemQueues i, Admit a) {
    	EventsSet<TandemEvent> eventSet = new EventsSet<TandemEvent>();
        eventSet.add(new TandemEvent(TEvent.Arrival));
        eventSet.add(new TandemEvent(TEvent.Service1));
        eventSet.add(new TandemEvent(TEvent.Service2));
        return eventSet;
    }


    @Override
    public States<TandemQueues> reachable(TandemQueues i, Admit a,
            TandemEvent e) {
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
    public Actions<Admit> feasibleActions(TandemQueues i) {
        ActionsSet<Admit> set = new ActionsSet<Admit>();
        if (i.getQ1() < b1) {// is admit possible?
           set.add(new Admit(1));
        }
        if (i.getQ1() +i.getQ2() >= b2){
        	set.add(new Admit(0));
        }
        return set;
    }

    @Override
    public double immediateCost(TandemQueues i, Admit a,
            TandemEvent e) {
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
            TandemEvent e) {
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
    public double prob(TandemQueues i, TandemEvent e) {
        double d = lambda + mu1 +mu2;
    	switch (e.get()) {
        case Arrival: {
        	return lambda/d;
        }
        case Service1: {
            return mu1/d;
        }
        case Service2: {
            return mu2/d;
        }
         default: {
        	return 0;
        }
        }
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
    		int b1 = 4, b2 = 4;
            double c1 = 461.995, c2 = 841.059, lambda = 91.524, mu1 = 98.421, mu2 = 98.421;
            int[] initState = { 0, 0 };
            States<TandemQueues> init = new StatesSet<TandemQueues>(new TandemQueues(
                    initState));

            AccessControl AC = new AccessControl(init, c1,c2, lambda, mu1, mu2, b1, b2);

            PolicyIterationSolverAvg<TandemQueues, Admit> theSolver = 
                    	new PolicyIterationSolverAvg<TandemQueues, Admit>(AC);
                    
            AC.setDebugLevel(4);
            
            theSolver.setPrintBias(true);
            theSolver.setPrintGain(true);
            theSolver.setPrintValueFunction(true);
            AC.setSolver(theSolver);
            
            AC.solve();
            AC.printSolution();

    	}
//    	else{
//    		int b1 = 3, b2 = 3;
//            int[] initState = { 0, 0 };
//            int[] s1 = { 2, 1 };
//            int[] s2 = { 1, 2 };
//            int[] s3 = { 2, 2 };
//            int[] s4 = { 2, 3 };
//            int[] s5 = { 1, 3 };
//            int[] s6 = { 0, 3 };
//            
//            TandemQueues S1 = new TandemQueues(s1);
//            TandemQueues S2 = new TandemQueues(s2);
//            TandemQueues S3 = new TandemQueues(s3);
//            TandemQueues S4 = new TandemQueues(s4);
//            TandemQueues S5 = new TandemQueues(s5);
//            TandemQueues S6 = new TandemQueues(s6);
//            
//            States<TandemQueues> init = new StatesSet<TandemQueues>(new TandemQueues(
//                    initState));
//	    	String stg = "Lam \t Mu1 \t Mu2 \t C1 \t C2 \t g \t (2,1) \t (1,2) \t (2,2) \t (2,3) \t (1,3) \t (0,3) \n";
//	    	Random randGen = new Random(123456);
//    		for (int i=1; i<=50002; i++){
//	    		double c1= randGen.nextDouble()*999+1;
//	    		double c2= randGen.nextDouble()*(999-c1)+c1;
//	    		double l= randGen.nextDouble()*99+1;
//	    		double m1= randGen.nextDouble()*99+1;
//	    		double m2= randGen.nextDouble()*99+1;
//	    		//if(Math.pow((m1/(m1+m2)),2) <= Math.pow((m1/(m1+m2)),3)*(1 + 2*(m2/(m1+m2)) + 2*Math.pow((m2/(m1+m2)),2))) {i--; continue;}
//	    		AccessControl AC = new AccessControl(init, c1,c2, l, m1, m2, b1, b2);
//
//	            PolicyIterationSolver2<TandemQueues, Admit> theSolver = 
//	                    	new PolicyIterationSolver2<TandemQueues, Admit>(AC);
//	                    
//	            System.out.println(i);
//	            //if (i==2078||i==4758||i==4923||i==5779) continue;
//	            AC.setDebugLevel(0);
//	            AC.setSolver(theSolver);
//	            AC.solve();
//	            Policy<TandemQueues, Admit> pol = AC.getOptimalPolicy();
//	            DecisionRule<TandemQueues, Admit> rule = pol.getDecisionRule();
//	            DecimalFormat df1 = new DecimalFormat("#.###");
//	            double g= theSolver.getGain();
//	            Admit r1 = rule.getAction(S1);
//	            Admit r2 = rule.getAction(S2);
//	            Admit r3 = rule.getAction(S3);
//	            Admit r4 = rule.getAction(S4);
//	            Admit r5 = rule.getAction(S5);
//	            Admit r6 = rule.getAction(S6);
//	            stg+=df1.format(l)+"\t "+df1.format(m1)+"\t "+df1.format(m2)+"\t "+df1.format(c1)+"\t "+df1.format(c2)+"\t "+df1.format(g)+"\t "+r1+"\t "+r2+" \t "+r3+" \t "+r4+" \t "+r5+"\t "+r6+"\n"; 
//	    		
//    		}// end for
//    		
//	    	FileWriter outFile;
//			try {
//				outFile = new FileWriter("C:/Users/Owner/Documents/HomePC Files/Daniel Universidad/MIIND/Tesis/jMarkov2_0a/examples/jmarkov/AccessControlFiles/result2.txt");
//				PrintWriter out = new PrintWriter(outFile);
//		        StringWriter sw = new StringWriter();
//		        stg += sw.toString();
//				out.print(stg);
//				out.close();
//			} catch (IOException EE) {
//				EE.printStackTrace();
//			}
//
//    	}//end else

    	else{
    		int b1 = 4, b2 = 4;
            int[] initState = { 0, 0 };
            int[] s1 = { 0, 4 };
            int[] s2 = { 1, 3 };
            int[] s3 = { 1, 4 };
            int[] s4 = { 2, 2 };
            int[] s5 = { 2, 3 };
            int[] s6 = { 2, 4 };
            int[] s7 = { 3, 1 };
            int[] s8 = { 3, 2 };
            int[] s9 = { 3, 3 };
            int[] s10 = { 3, 4 };
            
            TandemQueues S1 = new TandemQueues(s1);
            TandemQueues S2 = new TandemQueues(s2);
            TandemQueues S3 = new TandemQueues(s3);
            TandemQueues S4 = new TandemQueues(s4);
            TandemQueues S5 = new TandemQueues(s5);
            TandemQueues S6 = new TandemQueues(s6);
            TandemQueues S7 = new TandemQueues(s7);
            TandemQueues S8 = new TandemQueues(s8);
            TandemQueues S9 = new TandemQueues(s9);
            TandemQueues S10 = new TandemQueues(s10);
            
            States<TandemQueues> init = new StatesSet<TandemQueues>(new TandemQueues(
                    initState));
	    	String stg = "Lam \t Mu1 \t Mu2 \t C1 \t C2 \t g \t (0,4) \t (1,3) \t (1,4) \t (2,2) \t (2,3) \t (2,4) \t (3,1) \t (3,2) \t (3,3) \t (3,4) \t C2 \t g \t (0,4) \t (1,3) \t (1,4) \t (2,2) \t (2,3) \t (2,4) \t (3,1) \t (3,2) \t (3,3) \t (3,4)\t C2 \t g \t (0,4) \t (1,3) \t (1,4) \t (2,2) \t (2,3) \t (2,4) \t (3,1) \t (3,2) \t (3,3) \t (3,4)\t C2 \t g \t (0,4) \t (1,3) \t (1,4) \t (2,2) \t (2,3) \t (2,4) \t (3,1) \t (3,2) \t (3,3) \t (3,4)\t C2 \t g \t (0,4) \t (1,3) \t (1,4) \t (2,2) \t (2,3) \t (2,4) \t (3,1) \t (3,2) \t (3,3) \t (3,4)\t C2 \t g \t (0,4) \t (1,3) \t (1,4) \t (2,2) \t (2,3) \t (2,4) \t (3,1) \t (3,2) \t (3,3) \t (3,4)\t C2 \t g \t (0,4) \t (1,3) \t (1,4) \t (2,2) \t (2,3) \t (2,4) \t (3,1) \t (3,2) \t (3,3) \t (3,4)\t C2 \t g \t (0,4) \t (1,3) \t (1,4) \t (2,2) \t (2,3) \t (2,4) \t (3,1) \t (3,2) \t (3,3) \t (3,4)\t C2 \t g \t (0,4) \t (1,3) \t (1,4) \t (2,2) \t (2,3) \t (2,4) \t (3,1) \t (3,2) \t (3,3) \t (3,4)\t C2 \t g \t (0,4) \t (1,3) \t (1,4) \t (2,2) \t (2,3) \t (2,4) \t (3,1) \t (3,2) \t (3,3) \t (3,4)\n";
	    	Random randGen = new Random(123456);
    		for (int i=1; i<=10000; i++){
	    		double c1= randGen.nextDouble()*999+1;
	    		double c2= c1;
	    		double l= randGen.nextDouble()*99+1;
	    		double m1= randGen.nextDouble()*99+1;
	    		double m2= randGen.nextDouble()*99+1;
	    		//if(Math.pow((m1/(m1+m2)),2) <= Math.pow((m1/(m1+m2)),3)*(1 + 2*(m2/(m1+m2)) + 2*Math.pow((m2/(m1+m2)),2))) {i--; continue;}
	            if(i%100==0)System.out.println(i+",");
	            else System.out.print(i+",");
	            if (i==73||i==2283||i==2375||i==2446||i==2574||i==3305||i==4021||i==4751||i==8819||i==4348||i==3412||i==6977) continue;	    		
	            for(int j=1; j <= 10; j++){
	    			c2+=randGen.nextDouble()*c1;
	    			AccessControl AC = new AccessControl(init, c1,c2, l, m1, m2, b1, b2);

		            PolicyIterationSolverAvg<TandemQueues, Admit> theSolver = 
		                    	new PolicyIterationSolverAvg<TandemQueues, Admit>(AC);
		                    

		            AC.setDebugLevel(0);
		            AC.setSolver(theSolver);
		            AC.solve();
		            Policy<TandemQueues, Admit> pol = AC.getOptimalPolicy();
		            DecisionRule<TandemQueues, Admit> rule = pol.getDecisionRule();
		            DecimalFormat df1 = new DecimalFormat("#.###");
		            double g= theSolver.getGain();
		            int r1 = (rule.getAction(S1).toString()=="Admit")?1:0;
		            int r2 = (rule.getAction(S2).toString()=="Admit")?1:0;
		            int r3 = (rule.getAction(S3).toString()=="Admit")?1:0;
		            int r4 = (rule.getAction(S4).toString()=="Admit")?1:0;
		            int r5 = (rule.getAction(S5).toString()=="Admit")?1:0;
		            int r6 = (rule.getAction(S6).toString()=="Admit")?1:0;
		            int r7 = (rule.getAction(S7).toString()=="Admit")?1:0;
		            int r8 = (rule.getAction(S8).toString()=="Admit")?1:0;
		            int r9 = (rule.getAction(S9).toString()=="Admit")?1:0;
		            int r10 = (rule.getAction(S10).toString()=="Admit")?1:0;
		            if (j==1) stg+=df1.format(l)+"\t "+df1.format(m1)+"\t "+df1.format(m2)+"\t "+df1.format(c1)+"\t "+df1.format(c2)+"\t "+df1.format(g)+"\t "+r1+"\t "+r2+" \t "+r3+" \t "+r4+" \t "+r5+"\t "+r6+"\t "+r7+"\t "+r8+"\t "+r9+"\t "+r10+"\t"; 
		            else stg+= df1.format(c2)+"\t "+df1.format(g)+"\t "+r1+"\t "+r2+" \t "+r3+" \t "+r4+" \t "+r5+"\t "+r6+"\t "+r7+"\t "+r8+"\t "+r9+"\t "+r10+"\t";
		            if (j==10) stg+="\n";
		            	
		            }
	            
    		}// end for
    		
	    	FileWriter outFile;
			try {
				outFile = new FileWriter("C:/Users/Owner/Documents/HomePC Files/Daniel Universidad/MIIND/Tesis/jMarkov2_0a/examples/jmarkov/AccessControlFiles/result2.txt");
				PrintWriter out = new PrintWriter(outFile);
		        StringWriter sw = new StringWriter();
		        stg += sw.toString();
				out.print(stg);
				out.close();
			} catch (IOException EE) {
				EE.printStackTrace();
			}

    	}//end else
    	
    }// end main

}
    

