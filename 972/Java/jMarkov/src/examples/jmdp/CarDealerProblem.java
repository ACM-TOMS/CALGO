package examples.jmdp;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.DecimalFormat;

import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.Event;
import jmarkov.basic.Events;
import jmarkov.basic.EventsSet;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDPEv;
import jmarkov.jmdp.solvers.PolicyIterationSolver;
import jmarkov.jmdp.solvers.PolicyIterationSolverAvg;
import jmarkov.jmdp.solvers.RelativeValueIterationSolver;
import jmarkov.jmdp.solvers.ValueIterationSolver;
import Jama.Matrix;

/**
 * This example is intended to illustrate jMDP's capacity for modeling Markov Decision
 * Processes using events. It uses the classes InvLevel, Order and DemandEvent. 
 * This example models a car dealer's ordering it is a single item, periodic review, 
 * stochastic demand inventory problem. It is modeled like a long-run average cost, 
 * infinite horizon, Markov Decision Problem. Demand is assumed to be random according to 
 * a Poisson Process distribution with given rate per period. The system is a periodic
 * review problem in which an entity periodically checks the inventory level and
 * takes decisions according to the states he finds. There is a price of selling
 * each item and a cost for buying it. Besides, there is a holding cost incurred
 * when holding one item in stock from one period to another. There is also a
 * truck cost for ordering as units must be placed in trucks and only whole trucks can 
 * be sent, even if the 'last' truck is not full. We assume that backorders are not allowed.
 * The objective is to minimize the expected long-run average cost. We include a switch that 
 * allows modelling the problem with the objective of minimizing total discounted cost.
 * 
 * @author Daniel F. Silva
 */
public class CarDealerProblem extends DTMDPEv<InvLevel, Order, DemandEvent> {
    // Problem parameters:
    private int maxInv, truckSize;
    // Cost and demand parameters:
    private double truckCost, holdingCost, expDemand, price, cost, intRate;
    // Demand distribution information:
    private double[] demPMF, demCCDF, demandLoss1;
    private boolean isdisc = false;


    // Constructor
    /**
     * @param maxInv
     *            maximum physical capacity in inventory, warehouse size.
     * @param truckSize
     *            maximum items in each fixed cost order. Orders can be greater
     *            than this value, but will be charged more than one fixed cost.
     * @param truckCost
     *            fixed cost per order
     * @param price
     *            unit price
     * @param cost
     *            unit aquistion cost
     * @param holdingCost
     *            non-finantial holding cost (it does NOT include finantial cost)
     * @param discountFactor
     *            interest per period
     * @param expDemand
     *            demand mean
     * @param discounted
     *            Whether a discounted model (rather than average) is to be used.
     */

    @SuppressWarnings("unchecked")
    public CarDealerProblem(int maxInv, int truckSize, double truckCost, double price, double cost,
            double holdingCost, double discFactor, double expDemand, boolean discounted) {
        super(new StatesSet<InvLevel>(new InvLevel(0)));
        this.maxInv = maxInv;
        this.truckSize = truckSize;
        this.truckCost = truckCost;
        this.price = price;
        this.cost = cost;
        this.holdingCost = holdingCost;
        this.expDemand = expDemand;
        initializeProbabilities();
        this.isdisc = discounted;
        this.intRate = 1/discFactor-1;
        if (discounted)
            setSolver(new ValueIterationSolver(this, intRate));
        else
            setSolver(new RelativeValueIterationSolver(this));
    }
    
    /**
     * Creates vectors with the probability density function (pdf), the 
     * complimentary cumulative distribution function (CCDF) and the
     * first order Loss function for a Poisson distribution with parameter 
     * expDemand. A Poisson has infinite support, but the model only requires 
     * this information fro k=0,...,maxInv, so that is all we calculate.
     */
    private void initializeProbabilities() {
        demPMF = new double[maxInv + 1];
        demCCDF = new double[maxInv + 1];
        demandLoss1 = new double[maxInv + 1];
        double cdf, p = Math.exp(-expDemand);
        demPMF[0] = p;
        cdf = p; 
        demCCDF[0]=1;
        demandLoss1[0] = expDemand;
        int maxlevel = maxInv;
        for (int i = 1; i <= maxlevel; i++) {
            demCCDF[i] = 1-cdf; // P{demand >= i}
        	demPMF[i] = (p *= expDemand / i); // P{demand = i}
            cdf += p; // P{demand <= i}
            demandLoss1[i] = (expDemand - i) * (1 - cdf) + expDemand * p;  // = E[(D-i)^+]
        }
    }

    /**
     * Determines the feasible actions for state i. Namely, the actions
     * are all order amounts between 0 and maxInv-i.
     * @param i the inventory level at the end of the week
     */
    @Override
    public Actions<Order> feasibleActions(InvLevel i) {
        int max = maxInv - i.getLevel();
        Order[] vec = new Order[max + 1];
        for (int n = 0; n <= max; n++) {
            vec[n] = new Order(n);
        }
        return new ActionsSet<Order>(vec);
    }
    
    /**
     * Determines the active events for state i AFTER action a is taken. 
     * Namely, we have 2 types of events. Those where demand realizations  
     * less than the available, so the demand is between 0 and i+a-1, 
     * and those where demand is  i+a or greater, we define whether an event 
     * is of each type by a boolen variable, and the sales amount by an integer,
     * where sales are (i+a-demand)^+. 
     * @param i the inventory level before ordering
     * @param a the order amount
     */
    @Override
    public Events<DemandEvent> activeEvents ( InvLevel i , Order a) {
    	EventsSet<DemandEvent> eventSet = new EventsSet<DemandEvent>();
    	eventSet .add(new DemandEvent( i . getLevel () + a. getSize () , true ) ) ;
    	for ( int n = 0; n < i . getLevel () + a. getSize ( ) ; n++) {
    		eventSet .add(new DemandEvent(n, false ) ) ;
    	}
    	return eventSet ;
    }
    
    /**
     * Determines the states that are reacheble from state i, when
     * action a is taken and event e occurs. Namely, the only reachable state is
     * j= i + a - e.
     * @param i the inventory level before ordering
     * @param a the order amount
     * @param e the demand event
     */
    @Override
    public States<InvLevel> reachable ( InvLevel i , Order a, DemandEvent e) {
    	StatesSet<InvLevel> stSet = new StatesSet<InvLevel >();
    	if (e . getGreaterThan ( ) )
    		stSet .add(new InvLevel (0));
    	else
    		stSet .add(new InvLevel ( i . getLevel () + a. getSize () - e .getDemand ( ) ) ) ;
    	return stSet ;    
    }
    
    /**
     * Determines the probability that event e happens in state i.
     * Namely, this is the pmf of a Poission for e<i+a and the CCDF of
     * a Poisson for e=i+a. And zero for anything else.
     * @param i the inventory level
     * @param e the demand event
     */
	@Override
	public double prob(InvLevel i, DemandEvent e) {
    	if (e . getGreaterThan ( ) )
    	return demCCDF[e .getDemand ( ) ] ;
    	return demPMF[e .getDemand ( ) ] ; 
    }
    
	@Override
	public double prob(InvLevel i, InvLevel j, Order a, DemandEvent e) {
		if (j.getLevel()==i.getLevel()+a.getSize()-e.getDemand())
			return 1;
		else
			return 0;
	}

	
    /**
     * Determines the cost of ordering a units when there are i in inventory.
     * We refer the user to the user manual for a detailed explanation of 
     * these formulas.
     * @param i the inventory level
     * @param a the order amount
     */    
	@Override
    public double immediateCost(InvLevel i, Order a, DemandEvent e) {
    	int maxSale = i . getLevel () + a. getSize ( ) ;
    	double expectedSales = expDemand - demandLoss1[maxSale ];
    	double netProfit = price*expectedSales - orderCost (a. getSize()) - holdingCost * i . getLevel ( ) ;
    	return -netProfit ;
    	}
    double orderCost ( int x) {
    	return truckCost*Math. ceil ((double) x /truckSize ) + x * cost ;
    }
    
    /**
     * Returns the full transition probability matrix
     * @return 3-dimensional array of probabilities.
     */
    public final double[][][] getTheP(){
    	StatesSet<InvLevel> stts = getAllStates();
    	int numAct = getNumActions();
        double[][][] P = new double[numStates][numStates][numAct];
    	for(InvLevel i: stts){
    		for(Order a:feasibleActions(i)){
    			for(InvLevel j: reachable(i,a)){
    				P[i.getLevel()][j.getLevel()][a.getSize()]=prob(i,j,a);
    			}
    		}
    	}
        return P;
    }
    
    /**
     * Returns the full transition probability matrix
     * @return 3-dimensional array of probabilities.
     */
    public final double[][] getTheR(){
    	StatesSet<InvLevel> stts = getAllStates();
    	int numAct = getNumActions();
    	double[][] R = new double[numStates][numAct];
    	for (int k=0; k< numStates; k++)
    		for (int h=0; h< numAct; h++)
    			R[k][h]=100000000;
    	for(InvLevel i: stts){
    		for(Order a:feasibleActions(i)){
    			R[i.getLevel()][a.getSize()]=immediateCost(i,a);
    		}
    	}
        return R;
    }


    /**
     * Test Program.
     * 
     * @throws SolverException
     */
    public static void main( String a [ ] ) throws SolverException {
    	boolean runOnce=false;
    	if (runOnce){
        	int maxInventory = 10; int truckSize = 4; double lambda = 7; double price = 1100; double cost = 500;
        	double holdCost = 50; double truckCost = 800; double discFactor=0.9; 
	    	CarDealerProblem prob = new CarDealerProblem(maxInventory , truckSize , truckCost , 
	    			price , cost , holdCost, discFactor, lambda, true);
	    	prob.setDebugLevel(2);
	    	//double [][][] P= prob.getTheP();
	    	//double [][] R= prob.getTheR();
	    	PolicyIterationSolver<InvLevel,Order> polDiscSolver= new PolicyIterationSolver<InvLevel,Order>(prob, discFactor);
	    	prob.setSolver(polDiscSolver);
	    	prob.solve ( ) ;
	    	prob.printSolution ( ) ;
    	}
    	else{
        	int [] maxInventory = {10,50,100,150,200};//,1000,2000}; 
        	int [] truckSize = {4,20,40,60,80};//,400,800}; 
        	double [] lambda = {7,35,70,105,140};//,700,1400}; 
        	double [] truckCost = {800,4000,8000,12000,16000};//,80000,160000};
        	double price = 1100, cost = 500, holdCost = 50, discFactor=0.9, intRate=1/discFactor-1;
        	double solValueDisc, solRelValue, solPolicyDisc, solPolicyAvg;
        	double timeValueDisc, timeRelValue, timePolicyDisc, timePolicyAvg;
        	String stg = "i \t maxInventory \t truckSize \t lambda \t truckCost \t solValueDisc \t solPolicyDisc \t solRelValue \t solPolicyAvg \t timeValueDisc \t timePolicyDisc \t timeRelValue \t timePolicyAvg\n";
        	for (int i=0; i<=4; i++){
    	    	CarDealerProblem prob = new CarDealerProblem(maxInventory[i], truckSize[i], truckCost[i], 
    	    			price , cost , holdCost, discFactor, lambda[i], true  );
    	    	ValueIterationSolver<InvLevel,Order> discValSolver= new ValueIterationSolver<InvLevel,Order>(prob, intRate);
    	    	prob.setSolver(discValSolver);
    	    	prob.solve();
    	    	solValueDisc= -discValSolver.getValueFunction().get()[0];
    	    	timeValueDisc= discValSolver.getProcessTime();
    	    	PolicyIterationSolver<InvLevel,Order> polDiscSolver= new PolicyIterationSolver<InvLevel,Order>(prob, intRate);
    	    	prob.setSolver(polDiscSolver);
    	    	prob.solve();
    	    	solPolicyDisc= -polDiscSolver.getValueFunction().get()[0];
    	    	timePolicyDisc= polDiscSolver.getProcessTime();
    	    	prob = new CarDealerProblem(maxInventory[i], truckSize[i], truckCost[i], 
    	    			price , cost , holdCost, discFactor, lambda[i], false);
       	    	RelativeValueIterationSolver<InvLevel,Order> relValSolver= new RelativeValueIterationSolver<InvLevel,Order>(prob);
    	    	prob.setSolver(relValSolver);
    	    	prob.solve();
    	    	solRelValue= -relValSolver.getGain();
    	    	timeRelValue= relValSolver.getProcessTime();
    	    	PolicyIterationSolverAvg<InvLevel,Order> polAvgSolver= new PolicyIterationSolverAvg<InvLevel,Order>(prob);
    	    	prob.setSolver(polAvgSolver);
    	    	prob.solve();
    	    	solPolicyAvg= -polAvgSolver.getGain();
    	    	timePolicyAvg= polAvgSolver.getProcessTime();
    	    	DecimalFormat df1 = new DecimalFormat("##.###");
    	    	stg = stg + i +" \t" + maxInventory[i] +" \t"+ truckSize[i] + " \t" + lambda[i] +" \t"+ truckCost[i]+ " \t" +df1.format(solValueDisc) + " \t" + df1.format(solPolicyDisc) + " \t" + df1.format(solRelValue) + " \t" + df1.format(solPolicyAvg) + " \t" + df1.format(timeValueDisc/1000) +  " \t" + df1.format(timePolicyDisc/1000) + " \t" + df1.format(timeRelValue/1000) +" \t" + df1.format(timePolicyAvg/1000) + "\n";
    	    	System.out.print(stg);
    	    }
    	System.out.print(stg); 
		}    	
    }
}