/*
 * Created on 7/08/2005
 */
package examples.jmdp;

import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Map.Entry;

import jmarkov.MarkovProcess;
import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.PropertiesAction;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.CTMDP;
import jmarkov.jmdp.solvers.ValueIterationSolver;

/**
 * The present example describes a Bank office with "maxServers" tellers. There
 * is a queue for clients for whom the bank wants to give the best service, and
 * another queue for users who have no commercial relation with the bank and are
 * not the service priority. At every transition, the bank can add or remove a
 * teller serving one of the queues and remove him or add him to the other
 * queue.
 * 
 * @author sarmiea
 * 
 */

public class Bank2Queues extends CTMDP<BankQueues, BankServers> {

    double clientMinuteCost, userMinuteCost, clientProbability, lambda, mu;
    int maxCapacity, maxServers;

    /**
     * Builds a Bank2Queues
     * 
     * @param initial
     *            Initial state in the system
     * @param clientMinuteCost
     *            The cost that is charged to clients with commercial relation
     *            with the bank per minute they stay in queue
     * @param userMinuteCost
     *            The cost that is charged to users with no commercial relation
     *            with the bank per minute they stay in queue
     * @param clientProbability
     *            Probability that the next client has a commercial relation
     *            with the bank
     * @param lambda
     *            Arrival rate to the bank
     * @param mu
     *            Service rate of tellers
     * @param maxCapacity
     *            Queues capacity
     * @param maxServers
     *            Maximum number of tellers available in the bank
     */
    public Bank2Queues(States<BankQueues> initial, double clientMinuteCost,
            double userMinuteCost, double clientProbability, double lambda,
            double mu, int maxCapacity, int maxServers) {
        super(initial);
        this.clientMinuteCost = clientMinuteCost;
        this.userMinuteCost = userMinuteCost;
        this.clientProbability = clientProbability;
        this.lambda = lambda;
        this.mu = mu;
        this.maxCapacity = maxCapacity;
        this.maxServers = maxServers;
    }

    // public enum Measure {
    // CLIENT_QUEUE_LENGTH {
    // @Override
    // double eval(BankQueues i) {
    // return i.getClients()-Math.min(i.getClients(),i.getClientServers()); } },
    // USER_QUEUE_LENGTH {
    // @Override
    // double eval(BankQueues i) {
    // return
    // i.getUsers()-Math.min(i.getUsers(),maxServers-i.getClientServers()); } },
    // EFFECTIVE_LAMBDA {
    // @Override
    // double eval(BankQueues i) {
    // return (i.getClients()+i.getUsers()<maxCapacity)?1:0; } },
    // CLIENT_SERVER_UTILIZATION {
    // @Override
    // double eval(BankQueues i) {
    // return (i.getClients()>i.getClientServers())?1:0; } },
    // USER_SERVER_UTILIZATION {
    // @Override
    // double eval(BankQueues i) {
    // return (i.getUsers()>maxServers-i.getClientServers())?1:0; } };
    //	
    // // Do arithmetic op represented by this constant
    // abstract double eval(BankQueues i);
    // private boolean calculated = false;
    // private double[] MOPs;
    // }

    // @Override
    /**
     * @author German Riano. Universidad de los Andes. (C) 2006
     * 
     */
    protected enum Measure {
        /** Average clients queue */
        CLIENT_QUEUE_LENGTH,
        /** Average users queue */
        USER_QUEUE_LENGTH,
        /** The rate of people how really enter to the bank */
        EFFECTIVE_LAMBDA,
        /** Average utilization of the clients server */
        CLIENT_SERVER_UTILIZATION,
        /** Average utilization of the users server */
        USER_SERVER_UTILIZATION;
    }

    /**
     * Computes the measure of performance m of queue i
     * 
     * @param i
     *            Specify which queue bank is wanted to measure
     * @param m
     *            Measure of performance required
     * @return The m measure of performance of queue i
     */
    public double computeMOPs(BankQueues i, Measure m) {
        switch (m) {
        case CLIENT_QUEUE_LENGTH:
            return i.getClients()
                    - Math.min(i.getClients(), i.getClientServers());
        case USER_QUEUE_LENGTH:
            return i.getUsers()
                    - Math.min(i.getUsers(), maxServers - i.getClientServers());
        case EFFECTIVE_LAMBDA:
            return (i.getClients() + i.getUsers() < maxCapacity) ? 1 : 0;
        case CLIENT_SERVER_UTILIZATION:
            return (i.getClients() > i.getClientServers()) ? 1 : 0;
        case USER_SERVER_UTILIZATION:
            return (i.getUsers() > maxServers - i.getClientServers()) ? 1 : 0;
        }
        return 0.0;
    }

    @Override
    public double lumpCost(BankQueues i, BankServers a) {
        return 0;
    }

    @Override
    public double continuousCost(BankQueues i, BankServers a) {
        int waitingClients = Math.max(i.getClients() - i.getClientServers(), 0);
        int waitingUsers = Math.max(i.getUsers()
                - (maxServers - i.getClientServers()), 0);
        return clientMinuteCost * waitingClients + userMinuteCost
                * waitingUsers;
    }

    @Override
    public States<BankQueues> reachable(BankQueues i, BankServers a) {
        StatesSet<BankQueues> set = new StatesSet<BankQueues>();
        int total = i.getClients() + i.getUsers();
        int servers = i.getClientServers();
        int ac = a.getClientServers();
        if ((ac <= 0) && (i.getClients() > 0) && (servers > 1)) { // client
            // departure
            int[] temp = { i.getClients() - 1, i.getUsers(), servers + ac };
            set.add(new BankQueues(temp));
        }
        if ((ac >= 0) && (i.getUsers() > 0) && (servers < maxServers - 1)) { // user
            // departure
            int[] temp = { i.getClients(), i.getUsers() - 1, servers + ac };
            set.add(new BankQueues(temp));
        }
        if ((ac == 0) && (total < maxCapacity) && (servers > 0)) { // client
            // arrival
            int[] temp = { i.getClients() + 1, i.getUsers(), servers };
            set.add(new BankQueues(temp));
        }
        if ((ac == 0) && (total < maxCapacity) && (servers < maxServers)) { // user
            // arrival
            int[] temp = { i.getClients(), i.getUsers() + 1, servers };
            set.add(new BankQueues(temp));
        }

        return set;
    }

    // @Override
    // public States<BankQueues> reached(BankQueues i, BankServers a) {
    // StatesCollection<BankQueues> set = new StatesCollection<BankQueues>();
    // int total = i.getClients()+i.getUsers();
    // int servers = i.getClientServers();
    // int ac = a.getClientServers();
    // if( (ac<=0) && (i.getClients()>0) && (servers>0) ){ //client departure
    // int[] temp = {i.getClients()-1, i.getUsers(), servers+ac};
    // set.add(new BankQueues(temp));
    // }if( (ac>=0) && (i.getUsers()>0) && (servers<maxServers) ){ //user
    // departure
    // int[] temp = {i.getClients(), i.getUsers()-1, servers+ac};
    // set.add(new BankQueues(temp));
    // }if( (ac==0) && (total<maxCapacity) && (servers>0)){ //client arrival
    // int[] temp = {i.getClients()+1, i.getUsers(), servers};
    // set.add(new BankQueues(temp));
    // }if( (ac==0) && (total<maxCapacity) && (servers<maxServers)){ //user
    // arrival
    // int[] temp = {i.getClients(), i.getUsers()+1, servers};
    // set.add(new BankQueues(temp));
    // }
    //
    // return set;
    // }

    @Override
    public double rate(BankQueues i, BankQueues j, BankServers a) {
        int clients = i.getClients();
        int users = i.getUsers();
        int clientServers = i.getClientServers();
        int jClients = j.getClients();
        int jUsers = j.getUsers();
        if (jClients - clients > 0) { // client arrival
            if (a.getClientServers() == 0) // no change in servers
                return lambda * clientProbability;
        } else if (jUsers - users > 0) // user arrival
            if (a.getClientServers() == 0) // no change in servers
                return lambda * (1 - clientProbability);
        if ((clients - jClients > 0) && (clientServers > 0)) { // client
            // departure ->
            // client server
            // decrease
            if (a.getClientServers() <= 0)
                return Math.min(clients, clientServers) * mu;
        } else if ((users - jUsers > 0) && (clientServers < maxServers)) // user
            // departure
            // ->
            // client
            // server
            // increase
            if (a.getClientServers() >= 0)
                return Math.min(users, maxServers - clientServers) * mu;
        return 0.0;
    }

    @Override
    public Actions<BankServers> feasibleActions(BankQueues i) {
        ActionsSet<BankServers> set = new ActionsSet<BankServers>();
        if (i.getClients() > 0) {// client departures are possible
            if (i.getClientServers() > 1) { // possible to decrease client
                // servers
                int[] temp = { -1 };
                set.add(new BankServers(temp));
            }
        }
        if (i.getUsers() > 0) {// user departures are possible
            if (i.getClientServers() < maxServers - 1) { // possible to
                // increase client
                // servers
                int[] temp = { +1 };
                set.add(new BankServers(temp));
            }
        }// do nothing
        int[] temp = { 0 };
        set.add(new BankServers(temp));

        return set;
    }

    /**
     * @return Average queue length
     * @throws SolverException
     */
    public double length() throws SolverException {
        double queue = 0;
        ValueFunction<BankQueues> probs = getSteadyStateProbabilities();
        Iterator<Entry<BankQueues, Double>> it = probs.iterator();
        Entry<BankQueues, Double> e;
        for (; it.hasNext();) {
            e = it.next();
            queue += e.getValue()
                    * computeMOPs(e.getKey(), Measure.CLIENT_QUEUE_LENGTH);
        }
        return queue;
    }

    /**
     * This method just tests the class.
     * 
     * @param args
     *            Not used
     * @throws SolverException
     */
    public static void main(String[] args) throws SolverException {
        int maxCapacity = 5, maxServers = 4;
        double clientMinuteCost = 10, userMinuteCost = 4, clientProbability = 0.4, lambda = maxServers * 2, mu = 2;
        int[] initState = { 0, 0, 1 };
        States<BankQueues> init = new StatesSet<BankQueues>(new BankQueues(
                initState));

        Bank2Queues prob = new Bank2Queues(init, clientMinuteCost,
                userMinuteCost, clientProbability, lambda, mu, maxCapacity,
                maxServers);

        ValueIterationSolver<BankQueues, BankServers> solv = new ValueIterationSolver<BankQueues, BankServers>(
                prob, 0.06);
        // RelativeValueIterationSolver<BankQueues,BankServers> solv = new
        // RelativeValueIterationSolver<
        // BankQueues,BankServers>(prob,1.1);

        solv.useGaussSeidel(false);
        // try{
        prob.solve();
        // prob.printSolution(new PrintWriter("salida.txt"));
        prob.printSolution();
        prob.getOptimalValueFunction().print(new PrintWriter(System.out, true));
        // } catch (FileNotFoundException ex) {
        // throw ex;
        // }
        // System.out.println(solv.getProcessTime()+" milliseconds");
        // PrintWriter pw = new PrintWriter(System.out,true);
        // ProbabilitySolver<BankQueues,BankServers> solf =new
        // ProbabilitySolver<BankQueues,BankServers>(prob);
        // solf.setJacobi(true);
        // prob.setProbabilitySolver(solf);
        // ValueFunction<BankQueues> probs = prob.getSteadyStateProbabilities();
        // probs.print(pw,"%-12S", "%10.5f");
        // System.out.println("\n"+prob.length() +" average client queue
        // length");

        // solf = new ProbabilitySolver<BankQueues,BankServers>(prob);
        // solf.setJacobi(false);
        // prob.setProbabilitySolver(solf);
        // probs = prob.getSteadyStateProbabilities();
        // probs.print(pw,"%-12S", "%10.5f");
        // System.out.println("\n"+prob.length() +" average client queue
        // length");
        //		
        // solf = new ProbabilitySolver<BankQueues,BankServers>(prob);
        // solf.setJacobi(true);
        // solf.setGaussSeidel(true);
        // prob.setProbabilitySolver(solf);
        // probs = prob.getSteadyStateProbabilities();
        // probs.print(pw,"%-12S", "%10.5f");
        // System.out.println("\n"+prob.length() +" average client queue
        // length");
    }

}

// class QueueLength extends MOPs<BankQueues>{
//	
// public QueueLength(){
// super(1);
// }
// @Override
// public double eval(BankQueues i){
// return i.getClients()-Math.min(i.getClients(),i.getClientServers());
// }
// }

class BankQueues extends PropertiesState {
    /**
     * @param people
     *            State of the bank
     */
    public BankQueues(int[] people) {
        super(people);
    }

    /**
     * @return Number of clients in the system
     */
    public int getClients() {
        return prop[0];
    }// in system

    /**
     * @return Number of users in the system
     */
    public int getUsers() {
        return prop[1];
    }// in system

    /**
     * @return Number of servers serving in clients queue
     */
    public int getClientServers() {
        return prop[2];
    }// servers for queue1

    /**
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        // TODO Auto-generated method stub
        return true;
    }

    /**
     * @see jmarkov.basic.State#computeMOPs(MarkovProcess)
     */
    @Override
    public void computeMOPs(MarkovProcess mp) {
        // TODO Auto-generated method stub
    }

}

class BankServers extends PropertiesAction {
    /**
     * @param clientServers
     *            Number of clients server available
     */
    public BankServers(int[] clientServers) {
        super(clientServers);
    }

    /**
     * @return Number of servers serving in clients queue
     */
    public int getClientServers() {
        return getProperty(0);
    }
}
