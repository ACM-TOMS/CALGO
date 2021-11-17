package examples.jmarkov;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

import jmarkov.GeomProcess;
import jmarkov.GeomRelState;
import jmarkov.MarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.StatesSet;

// The main class
/**
 * This class implements a system with N servers with erlang services.
 */
public class ErlangQueue extends
        GeomProcess<ErlangQueueState, ErlangQueueEvent> {

    double lambda;// Arrival rate
    double mu;// Phase change rate
    int N;// number of servers
    int M;// Number of phases
    static final int ARRIVAL = 1;
    static final int DEPARTURE = 2;
    static final int PHASECHANGE = 3;


    /**
     * @param lambda
     * @param mu
     * @param N
     * @param M
     */
    public ErlangQueue(double lambda, double mu, int N, int M) {
        super(new ErlangQueueState(M, N), //
                ErlangQueueEvent.getAllEvents(N));
        this.lambda = lambda;
        this.mu = mu;
        this.N = N;
        this.M = M;

    }

    @Override
    public boolean active(ErlangQueueState i, int iLevel, ErlangQueueEvent e) {

        boolean result = false;

        switch (e.type) {
        case ARRIVAL:// arrivals can always happen

            result = true;
            break;

        case DEPARTURE: // A departure can occur if a server
            // is in the last phase

            if (i.getProp(e.server) == M) {
                result = true;
            }
            break;

        case PHASECHANGE:// can occur if a server is busy but
            // not in the last phase

            if (i.getProp(e.server) != 0 && i.getProp(e.server) != M) {
                result = true;
            }
            break;
        }
        return result;
    }

    /*
     * This fuction defines the relative destinations. We defined the
     * levels with respect to the queue length, so if there is 1
     * person in queue, that is the level 1, and so on.
     */
    @Override
    public GeomRelState<ErlangQueueState>[] dests(ErlangQueueState i,
            int absLevel, ErlangQueueEvent e) {

        StatesSet<GeomRelState<ErlangQueueState>> destStates = new StatesSet<GeomRelState<ErlangQueueState>>();
        int[] props = i.getProperties();
        int rLevel = 0;
        switch (e.type) {
        case ARRIVAL:

            if (i.getIdle() != -1) {// if there are idle
                // servers, the client can
                // pick any of them
                for (int k = 0; k < N; k++) {
                    if (props[k] == 0) {
                        props[k] = 1;
                        destStates.add(new GeomRelState<ErlangQueueState>(
                                new ErlangQueueState(props, N, M), rLevel));
                        props[k] = 0;
                    }
                }
            }

            if (i.getIdle() == -1) {// if there is no idle server,
                // the client must stay in
                // queue
                rLevel = +1;
                destStates.add(new GeomRelState<ErlangQueueState>(
                        new ErlangQueueState(props, N, M), rLevel));
            }

            break;

        case DEPARTURE:

            if (absLevel == 0) {// if there is no queue, the
                // server state is set to idle

                props[e.server] = 0;

                destStates.add(new GeomRelState<ErlangQueueState>(
                        new ErlangQueueState(props, N, M), rLevel));
            }

            else {// if there is queue, the queue length
                // decreases by one, and the client enters the
                // server in phase 1
                rLevel = -1;
                props[e.server] = 1;

                destStates.add(new GeomRelState<ErlangQueueState>(
                        new ErlangQueueState(props, N, M), rLevel));
            }

            break;

        case PHASECHANGE:
            rLevel = 0;

            props[e.server]++;// if phase change occurs, the
            // phase is increased by one
            destStates.add(new GeomRelState<ErlangQueueState>(
                    new ErlangQueueState(props, N, M), rLevel));
            break;
        }

        return destStates.toStateArray();

    }// end of dest

    /**
     * This method calculates the rate of transition from i to j when
     * the event e occurs.
     * @param i initial state.
     * @param j final state.
     * @param e event.
     */
    @Override
    public double rate(ErlangQueueState i, int iLevel, ErlangQueueState j,
            int jLevel, ErlangQueueEvent e) {

        double rate = -1;

        switch (e.type) {
        case ARRIVAL:

            rate = lambda * i.prob(); // a client chooses
            // between the servers
            // with the same
            // probability

            break;

        case DEPARTURE:
            rate = mu;// Departure rate is the rate of the last
            // phase change, so it is the same phase
            // change rate
            break;
        case PHASECHANGE:
            rate = mu;// phase change rate

            break;
        }
        return rate;
    }

    @Override
    public String description() {
        return "Erlang queue with " + N + " Servers and " + M + " Phases";
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        double lambda = 0;
        double mu = 0;
        int M = 0;
        int N = 0;

        BufferedReader rdr = // the model asks the user for the
        // values
        new BufferedReader(new InputStreamReader(System.in));
        try {
            System.out.println("Arrival Rate: ");
            lambda = Double.parseDouble(rdr.readLine());
            System.out.println("Service Mean time in minutes: ");
            mu = Double.parseDouble(rdr.readLine());
            System.out.println("Total Serves: ");
            N = Integer.parseInt(rdr.readLine());
            System.out.println("Total Phases: ");
            M = Integer.parseInt(rdr.readLine());

        } catch (IOException e) {
        }
        ;

        mu = ((M * 60) / mu);// the user provides the whole service
        // mean
        // time in minutes, so it is necessary to
        // tranform it to the rate of the phase
        // change,
        // since the distribution is an erlang all
        // phases have the same rate.

        
        ErlangQueue TheModel = new ErlangQueue(lambda, mu, N, M);
        TheModel.showGUI();
        TheModel.setDebugLevel(4);
        TheModel.generate();
        TheModel.printDenseMatrix(new PrintWriter(System.out));
        TheModel.printAll();

    }

}

/** 
 * Model of an Erlang Queue with N servers. The services follow an 
 * Erlang distribution with M phases. 
 * 
 * @author Leonardo Lozano,  Laura Vielma
 */
class ErlangQueueState extends PropertiesState {

    /**
     * We identify the states with the current phase of the n server
     * or 0 if idle and the Queue in the last position (1,..,n, Q).
     */
    int N;  
    int M; 

    /**
     * @param M Number of servers
     * @param N Number of phases
     */
    
    public ErlangQueueState(int M, int N) {

        this(new int[N], N, M);
    }

    // Constructor
    /**
     * @param status
     * @param N
     * @param M
     */
    public ErlangQueueState(int[] status, int N, int M) {

        super(status);
        this.N = N;
        this.M = M;
    }

    @Override
    public void computeMOPs(MarkovProcess mp) {

        setMOP(mp, "Queue Length", 1);

    }

    /**
     * @param i is position in the vector
     * @return an integer that shows the current state
     */
    public int getProp(int i) {
        return prop[i];
    }

    // Search for idle serves. Returns the position of the idle
    // server or -1 if all servers are busy

    /**
     * @return Search for idle serves. Returns the position of the idle server or -1 if all servers are busy
     */
    
    public int getIdle() {
        int aux = -1;
        for (int i = 0; i < N; i++) {
            if (prop[i] == 0) {
                aux = i;
            }
        }
        return aux;
    }

   

    /**
     * @return the probability of picking an idle server.
     */
    public double prob() {
        double aux = 0;
        double prob = 1;
        for (int i = 0; i < N; i++) {
            if (prop[i] == 0) {
                aux++;
            }
        }
        if (aux != 0) {
            prob = 1 / aux;
            return prob;
        }

        else {
            return prob;
        }
    }

    @Override
    public boolean isConsistent() {
        return true;
    }

    @Override
    public String description() {
        String stg = "";

        for (int k = 0; k < N; k++) {
            stg += " Server: ";
            stg += (k + 1) + " Phase: " + prop[k] + " ";
        }
        return stg;
    }

    @Override
    public String label() {
        String stg = "";

        for (int k = 0; k < N; k++) {
            stg += prop[k] + ((k != N-1) ? "," : "");
        }

        return stg;
    }

}

/*
 * There can be 3 kinds of events. An arrival, that can always occur,
 * a phase change that can occur if the server is in an inner phase,
 * or a departure that occurs when a server is in the last phase.
 */
class ErlangQueueEvent extends Event {

    static final int ARRIVAL = 1;
    static final int DEPARTURE = 2;
    static final int PHASECHANGE = 3;
    int type; // ARRIVAL, DEPARTURE, PHASECHANGE
    int server;// server number

    ErlangQueueEvent(int type, int server) {
        this.type = type;
        this.server = server;
    }

    static EventsSet<ErlangQueueEvent> getAllEvents(int N) {
        EventsSet<ErlangQueueEvent> eSet = new EventsSet<ErlangQueueEvent>();

        eSet.add(new ErlangQueueEvent(ARRIVAL, 0)); // arrivals
        // are
        // independent
        // of the
        // server
        // number, so
        // we let the
        // arrivals
        // happen only
        // for the
        // first
        // server

        for (int i = 0; i < N; i++) {
            eSet.add(new ErlangQueueEvent(PHASECHANGE, i));// phase
            // change
            // for
            // each
            // server
            eSet.add(new ErlangQueueEvent(DEPARTURE, i));
        }// Departure
        // for each
        // server

        return eSet;
    }

    /*
     * (non-Javadoc)
     * @see java.lang.Object#toString()
     */
    @Override
    public String label() {
        String stg = "";
        switch (type) {
        case (ARRIVAL):
            stg += "Arrival ";
            break;
        case (PHASECHANGE):
            stg += "Phase change in server" + (server + 1);
            break;
        case (DEPARTURE):
            stg += "Departure from server " + (server + 1);
            break;

        }
        return stg;
    }
}
