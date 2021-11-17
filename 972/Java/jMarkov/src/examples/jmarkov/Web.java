package examples.jmarkov;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import jmarkov.MarkovProcess;
import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.NotUnichainException;

/**
 * This is an example of a web server that can receive different M
 * types of items and can serve N simultaneous items with a maximum
 * queue capacity of MaxQueue. 
 * @author Laura Vielma and Leonardo Lozano
 */
public class Web extends SimpleMarkovProcess<WebState, WebEvent> {
    private double Gamma; // bytes per second that the server works
    // for each item served (no matter the
    // amount of items)
    private double[] mu; // mu[i] is the mean size in bytes for the
    // item i
    double[] lambda; // lambda[i] is the Arrival rate for an item i
    private double[] Prob; // Probability for an item in queue being
    // of a certain type
    private int M;// number of different items types
    private int N;// Max number of items that can be served
    // simultaneously
    private int MaxQueue;// Max number in queue
    private static final int ARRIVAL = WebEvent.ARRIVAL;
    private static final int DEPARTURE = WebEvent.DEPARTURE;

    /**
     * @param lambda
     * @param mu
     * @param M
     * @param Gamma
     * @param MaxQueue
     * @param N
     */
    public Web(double[] lambda, double[] mu, int M, double Gamma, int MaxQueue,
            int N) {
        super(new WebState(M), WebEvent.getAllEvents(M));
        this.M = M;
        this.lambda = lambda;
        this.mu = mu;
        this.Gamma = Gamma;
        this.N = N;
        this.MaxQueue = MaxQueue;
        int aux = 0;
        this.Prob = new double[M];
        for (int i = 0; i < M; i++) {
            aux += lambda[i]; // this is the sum of all rates
        }
        for (int i = 0; i < M; i++) {
            this.Prob[i] = lambda[i] / aux; // each probability is the
            // rate divided by the sum
            // of all rates
        }
    }

    /**
     * Determines the active events.
     */
    @Override
    public boolean active(WebState i, WebEvent e) {
        boolean result = false;
        switch (e.type) {
        case (ARRIVAL): // ARRIVAL occurs only if there is room in the
            // queue
            if (i.getQSize() < MaxQueue) {
                result = true;
            }
            break;

        case (DEPARTURE):// DEPARTURE occurs if there are items of
            // the type e.item in service
        {
            if (i.getNumber(e.item) > 0) {
                result = true;
            }

        }
        }

        return result;
    }

    /*
     * Determines the possible destinations
     */

    @Override
    public States<WebState> dests(WebState i, WebEvent e) {
        StatesSet<WebState> set = new StatesSet<WebState>();
        int[] NumItems = new int[M];
        for (int k = 0; k < M; k++)
            NumItems[k] = i.getNumber(k); // copy current values
        int Q = i.getQSize();
        switch (e.type) {
        case (ARRIVAL):
            if (i.getNumServe() == N) {
                Q++;
                set.add(new WebState(NumItems, Q, M));
            } // if there are N items being served, the new item
            // must stay in queue
            else {
                NumItems[e.item]++;
                set.add(new WebState(NumItems, Q, M));
            } // if there are less than N items being served, the
            // new item can be served immediatly

            break;

        case (DEPARTURE):
            if (i.getQSize() == 0) {
                NumItems[e.item]--;
                set.add(new WebState(NumItems, Q, M));
            } // If there is no queue, the item that departs leave
            // an empty space in the server
            else {
                Q--;
                NumItems[e.item]--;
                for (int k = 0; k < M; k++) {// if the item type i
                    // is the next to be
                    // served we must
                    // create a possible
                    // destination state
                    // for each item
                    NumItems[k]++; // considering the item k is the
                    // next to be served
                    set.add(new WebState(NumItems, Q, M)); // if
                    // there
                    // is
                    // queue,
                    // the
                    // next
                    // item
                    // that
                    // will be
                    // served
                    // can be
                    // any
                    // item
                    // depending
                    // on its
                    // probability.
                    NumItems[k]--;
                }// return to original value to consider another
                // type k
            }
        }

        return set;
    }

    @Override
    public double rate(WebState i, WebState j, WebEvent e) {
        double result = 0;
        int aux = e.type;// aux used to calculate the rate
        switch (e.type) {
        case (ARRIVAL):
            result = lambda[e.item]; // Arrival Rate is lambda
            break;
        case (DEPARTURE):
            if (i.getQSize() > 0) {
                for (int k = 0; k < M; k++) {
                    if (i.getNumber(k) != j.getNumber(k)) {
                        aux = k;
                    }
                }
                result = Gamma / mu[e.item] * i.getNumber(e.item) * Prob[aux];
            }// if there is queue, the rate for departure is
            // Gamma/Mu times the number of items type i being
            // served times the probability that the next item to
            // be served is of type i.
            else {
                result = Gamma / mu[e.item] * i.getNumber(e.item);
            }// if there is no queue, the rate is the max of X
            // exponencial services of type i.
            break;
        }
        return result;
    }

    /**
     * @param a
     */
    public static void main(String[] a) {
        double[] lambda = new double[0];
        double[] mu = new double[0];
        double Gamma = 0;
        ;
        int M = 0;
        int N = 0;
        int MaxQueue = 0;

        BufferedReader rdr = new BufferedReader(
                new InputStreamReader(System.in));
        try {
            System.out.println("Number of different items types:");
            M = Integer.parseInt(rdr.readLine());
            System.out
                    .println("Max number of items that can be served simultaneously:");
            N = Integer.parseInt(rdr.readLine());
            System.out.println("Bytes per second that the server works: ");
            Gamma = Double.parseDouble(rdr.readLine());
            System.out.println("Queue capacity: ");
            MaxQueue = Integer.parseInt(rdr.readLine());
            mu = new double[M];
            for (int k = 0; k < M; k++) {
                System.out.println("Item " + (k + 1) + " mean size in bytes:");
                mu[k] = Double.parseDouble(rdr.readLine());
            }
            lambda = new double[M];
            for (int k = 0; k < M; k++) {
                System.out.println("Item " + (k + 1) + " arrival rate:");
                lambda[k] = Double.parseDouble(rdr.readLine());
            }

        } catch (IOException e) {
        }
        ;

        Web TheModel = new Web(lambda, mu, M, Gamma, MaxQueue, N);
        TheModel.printAll();

        TheModel.showGUI();

    }

    /**
     * 
     */
    @Override
    public String description() {
        double lambda = 0;
        for (int i = 0; i < this.M; i++) {
            lambda += this.lambda[i];
        }
        String stg = "Web Server SYSTEM\n\n";
        stg += "One server queue with " + this.M
                + " different kinds of items\n";
        try {
            stg += "Avarge waiting time is "
                    + ((this.getMOPsAvg("Queue Length") / lambda) * 60)
                    + " Minutes";
        } catch (NotUnichainException e) {
            e.printStackTrace();
        }
        return stg;
    }

}

/**
 * Each state is defined by the number of items type i in service plus
   the number of items in queue. 
 * @author Laura Vielma and Leonardo Lozano
 */
class WebState extends PropertiesState {

    private int M; // Different items available.

    // Constructor for an empty state
    WebState(int M) {
        this(new int[M], 0, M);
    }

    // Generic Constructor
    WebState(int[] NumItems, int Qsize, int M) {
        super(M + 1);
        this.M = M;
        for (int i = 0; i < M; i++) {
            prop[i] = NumItems[i];
        }
        prop[M] = Qsize;
    }

    /**
     * Computes the MOPs
     * @see jmarkov.basic.State#computeMOPs(MarkovProcess)
     */
    @Override
    public void computeMOPs(MarkovProcess mp) {
        setMOP(mp, "Queue Length", getQSize());

    }

    /**
     * @param M is the position
     * @return the number of items for the M type
     */
    public int getNumber(int M) {
        return prop[M];
    }

    /**
     * @return the size of the queue
     */
    public int getQSize() {
        return prop[M];
    }

    @Override
    public boolean isConsistent() {

        return true;
    }

    /**
     * Returns a label for the state
     */
    @Override
    public String label() {
        String stg = "State";
        for (int k = 0; k < M; k++) {
            stg += " Item " + (k + 1) + ": " + prop[k];
        }
        stg += " Queue:" + prop[M];

        return stg;
    }

    /**
     * @return the total number of items being served in a moment
     */
    public int getNumServe() {
        int aux = 0;
        for (int i = 0; i < M; i++) {
            aux += prop[i];
        }
        return aux;
    }

}

/**
 * This class defines the events.
 */
class WebEvent extends Event {
    final static int DEPARTURE = 0;
    final static int ARRIVAL = 1;
    int type;// type of event
    int item;// type of item

    // Generic Constructor
    /**
     * @param type
     * @param item
     */
    public WebEvent(int type, int item) {
        super();
        this.type = type;
        this.item = item;
    }

    // This method gets all the possible events
    static EventsSet<WebEvent> getAllEvents(int items) {
        EventsSet<WebEvent> eSet = new EventsSet<WebEvent>();
        for (int i = 0; i < items; i++) {
            eSet.add(new WebEvent(ARRIVAL, i));
        }
        for (int i = 0; i < items; i++) {
            eSet.add(new WebEvent(DEPARTURE, i));
        }
        return eSet;
    }

    // Label for each event
    @Override
    public String label() {
        String stg = "";
        switch (type) {
        case (ARRIVAL):
            stg += "arrival item " + (this.item + 1);
            break;

        case (DEPARTURE):
            stg += "Departure item " + (this.item + 1);
            break;
        }
        return stg;
    }
}
