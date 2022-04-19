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

/**
 * This class represents a system with 2 different exponential
 * servers with rates mu1 and mu2, respectively, and arrival rate
 * lambda.
 * @author German Riano. Universidad de los Andes.
 */

public class QueueMM2dN extends SimpleMarkovProcess<MM2dNState, QMM2dNEvent> {
    // Events
    final int ARRIVAL = 0;
    final int ARRIVAL1 = 1; // only for empty system
    final int ARRIVAL2 = 2; // only for empty system
    final int DEPARTURE1 = 3;
    final int DEPARTURE2 = 4;
    private double lambda;
    private double mu1, mu2, alpha;
    private int N;

    /**
     * Constructs a M/M/2d queue with arrival rate lambda and service
     * rates mu1 and mu 2.
     * @param lambda Arrival rate
     * @param mu1 Server 1 rate
     * @param mu2 Server 2 rate
     * @param alpha Probability of an arriving customer choosing
     *        server 1 (if both idle)
     * @param N Max number in the system
     */
    public QueueMM2dN(double lambda, double mu1, double mu2, double alpha, int N) {
        super((new MM2dNState(0, 0, 0)), //
                QMM2dNEvent.getAllEvents()   ); // num Events
        this.lambda = lambda;
        this.mu1 = mu1;
        this.mu2 = mu2;
        this.alpha = alpha;
        this.N = N;
    }

    /**
     * Returns an QueueMM2N object with arrival rate 4.0, service rate
     * of the first server 2.0, service rate of the second server 3.0,
     * probability of choose the first server 0.3 and capacity of 8
     * customers in the system. Used by GUI
     */
    public QueueMM2dN() {
        this(1.0, 2.0, 3.0, 0.3, 8);
    }

    /**
     * Determines the active events
     */
    public @Override boolean active(MM2dNState i, QMM2dNEvent e) {
        boolean result = false;
        switch (e.getType()) {
        case ARRIVAL:
            result = ((i.getQSize() < N - 2) && (!i.isEmpty()));
            break;
        case ARRIVAL1:
            result = i.isEmpty();
            break;
        case ARRIVAL2:
            result = i.isEmpty();
            break;
        case DEPARTURE1:
            result = (i.getStatus1() > 0);
            break;
        case DEPARTURE2:
            result = (i.getStatus2() > 0);
            break;
        }
        return result;
    }

    public @Override States<MM2dNState> dests(MM2dNState i, QMM2dNEvent e) {
        int newx = i.getStatus1();
        int newy = i.getStatus2();
        int newz = i.getQSize();

        switch (e.getType()) {
        case ARRIVAL:
            if (i.getStatus1() == 0) {
                newx = 1;
            } // serv 1 desocupado
            else if (i.getStatus2() == 0) {
                newy = 1;
            } // serv 2 desocupado
            else { // ambos ocupados
                newz = i.getQSize() + 1;
            }
            break;
        case ARRIVAL1:
            newx = 1;
            break;
        case ARRIVAL2:
            newy = 1;
            break;
        case DEPARTURE1:
            if (i.getQSize() != 0) {
                newx = 1;
                newz = i.getQSize() - 1;
            } else {
                newx = 0;
            }
            break;
        case DEPARTURE2:
            if (i.getQSize() != 0) {
                newy = 1;
                newz = i.getQSize() - 1;
            } else {
                newy = 0;
            }
            break;
        }
        return new StatesSet<MM2dNState>( new MM2dNState(newx, newy, newz));
    }

    public @Override double rate(MM2dNState i,MM2dNState j, QMM2dNEvent e) {
        double res = 0;
        switch (e.getType()) {
        case ARRIVAL:
            res = lambda;
            break;
        case ARRIVAL1:
            res = lambda * alpha;
            break;
        case ARRIVAL2:
            res = lambda * (1 - alpha);
            break;
        case DEPARTURE1:
            res = mu1;
            break;
        case DEPARTURE2:
            res = mu2;
            break;
        }
        return res;
    }

    @Override
    public String description() {
        return "M/M/2/N SYSTEM\nQueueing System with two servers, with rates "
                + mu1 + " and " + mu2 + ".\nArrivals are Poisson with rate "
                + lambda + ",\nand the maximum number in the system is " + N;

    }

    /**
     * This method just tests the class.
     * @param a Not used
     */
    public static void main(String[] a) {
        String stg;
        BufferedReader rdr = new BufferedReader(
                new InputStreamReader(System.in));
        try {
            System.out.println("Input rate ");
            stg = rdr.readLine();
            double lda = Double.parseDouble(stg);
            System.out.println("Service rate 1  ");
            stg = rdr.readLine();
            double mu1 = Double.parseDouble(stg);
            System.out.println("Service rate 2  ");
            stg = rdr.readLine();
            double mu2 = Double.parseDouble(stg);
            System.out.println("Provide alpha  ");
            stg = rdr.readLine();
            double alpha = Double.parseDouble(stg);
            System.out.println("Max in the system ");
            stg = rdr.readLine();
            int N = Integer.parseInt(stg);
            QueueMM2dN theQueue = new QueueMM2dN(lda, mu1, mu2, alpha, N);
            theQueue.showGUI();
            theQueue.printAll();
        } catch (IOException e) {
        }
        ;
    }

} // class end

/**
 * This is a particular case of propertiesState, whith three
 * properties, namely the server 1 and 2 status, plus the queue level.
 * @author German Riano. Universidad de los Andes.
 */

class MM2dNState extends PropertiesState {

    /**
     * We identify each State with the triplet (x,y,z), where x and y
     * are the status of the servers and z the number in queue (0,1,
     * ..,N-2).
     */

    MM2dNState(int x, int y, int z) {
        super(3); // Creates a PropertiesState with 3 properties.
        this.prop[0] = x;
        this.prop[1] = y;
        this.prop[2] = z;
    }

    @Override
    public void computeMOPs(MarkovProcess mp) {
        setMOP(mp, "Status Server 1", getStatus1());
        setMOP(mp, "Status Server 2", getStatus2());
        setMOP(mp, "Queue Length", getQSize());
        setMOP(mp, "Number in System", getStatus1() + 
        		getStatus2() + getQSize());
    }

    /**
     * Returns the status of the first Server
     * @return Status of the first Server
     */
    public int getStatus1() {
        return prop[0];
    }

    /**
     * Returns the status of the second Server
     * @return Status of the second Server
     */
    public int getStatus2() {
        return prop[1];
    }

    /**
     * Returns the size of the queue
     * @return Status of the size of the queue
     */
    public int getQSize() {
        return prop[2];
    }

    /*
     * isEmpty detects is the system is empty. It comes handy when
     * checking whether the events ARRIVAL1 and ARRIVAL2 are active.
     */
    boolean isEmpty() {
        return (getStatus1() + getStatus2() + getQSize() == 0);
    }

    /**
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        // TODO Complete
        return true;
    }

    /*
     * We implement label so that States are labeld 1, 1A, 1B, 2, 3,
     * ..., N-2
     */
    @Override
    public String label() {
        String stg = "0";
        if ((getStatus1() == 1) && (getStatus2() == 0))
            stg = "1A";
        if ((getStatus2() == 1) && (getStatus1() == 0))
            stg = "1B";
        if ((getStatus2() == 1) && (getStatus1() == 1))
            stg = "" + (2 + getQSize());
        return stg;
    }

    /*
     * This method gives a verbal description of the State.
     */
    @Override
    public String description() {
        String stg = "";
        stg += "Server 1 is " + ((getStatus1() == 1) ? "busy" : "idle");
        stg += ". Server 2 is " + ((getStatus2() == 1) ? "busy" : "idle");
        stg += ". There are " + getQSize() + " customers waiting in queue.";
        return stg;
    }

}

class QMM2dNEvent extends Event {
    /** Event types */
    public enum Type {
        /** An arrival */
        ARRIVAL,
        /** Arrival to server 1 (only for emtpy system) */
        ARRIVAL1,
        /** Arrival to server 2 (only for emtpy system) */
        ARRIVAL2,
        /** departure from server 1 */
        DEPARTURE1,
        /** departure from server 2 */
        DEPARTURE2;
    }

    private Type type;

    /**
     * @param type
     */
    public QMM2dNEvent(Type type) {
        super();
        this.type = type;
    }

    /**
     * @return Returns the type.
     */
    public final Type getType() {
        return type;
    }

    /**
     * @return the set of all events.
     */
    public static EventsSet<QMM2dNEvent> getAllEvents() {
        EventsSet<QMM2dNEvent> evSet = new EventsSet<QMM2dNEvent>();
        for (Type type : Type.values())
            evSet.add(new QMM2dNEvent(type));
        return evSet;
    }
}

// Now we define main the class

