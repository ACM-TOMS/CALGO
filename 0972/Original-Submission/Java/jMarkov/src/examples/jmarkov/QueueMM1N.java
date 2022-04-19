package examples.jmarkov;

import jmarkov.MarkovProcess;
import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import examples.jmarkov.QueueMM1NEvent.Type;

/**
 * This class implements an M/M/1/N queue. That is a single server
 * exponential queue with Poisson Arrivals. It extends SimpleMarkovProcess.
 * @author Germán Riaño. Universidad de los Andes.
 * @version 1.0a
 */

public class QueueMM1N extends
        SimpleMarkovProcess<QueueMM1NState, QueueMM1NEvent> {

    private double lambda;// arrival rate
    private double mu; // service rate
    private int N;// max in teh system

    /**
     * Constructs a M/M/1/N queue.
     * @param lambda Arrival rate
     * @param mu Service Rate
     * @param N Max number in the system.
     */
    public QueueMM1N(double lambda, double mu, int N) {
        super(//
                new QueueMM1NState(0),// 
                QueueMM1NEvent.getAllEvents()//
        );
        this.lambda = lambda;
        this.mu = mu;
        this.N = N;
    }

    /**
     * Returns an QueueMM1N object with arrival rate 4.0, service rate
     * 2.0 and capacity of 4 customers in the system. Used by GUI
     */
    public QueueMM1N() {
        this(4.0, 2.0, 4);
    }

    /**
     * @return Returns lambda.
     */
    public final double getLambda() {
        return lambda;
    }

    /**
     * @return Returns mu.
     */
    public final double getMu() {
        return mu;
    }

    /**
     * @return Returns max number of the System.
     */
    public final int getMax() {
        return N;
    }

    @Override
    public boolean active(QueueMM1NState i, QueueMM1NEvent e) {

        boolean result = false;
        switch (e.getType()) {
        case ARRIVAL:
            result = (i.getLevel() < N);
            break;
        case DEPARTURE:
            result = (i.getLevel() > 0);
            break;
        }
        return result;
    }

    public @Override States<QueueMM1NState> dests(QueueMM1NState i, QueueMM1NEvent e) {
        int newx= i.getLevel() ;

        switch (e.getType()) {
        case ARRIVAL:
            newx ++;
            break;
        case DEPARTURE:
            newx --;
            break;
        }
        return new StatesSet<QueueMM1NState>(new QueueMM1NState(newx));
    }

    /**
     * The rate is lambda ore mu depending on whether the event i s
     * arrival or departure.
     * @see SimpleMarkovProcess#rate(State, State, Event)
     */
    @Override
    public double rate(QueueMM1NState i, QueueMM1NState j, QueueMM1NEvent e) {
        return (e.getType() == Type.ARRIVAL) ? lambda : mu;
    }

    /**
     * This method just tests the class.
     * @param a Not used
     */
    public static void main(String[] a) {
        QueueMM1N theQueue = new QueueMM1N(4.0, 2.0, 4);
        theQueue.showGUI();
        theQueue.printAll();
    }

    @Override
    public String description() {
        return "SINGLE SERVER QUEUE (M/M/1/N)\n" + "\nArrival Rate    = "
                + lambda + "\nService Rate    = " + mu + "\nSystem capacity = "
                + N;
    }

} // class end

/**
 * This is a particular case of propertiesState, whith a singlle
 * property, nemly the level.
 */

class QueueMM1NState extends PropertiesState {
    /**
     * Constructs a new State with level x
     * @param x Amount of customers in the system.
     */
    QueueMM1NState(int x) {
        super(1);
        this.prop[0] = x;
    }

    /**
     * Compute the MOP´s
     * @see jmarkov.basic.State#computeMOPs(MarkovProcess)
     */
    @Override
    public void computeMOPs(MarkovProcess mp) {
        setMOP(mp, "Number in the system", getLevel());
        setMOP(mp, "Number in queue", getNumberQ());
        setMOP(mp, "Server Utilization", getServerStatus());
    }

    /**
     * Gets the queue level
     * @return Number of customers in the system.
     */
    public int getLevel() {
        return prop[0];
    }

    /**
     * The number of busy servers.
     * @return 1 if server is busy, 0 otherwise
     */
    public int getServerStatus() {
        return (prop[0] > 0) ? 1 : 0;
    }

    /**
     * The number in queue
     * @return the number in queue.
     */
    public int getNumberQ() {
        return getLevel() - getServerStatus();
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
     * Returns the system level (number in the system)
     */
    @Override
    public String label() {
        return "" + getLevel();
    }

    @Override
    public String description() {
        return "Number in system = " + getLevel() + "  (" + getNumberQ()
                + " in queue).";
    }
}

/**
 * @author German Riano. Universidad de los Andes. (C) 2006
 */
class QueueMM1NEvent extends Event {
    /**
     * @author German Riano. Universidad de los Andes. (C) 2006
     */
    public enum Type {
        /** Arrival event */
        ARRIVAL,
        /** Service completion event */
        DEPARTURE
    }

    private Type type;

    /**
     * @param type
     */
    public QueueMM1NEvent(Type type) {
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
     * @return all events
     */
    public static EventsSet<QueueMM1NEvent> getAllEvents() {
        EventsSet<QueueMM1NEvent> evtsSet = new EventsSet<QueueMM1NEvent>();
        evtsSet.add(new QueueMM1NEvent(Type.ARRIVAL));
        evtsSet.add(new QueueMM1NEvent(Type.DEPARTURE));
        return evtsSet;
    }
}
