package examples.jmarkov;

import jmarkov.GeomProcess;
import jmarkov.GeomRelState;
import jmarkov.MarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.StatesSet;

/**
 * This example represents a queueing system with Poisson arrivals at rate
 * $\lambda’$ if the system is empty, and a rate $\lambda$ otherwise. Costumers
 * require two exponential stages of service: the first at rate $\mu_1$ and the
 * second at rate $\mu_2$. A Matrix Geometric Example (Nelson chapter 9).
 * 
 * @author Julio Goez. Universidad de los Andes.
 */
public class QueueMH2k1 extends GeomProcess<QueueMH2k1State, QueueMH2k1Event> {
    // only for emtpy system.
    final int ARRIVAL1 = QueueMH2k1Event.ARRIVAL1;
    // When there are one or more
    // costumers in the system.
    final int ARRIVAL2 = QueueMH2k1Event.ARRIVAL2;
    // the first stage of service.
    final int SERVICE1 = QueueMH2k1Event.SERVICE1;
    // The Final stage of service.
    final int DEPARTURE = QueueMH2k1Event.DEPARTURE;
    private double lambda1, lambda2; // Arrivals rates
    private double mu1, mu2; // service rates

    /**
     * Constructs a M/H2(k)/1 queue whit arrival sates lambda1 and lamdda2 and
     * service rates mu1 and mu2.
     * 
     * @param lambda1
     *            Arrival rate when the system is empty.
     * @param lambda2
     *            Arival rate otherwise.
     * @param mu1
     *            Service rate for the first service stage.
     * @param mu2
     *            Service rate for the second service stage.
     */
    public QueueMH2k1(int lambda1, int lambda2, double mu1, double mu2) {
        super(new QueueMH2k1State(0), QueueMH2k1Event.getAllEvents());
        this.lambda1 = lambda1;
        this.lambda2 = lambda2;
        this.mu1 = mu1;
        this.mu2 = mu2;
    }

    /**
     * Used by GUI
     */
    public QueueMH2k1() {
        this(3, 2, 5, 4);
    }

    /**
     * @see jmarkov.GeomProcess#active(jmarkov.basic.State, Event)
     */
    @Override
    public boolean active(QueueMH2k1State i, int absLevel, QueueMH2k1Event e) {

        boolean result = false;

        switch (e.type) {
        case ARRIVAL1:
            result = (absLevel == 0);
            break;
        case ARRIVAL2:
            result = (absLevel != 0);
            break;
        case SERVICE1:
            result = (i.getServerStatus() == 1);
            break;
        case DEPARTURE:
            result = (i.getServerStatus() == 2);
            break;
        }

        return result;
    }

    /*
     * (non-Javadoc)
     * 
     * @see jmarkov.GeomProcess#dest(jmarkov.State, Event)
     */

    @Override
    public GeomRelState<QueueMH2k1State>[] dests(QueueMH2k1State i,
            int absLevel, QueueMH2k1Event e) {
        int newy = i.getServerStatus();
        int rLevel = 0;

        switch (e.type) {
        case ARRIVAL1:
            newy = 1;
            rLevel = 1;
            break;
        case ARRIVAL2:
            rLevel = 1;
            break;
        case SERVICE1:
            newy++;
            break;
        case DEPARTURE:
            if (absLevel == 1) {
                newy = 0;
            } else {
                rLevel = -1;
                newy--;
            }
            break;

        }
        QueueMH2k1State newSubState = new QueueMH2k1State(newy);
        GeomRelState<QueueMH2k1State> s;
        if (e.type == DEPARTURE && absLevel == 1) {
            s = new GeomRelState<QueueMH2k1State>(newSubState);
        } else {
            s = new GeomRelState<QueueMH2k1State>(newSubState, rLevel);
        }
        StatesSet<GeomRelState<QueueMH2k1State>> statesSet //
        = new StatesSet<GeomRelState<QueueMH2k1State>>(s);
        return statesSet.toStateArray();
    }

    /*
     * (non-Javadoc)
     * 
     * @see jmarkov.GeomProcess#rate(jmarkov.State, int)
     */
    @Override
    public double rate(QueueMH2k1State i, int iLevel, QueueMH2k1State j,
            int jLevel, QueueMH2k1Event e) {
        double result = 0;
        switch (e.type) {
        case ARRIVAL1:
            result = lambda1;
            break;
        case ARRIVAL2:
            result = lambda2;
            break;
        case SERVICE1:
            result = mu1;
            break;
        case DEPARTURE:
            result = mu2;
            break;
        }
        return result;
    }

    /*
     * (non-Javadoc)
     * 
     * @see jmarkov.SimpleMarkovProcess#description()
     */
    @Override
    public String description() {
        String stg = "System M/E(k)/1";
        stg += "\nLambda 1 = " + lambda1;
        stg += "\nLambda 2 = " + lambda2;
        stg += "\nMu 1     = " + mu1;
        stg += "\nMu 2     = " + mu2;
        return stg;
    }

    /**
     * Main method
     * @param a Not used
     */
    public static void main(String[] a) {
        QueueMH2k1 theQueue = new QueueMH2k1(7, 4, 9, 9);
        theQueue.setDebugLevel(0);
        theQueue.showGUI();
        // GeomState[] estados = theQueue.getGeomStates();
        // System.out.println(estados);
        // DenseMatrix A = theQueue.GetInitialSol();
        theQueue.printAll();

    }

}

/**
 * This class extends PropertiesState and define the properties of a system whit
 * one server and a queue with infinite capacity, it has one property, namely
 * the server status. The status is 0 if the server is idle.
 * 
 * @author Julio Goez. Universidad de los andes.
 */
class QueueMH2k1State extends PropertiesState {

    /**
     * We identify the states with is the current service stage of the costumer
     * in service (1,2) or 0 if idle.
     * 
     * @param status
     *            Current stage of service of costumer.
     */

    QueueMH2k1State(int status) {
        super(1);
        setProperty(0, status);
    }

    @Override
    public void computeMOPs(MarkovProcess mp) {
        setMOP(mp,"Utilization", (getServerStatus() != 0) ? 1 : 0);
        setMOP(mp,"Stage 1 Utilization", (getServerStatus() == 1) ? 1 : 0);
        setMOP(mp,"Stage 2 Utilization", (getServerStatus() == 2) ? 1 : 0);
    }

    /**
     * @return the server status
     */

    public int getServerStatus() {
        return getProperty(0);
    }

    /**
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        // TODO Complete
        return true;
    }

    @Override
    public String label() {
        return "S" + getProperty(0);
    }

    /*
     * (non-Javadoc)
     * 
     * @see jmarkov.State#description()
     */
    @Override
    public String description() {
        return "Server status = " + getServerStatus();
    }

}

/**
 * This class define the events.
 * 
 * @author Julio To change the template for this generated type comment go to
 *         Window>Preferences>Java>Code Generation>Code and Comments
 */
class QueueMH2k1Event extends Event {

    final static int ARRIVAL1 = 0; // only for emtpy system.
    final static int ARRIVAL2 = 1; // When there are one or more costumers in
    // the
    // system.
    final static int SERVICE1 = 2; // the first stage of service.
    final static int DEPARTURE = 3; // the second stage of service.
    int type;

    QueueMH2k1Event(int type) {
        this.type = type;
    }

    static EventsSet<QueueMH2k1Event> getAllEvents() {
        EventsSet<QueueMH2k1Event> eSet = new EventsSet<QueueMH2k1Event>();
        eSet.add(new QueueMH2k1Event(ARRIVAL1));
        eSet.add(new QueueMH2k1Event(ARRIVAL2));
        eSet.add(new QueueMH2k1Event(SERVICE1));
        eSet.add(new QueueMH2k1Event(DEPARTURE));
        return eSet;

    }

    @Override
    public String label() {
        String stg = "";
        switch (type) {
        case (ARRIVAL1):
            stg = "Arrival on empty";
            break;
        case (ARRIVAL2):
            stg = "Arrival";
            break;
        case (SERVICE1):
            stg = "First stage";
            break;
        case (DEPARTURE):
            stg = "Second stage";
            break;
        }
        return stg;
    }
}
