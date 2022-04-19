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
 * lambda' if the system is empty, and a rate lambda otherwise. Costumers
 * require two exponential stages of service: the k-th stage has exponential
 * rate $\mu_k$. A Matrix Geometric Example (Nelson chapter 9).
 * 
 * @author Julio Goez. Universidad de los Andes. Test del repo
 */
public class QueueMEk1 extends GeomProcess<QueueMEk1State, QueueMEk1Event> {
    // only for emtpy system.
    final int ARRIVAL1 = QueueMEk1Event.ARRIVAL1;
    // When there are one or more
    // costumers in the system.
    final int ARRIVAL2 = QueueMEk1Event.ARRIVAL2;
    // the first stage of service.
    final int SERVICE1 = QueueMEk1Event.SERVICE1;
    // The Final stage of service.
    final int DEPARTURE = QueueMEk1Event.DEPARTURE;
    private double lambda1, lambda2; // Arrivals rates
    private double[] mu; // service rates

    /**
     * Constructs a M/H2(k)/1 queue whit arrival sates lambda1 and lamdda2 and
     * service rates mu1 and mu2.
     * 
     * @param lambda1
     *            Arrival rate when the system is empty.
     * @param lambda2
     *            Arival rate otherwise.
     * @param mu
     *            Service rates for the service stages.
     */
    public QueueMEk1(int lambda1, int lambda2, double[] mu) {
        super(new QueueMEk1State(0, mu.length),//
                QueueMEk1Event.getAllEvents(mu.length));
        this.lambda1 = lambda1;
        this.lambda2 = lambda2;
        this.mu = mu;
    }

    /**
     * Used by GUI
     */
    public QueueMEk1() {
        this(3, 2, new double[] { 9, 9 });
    }

    /**
     * The number of possible stages
     * 
     * @return number of stages.
     */
    public int getNumStages() {
        return mu.length;
    }

    /**
     * @see jmarkov.GeomProcess#active(jmarkov.basic.State, Event)
     */
    @Override
    public boolean active(QueueMEk1State i, int absLevel, QueueMEk1Event e) {

        boolean result = false;

        switch (e.type) {
        case ARRIVAL1:
            result = (absLevel == 0);
            break;
        case ARRIVAL2:
            result = (absLevel != 0);
            break;
        case SERVICE1:
            result = (i.getServerStatus() == e.stage);
            break;
        case DEPARTURE:
            result = (i.getServerStatus() == getNumStages());
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
    public GeomRelState<QueueMEk1State>[] dests(QueueMEk1State i, int absLevel,
            QueueMEk1Event e) {
        int newStatus = i.getServerStatus();
        int rLevel = 0;

        switch (e.type) {
        case ARRIVAL1:
            newStatus = 1;
            rLevel = 1;
            break;
        case ARRIVAL2:
            rLevel = 1;// and keeps the current status.
            break;
        case SERVICE1:
            newStatus++;
            break;
        case DEPARTURE:
            if (absLevel == 1) {
                newStatus = 0; // goes idle
            } else {
                rLevel = -1;
                newStatus = 1;// starts service
            }
            break;

        }
        QueueMEk1State newSubState = new QueueMEk1State(newStatus,
                getNumStages());
        GeomRelState<QueueMEk1State> s;
        if (e.type == DEPARTURE && absLevel == 1) {
            s = new GeomRelState<QueueMEk1State>(newSubState);
        } else {
            s = new GeomRelState<QueueMEk1State>(newSubState, rLevel);
        }
        StatesSet<GeomRelState<QueueMEk1State>> statesSet //
        = new StatesSet<GeomRelState<QueueMEk1State>>(s);
        return statesSet.toStateArray();
    }

    /*
     * (non-Javadoc)
     * 
     * @see jmarkov.GeomProcess#rate(jmarkov.State, int)
     */
    @Override
    public double rate(QueueMEk1State i, int iLevel, QueueMEk1State j,
            int jLevel, QueueMEk1Event e) {
        double result = 0;
        switch (e.type) {
        case ARRIVAL1:
            result = lambda1;
            break;
        case ARRIVAL2:
            result = lambda2;
            break;
        case SERVICE1:
            result = mu[e.stage - 1];// index 0 is stage, etc
            break;
        case DEPARTURE:
            result = mu[getNumStages() - 1];
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
        for (int i = 1; i <= getNumStages(); i++) {
            stg += "\nMu " + i + "     = " + mu[i - 1];
        }
        return stg;
    }

    /**
     * Main method
     * @param a Not used
     */
    public static void main(String[] a) {
        QueueMEk1 theQueue = new QueueMEk1(7, 4, new double[] { 20, 20, 20 });
        // QueueMEk1 theQueue = new QueueMEk1(7, 4, new double[] { 9, 9, 9
        // });//unstable
        theQueue.setDebugLevel(0);
        // theQueue.showGUI();
        // GeomState[] estados = theQueue.getGeomStates();
        // System.out.println(estados);
        // DenseMatrix A = theQueue.GetInitialSol();
        theQueue.printAll();

    }

}

/**
 * This class extends PropertiesState and define the properties of a system with
 * one server and a queue with infinite capacity, it has one property, namely
 * the server status. The server goes thru K different exponential stages before
 * completing the task. The status is 0 if the server is idle.
 * 
 * @author Julio Goez. Universidad de los andes.
 */
class QueueMEk1State extends PropertiesState {
    private int k;

    /**
     * We identify the states with is the current service stage of the costumer
     * in service (1,2) or 0 if idle.
     * 
     * @param status
     *            Current stage of service of costumer.
     */

    QueueMEk1State(int status, int numStages) {
        super(1);
        this.k = numStages;
        setProperty(0, status);
    }

    @Override
    public void computeMOPs(MarkovProcess mp) {
        setMOP(mp,"Server Utilization", (getServerStatus() != 0) ? 1 : 0);
        for (int i = 1; i <= k; i++) {
            setMOP(mp,"Stage " + i + " Utilization", //
                    (getServerStatus() == i) ? 1.0 : 0.0);
        }
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
class QueueMEk1Event extends Event {

    final static int ARRIVAL1 = 0; // only for emtpy system.
    final static int ARRIVAL2 = 1; // When there are one or more costumers in
    // the
    // system.
    final static int SERVICE1 = 2; // the first (k-1) stages of service.
    final static int DEPARTURE = 3; // the second stage of service.
    int type;
    int stage; // used only for SERVICE1

    QueueMEk1Event(int type) {
        this.type = type;
    }

    QueueMEk1Event(int type, int stage) {
        this.type = type;
        this.stage = stage;
    }

    static EventsSet<QueueMEk1Event> getAllEvents(int k) {
        EventsSet<QueueMEk1Event> eSet = new EventsSet<QueueMEk1Event>();
        eSet.add(new QueueMEk1Event(ARRIVAL1));
        eSet.add(new QueueMEk1Event(ARRIVAL2));
        for (int i = 1; i < k; i++)
            eSet.add(new QueueMEk1Event(SERVICE1, i));
        eSet.add(new QueueMEk1Event(DEPARTURE));
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
            stg = "End Stage " + stage;
            break;
        case (DEPARTURE):
            stg = "Second stage";
            break;
        }
        return stg;
    }
}
