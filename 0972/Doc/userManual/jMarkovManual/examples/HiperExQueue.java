package examples.jmarkov;

import static examples.jmarkov.HiperExQueueEvent.Type.ARRIVAL;
import static examples.jmarkov.HiperExQueueEvent.Type.FINISH_SERVICE;
import jmarkov.GeomProcess;
import jmarkov.GeomRelState;
import jmarkov.MarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.StatesSet;
import jphase.DenseContPhaseVar;
import jphase.PhaseVar;

/**
 * This class implements a system with one server with PH
 * service. The server has an infinite capacity.
 * @author Julio C. Goez - Germán Riaño. Universidad de los Andes.
 */
public class HiperExQueue extends
        GeomProcess<HiperExQueueState, HiperExQueueEvent> {
    private double lambda;
    private PhaseVar servDist = null; // service distributions.

    // private int[] capacities = null; // cpacities

    /**
     * This class model a series of stations with poisson arrival
     * rate, and Phase type service distributions. There is a limited
     * buffer capacity between stations
     * @param lambda Arrival rate
     * @param dists service distributions
     */
    public HiperExQueue(double lambda, PhaseVar dists) {
        super(new HiperExQueueState(0), //
                HiperExQueueEvent.getAllEvents(dists));
        this.lambda = lambda;
        this.servDist = dists;
    }

    /**
     * Used by GUI.
     */
    public HiperExQueue() {
        this(1.0, DenseContPhaseVar.HyperExpo(new double[] { 8.0,
                        8.0 }, new double[] { 0.5, 0.5 }));
    }

    /**
     * Add to the set all destinations generates by the server
     * starting a new service.
     * @param s
     * @param destStates
     */
    void addDestsFinishServer(int relLevel,
            StatesSet<GeomRelState<HiperExQueueState>> destStates) {
        PhaseVar v = servDist;
        double alpha[] = v.getVectorArray();// start probs
        int m = v.getNumPhases();
        for (int n = 0; n < m; n++) {
            if (alpha[n] > 0) {
                GeomRelState<HiperExQueueState> gs;
                gs = new GeomRelState<HiperExQueueState>(new HiperExQueueState(
                        n+1), relLevel);
                destStates.add(gs);
            }
        }
    }

    /*
     * (non-Javadoc)
     * @see jmarkov.GeomProcess#active(jmarkov.State, int)
     */
    @Override
    public boolean active(HiperExQueueState i, int iLevel, HiperExQueueEvent e) {

        boolean result = false;
        switch (e.type) {
        case ARRIVAL:
            result = true;
            break;
        case FINISH_SERVICE:
            result =  (i.getSrvPhase() == e.getCurPH());
            break;
        }
        return result;
    }

    /*
     * (non-Javadoc)
     * @see jmarkov.GeomProcess#dest(jmarkov.State, int)
     */
    @Override
    public GeomRelState<HiperExQueueState>[] dests(HiperExQueueState i,
            int absLevel, HiperExQueueEvent e) {
        StatesSet<GeomRelState<HiperExQueueState>> destStates //
        = new StatesSet<GeomRelState<HiperExQueueState>>();

        int newPhase = i.getSrvPhase();
        int rLevel = 0;

        switch (e.type) {
        case ARRIVAL:
            rLevel = +1;
            if (absLevel == 0)
                addDestsFinishServer(rLevel, destStates);
            else
                destStates.add(new GeomRelState<HiperExQueueState>(
                        new HiperExQueueState(newPhase), rLevel));
            break;
        case FINISH_SERVICE:
            rLevel = -1;
            if (absLevel == 1)
                destStates.add(new GeomRelState<HiperExQueueState>(
                        new HiperExQueueState(0)));
            else
                addDestsFinishServer(rLevel, destStates);
            break;
        }
        return destStates.toStateArray();
    }// end of dest

    /**
     * This method calculates the rate of transition from i to j when
     * occurs the event e.
     * @param i initial state.
     * @param j final state.
     * @param e event.
     * @return the rate of a transition from i to j when ocurr
     */
@Override
    public double rate(HiperExQueueState i, int iLevel, HiperExQueueState j, int jLevel,
            HiperExQueueEvent e) {
        double rate = -1;
        // gets the absotion rate vector vector
        double[] a = servDist.getMat0Array();
        switch (e.type) {
        case ARRIVAL:
            if(iLevel == 0){
                int newPhase = j.getSrvPhase();
                double alpha[] = servDist.getVectorArray();
                rate = lambda * alpha[newPhase - 1]; 
            }
            else{
                rate = lambda;               
            }
 
            break;
        case FINISH_SERVICE:
            if (iLevel > 1){
                int newPhase = j.getSrvPhase();
                double alpha[] = servDist.getVectorArray();
                rate = a[e.getCurPH() - 1] * alpha[newPhase - 1];
            }
            else{
                rate = a[e.getCurPH() - 1];    
            }
            break;
        }
        return rate;
    }// end of rate
    /*
     * (non-Javadoc)
     * @see jmarkov.SimpleMarkovProcess#description()
     */
    @Override
    public String description() {
        return "Line of PH stations";
    }// end of description.

    /**
     * Main method
     * @param a Not used
     */
    public static void main(String[] a) {
        double lambda = 1.0;
        PhaseVar v1 = DenseContPhaseVar.HyperExpo(new double[] { 8.0, 8.0 },
                new double[] { 0.5, 0.5 });

        HiperExQueue model = new HiperExQueue(lambda, v1);
        model.showGUI();

        model.generate();
        model.setDebugLevel(0);
        model.printAll();

    }// end of main

}// end of class

/**
 * This class defines the events in the queue.
 * @author Julio Goez - German Riaño
 */
class HiperExQueueEvent extends Event {

    /**
     * Enumeration of all posible events
     * @author Julio Goez - German Riaño. Universidad de los Andes.
     *         (C) 2006
     */
    public enum Type {
        /** Arrivals to the system. */
        ARRIVAL,
        /** Finished curPH of server. */
        FINISH_SERVICE
    }

    Type type;// type
    private int curPH;// current phase, which is finished.

    /** Arrival event */
    HiperExQueueEvent() {
        this.type = ARRIVAL;
    }

    /** general event */
    HiperExQueueEvent(Type type, int phase) {
        this.type = type;
        this.curPH = phase;
    }

    /**
     * @return Returns the currrent Phase.
     */
    public int getCurPH() {
        if (type == ARRIVAL)
            throw new IllegalArgumentException(
                    "Current phase is not defined for event " + ARRIVAL);
        return curPH;
    }

    /**
     * @return Returns the type.
     */
    public Type getType() {
        return type;
    }

    static EventsSet<HiperExQueueEvent> getAllEvents(PhaseVar phaseVar) {
        EventsSet<HiperExQueueEvent> E = new EventsSet<HiperExQueueEvent>();
        E.add(new HiperExQueueEvent());// arrival event
        int numPhases = phaseVar.getNumPhases();
        for (int n = 1; n <= numPhases; n++) {
            // finish in phase n
            E.add(new HiperExQueueEvent(FINISH_SERVICE, n));
        }
        return E;
    }

    @Override
    public String label() {
        String stg = "";
        switch (type) {
        case ARRIVAL:
            stg = "Arrival";
            break;
        case FINISH_SERVICE:
            stg = "Ph(" + curPH + ")";
        }
        return stg;
    }
}

/**
 *  * This class define the states in the queue.
 * @author Julio Goez - German Riano. Universidad de los Andes.
 */
class HiperExQueueState extends PropertiesState {

    /**
     * We identify the states with the curPH of server in station, (1,
     * ..,n) or 0 if idle.
     * @param servPahse Service current phase in station.
     */
    public HiperExQueueState(int servPahse) {
        super(1);
        setProperty(0, servPahse);
    }

    @Override
    public void computeMOPs(MarkovProcess mp) {
        setMOP(mp, "Server Utilization", (getSrvPhase() != 0) ? 1 : 0);
    }

    /**
     * Returns the service phase of process
     * @return Service phase
     */
    public int getSrvPhase() {
        return this.prop[0];
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
     * Returns the service status
     * @return Service status (1 = busy, 0 = free)
     */
    public int getSrvStatus() {
        return (getSrvPhase() == 0) ? 0 : 1;
    }

    @Override
    public HiperExQueueState clone() {
        return new HiperExQueueState(getSrvPhase());
    }

    @Override
    public String label() {
        String stg = "";
        stg += "F" + getSrvPhase();
        return stg;
    }

    /*
     * (non-Javadoc)
     * @see jmarkov.State#description()
     */
    @Override
    public String description() {
        return label();
    }

}

