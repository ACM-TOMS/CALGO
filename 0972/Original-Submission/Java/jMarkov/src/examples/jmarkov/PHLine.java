package examples.jmarkov;

import static examples.jmarkov.PHLineEvent.Type.ARRIVAL;
import static examples.jmarkov.PHLineEvent.Type.CHANGE_PHASE;
import static examples.jmarkov.PHLineEvent.Type.FINISH_AND_IDLE;
import static examples.jmarkov.PHLineEvent.Type.FINISH_AND_IDLE_MOVE;
import static examples.jmarkov.PHLineEvent.Type.FINISH_AND_RESTART;
import static examples.jmarkov.PHLineEvent.Type.FINISH_AND_RESTART_MOVE;

import java.io.PrintWriter;

import jmarkov.GeomProcess;
import jmarkov.GeomRelState;
import jmarkov.MarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.StatesSet;
import jphase.DenseContPhaseVar;
import jphase.PhaseVar;
import no.uib.cipr.matrix.Matrix;

/**
 * This class implements a system with two tandem servers, each one with a PH
 * service. The first server has an infinite capacity and the second server has
 * finite capacity.
 * 
 * @author Julio Goez. Universidad de los Andes.
 */
public class PHLine extends GeomProcess<PHLineState, PHLineEvent> {
    private double lambda;
    private PhaseVar[] servDist = null; // service distributions.
    private int[] capacities = null; // capacities

    /**
     * This class model a series of stations with poisson arrival rate, and
     * Phase type service distributions. There is a limited buffer capacity
     * between stations
     * 
     * @param lambda
     *            Arrival rate
     * @param dists
     *            service distributions
     * @param capacity
     *            capacities for buffer 1, 2, ...
     */
    public PHLine(double lambda, PhaseVar[] dists, int capacity[]) {
        super(new PHLineState(new int[dists.length - 1],//
                new int[dists.length]), //
                PHLineEvent.getAllEvents(dists));
        this.lambda = lambda;
        this.servDist = dists;
        this.capacities = capacity;
    }

    /**
     * Used by GUI.
     * 
     */
    public PHLine() {
        super(new PHLineState(new int[2], new int[2]), PHLineEvent
                .getAllEvents(new PhaseVar[] { DenseContPhaseVar.expo(4),
                        DenseContPhaseVar.expo(5) }));//
        this.lambda = 3;
        this.capacities = new int[] { 4 };
    }

    /**
     * The number of stations
     * 
     * @return number of stations.
     */
    public int getNumStations() {
        return servDist.length;
    }

    /**
     * Returns the capacity of the specified station
     * @param station station index
     * @return Capacity of the specified station
     */
    public int getCapacity(int station) {
        return capacities[station - 1];
    }

    /**
     * Add to the set all destinations generates by the server in station s
     * starting a new service.
     * 
     * @param s
     * @param destStates
     */
    void addDestsChangePhase(PHLineState i, int s, int relLevel,
            StatesSet<GeomRelState<PHLineState>> destStates) {
        PhaseVar v = servDist[s];
        double alpha[] = v.getVectorArray();// start probs for this station
        int m = v.getNumPhases();
        for (int n = 1; n <= m; n++) {
            if (alpha[n - 1] > 0) {
                GeomRelState<PHLineState> gs;
                gs = new GeomRelState<PHLineState>(i.changePhase(s, n),
                        relLevel);
                destStates.add(gs);
            }
        }

    }

    /**
     * 
     * Add to the set all destinations generates by the servers in station s1
     * and s2 starting a new service.
     * 
     * @param i
     * @param s1
     * @param s2
     * @param destStates
     */
    void addDestsChangePhase(PHLineState i, int s1, int s2, int relLevel,
            StatesSet<GeomRelState<PHLineState>> destStates) {
        PhaseVar v1 = servDist[s1];
        double alpha1[] = v1.getVectorArray();// start probs for this station
        int m1 = v1.getNumPhases();
        PhaseVar v2 = servDist[s1];
        double alpha2[] = v2.getVectorArray();// start probs for this station
        int m2 = v2.getNumPhases();
        for (int n1 = 1; n1 <= m1; n1++) {
            for (int n2 = 1; n2 <= m2; n2++) {
                if (alpha1[n1 - 1] > 0 && alpha2[n2 - 1] > 0) {
                    GeomRelState<PHLineState> gs;
                    gs = new GeomRelState<PHLineState>(i.changePhase(s1, n1,
                            s2, n2), relLevel);
                    destStates.add(gs);
                }
            }
        }

    }

    /*
     * (non-Javadoc)
     * 
     * @see jmarkov.GeomProcess#active(jmarkov.State, int)
     */
    @Override
    public boolean active(PHLineState i, int iLevel, PHLineEvent e) {

        boolean result = false;
        if (e.type == ARRIVAL)
            return true;
        int s = e.getStation();
        boolean canFinish = (i.getSrvPhase(s) == e.getCurPH());
        int bSize = (s == 0) ? iLevel : i.getBufferSize(s);
        boolean canRestart = (bSize > 0);
        boolean canMove = (s < getNumStations() - 1) ? // is lastStation?
                i.getBufferSize(s + 1) < getCapacity(s + 1) : false;

        switch (e.type) {
        case ARRIVAL:
            result = true;
            break;
        case CHANGE_PHASE:
            result = (i.getSrvPhase(s) == e.getCurPH());
            break;
        case FINISH_AND_IDLE:
            result = canFinish && !canRestart && !canMove;
            break;
        case FINISH_AND_RESTART:
            result = canFinish && canRestart && !canMove;
            break;
        case FINISH_AND_IDLE_MOVE:
            result = canFinish && !canRestart && canMove;
            break;
        case FINISH_AND_RESTART_MOVE:
            result = canFinish && canRestart && canMove;
            break;
        }

        return result;
    }

    /*
     * (non-Javadoc)
     * 
     * @see jmarkov.GeomProcess#dest(jmarkov.State, int)
     */
    @Override
    public GeomRelState<PHLineState>[] dests(PHLineState i, int absLevel,
            PHLineEvent e) {
        StatesSet<GeomRelState<PHLineState>> destStates //
        = new StatesSet<GeomRelState<PHLineState>>();
        int s = e.getStation();// station that finishes
        int rLevel = 0;
        if (e.type == ARRIVAL) {
            rLevel = +1;
        } else if (s == 0) {
            if (e.type == FINISH_AND_IDLE_MOVE
                    || e.type == FINISH_AND_RESTART_MOVE) {
                rLevel = -1;
            }
        }

        switch (e.type) {
        case ARRIVAL:
            if (absLevel == 0) {
                destStates.add(new GeomRelState<PHLineState>(i, rLevel));
            } else {
                addDestsChangePhase(i, s, rLevel, destStates);
            }
            break;
        case CHANGE_PHASE:
            addDestsChangePhase(i, s, rLevel, destStates);
            break;
        case FINISH_AND_IDLE:
            GeomRelState<PHLineState> gs;
            gs = new GeomRelState<PHLineState>(i.setIdle(s), 0);
            destStates.add(gs);
            break;
        case FINISH_AND_RESTART:
            PHLineState j = i.reduceBuffer(s);// reduce buffer
            addDestsChangePhase(j, s, rLevel, destStates);
            break;
        case FINISH_AND_IDLE_MOVE:
            j = i.move(s, s + 1);// move item
            addDestsChangePhase(j, s, s + 1, destStates);// add states
            break;
        case FINISH_AND_RESTART_MOVE:
            j = i.move(s, s + 1);// move item
            j = j.setIdle(s); // set server idle
            addDestsChangePhase(j, s, s + 1, destStates);
            break;
        }

        return destStates.toStateArray();
    }// end of dest

    /**
     * This method calculates the rate of transition from i to j when occurs the
     * event e.
     * 
     * @param i
     *            initial state.
     * @param j
     *            final state.
     * @param e
     *            event.
     * @return the rate of a transition from i to j when ocurr
     */
    @Override
    public double rate(PHLineState i, int iLevel, PHLineState j, int jLevel,
            PHLineEvent e) {
        double rate = -1;
        int s = e.getStation();
        // gets the absotion rate vector vector
        double[] a = servDist[s].getMat0Array();
        switch (e.type) {
        case ARRIVAL:
            rate = lambda;
            break;
        case CHANGE_PHASE:
            int curPhase = e.getCurPH();
            int newPhase = j.getSrvPhase(s);
            Matrix A = servDist[s].getMatrix();
            rate = A.get(curPhase - 1, newPhase - 1);
            break;
        case FINISH_AND_IDLE:
            rate = a[e.getCurPH() - 1];
            break;
        case FINISH_AND_RESTART: {
            newPhase = j.getSrvPhase(s);
            double alpha[] = servDist[s].getVectorArray();
            rate = a[e.getCurPH() - 1] * alpha[newPhase - 1];
            break;
        }
        case FINISH_AND_IDLE_MOVE: {
            double alpha2[] = servDist[s + 1].getVectorArray();
            rate = a[e.getCurPH() - 1] * alpha2[j.getSrvPhase(s + 1) - 1];
            break;
        }
        case FINISH_AND_RESTART_MOVE: {
            newPhase = j.getSrvPhase(s);
            double alpha[] = servDist[s].getVectorArray();
            double alpha2[] = servDist[s + 1].getVectorArray();
            rate = a[e.getCurPH()] * alpha[newPhase - 1]
                    * alpha2[j.getSrvPhase(s + 1) - 1];
            break;
        }
        }
        return rate;
    }// end of rate

    /*
     * (non-Javadoc)
     * 
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
        double lambda = 4.0;
        PhaseVar v1 = DenseContPhaseVar.expo(5);
        PhaseVar v2 = DenseContPhaseVar.expo(6);
        PhaseVar vars[] = { v1, v2 };
        int cap[] = { 3 };
        PHLine model = new PHLine(lambda, vars, cap);
         model.showGUI();

        model.generate();
        model.setDebugLevel(0);
        PrintWriter out = new PrintWriter(System.out, true);
        // theQueue.printDenseMatrix(out);
        out.println();
        model.printAll(out);

    }// end of main

}// end of class

/**
 * This class define the events on the queue.
 * 
 * @author Julio C. Goez
 * 
 */
class PHLineEvent extends Event {
    
    /**
     * Enumeration of all posible events 
     * @author German Riano. Universidad de los Andes. (C) 2006
     *
     */
    public enum Type {
        /** Arrivals to the system. */
        ARRIVAL,
        /** Finished curPH of server. */
        CHANGE_PHASE,
        /** Finished service, server gets idle or blocked */
        FINISH_AND_IDLE,
        /** Finished service, server re starts */
        FINISH_AND_RESTART,
        /** Finished service, server idle, piece moves */
        FINISH_AND_IDLE_MOVE,
        /** Finished service, server re-starts, piece moves */
        FINISH_AND_RESTART_MOVE
    }

    Type type;// type
    private int station; // station number
    private int curPH;// current phase, which is finished.

    /** Arrival event */
    PHLineEvent() {
        this.type = ARRIVAL;
        this.station = 0;// arrival ALWAYS is for station 0.
    }

    /** general event */
    PHLineEvent(Type type, int station, int phase) {
        this.type = type;
        this.curPH = phase;
        this.station = station;
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
     * @return Returns the station.
     */
    public int getStation() {
        return station;
    }

    /**
     * @return Returns the type.
     */
    public Type getType() {
        return type;
    }

    static EventsSet<PHLineEvent> getAllEvents(PhaseVar phaseVars[]) {
        EventsSet<PHLineEvent> E = new EventsSet<PHLineEvent>();
        E.add(new PHLineEvent());// arrival event
        int numStats = phaseVars.length;// num stations
        for (int s = 0; s < numStats; s++) {
            int numPhases = phaseVars[s].getNumPhases();
            for (int n = 1; n <= numPhases; n++) {
                // change form phase n
                E.add(new PHLineEvent(CHANGE_PHASE, s, n));
                // finish in phase n
                E.add(new PHLineEvent(FINISH_AND_IDLE, s, n));
                E.add(new PHLineEvent(FINISH_AND_RESTART, s, n));
                E.add(new PHLineEvent(FINISH_AND_IDLE_MOVE, s, n));
                E.add(new PHLineEvent(FINISH_AND_RESTART_MOVE, s, n));
            }
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
        case CHANGE_PHASE:
            stg = "Ph(" + (station+1) + "," + curPH + ")";
            break;
        case FINISH_AND_IDLE:
            stg = "End(" + (station+1) + "," + curPH + ")";
            break;
        case FINISH_AND_IDLE_MOVE:
            stg = "End(" + (station+1) + "," + curPH + ",M)";
        case FINISH_AND_RESTART:
            stg = "End(" + (station+1) + "," + curPH + ",R)";
        case FINISH_AND_RESTART_MOVE:
            stg = "End(" + (station+1) + "," + curPH + ",R & M)";
        }
        return stg;
    }
}

/**
 * This class defines the states of the queue
 * @author Julio Goez - Germ�n Ria�o. Universidad de los Andes.
 * 
 */
class PHLineState extends PropertiesState {

    /**
     * We identify the states with a vector with the number of costumers in
     * server 2 (0,1,2,3,...), the curPH of server in station 1 (1,2, .., n1), 0
     * if idle or blocked, and the curPH of server in station 2, (0, 1, ..,n2)
     * or 0 if idle.
     * 
     * @param bufferSize
     *            Costumers in each station.
     * @param servPahse
     *            Service current phase in each station 1.
     */
    public PHLineState(int[] bufferSize, int[] servPahse) {
        super(bufferSize.length + servPahse.length);
        // sets prop = [servPhase, bufferSize]
        int k = servPahse.length;//
        System.arraycopy(servPahse, 0, prop, 0, k);// first k posits
        System.arraycopy(bufferSize, 0, prop, k, bufferSize.length);
    }

    @Override
    public void computeMOPs(MarkovProcess mp) {
        int k = getNumStations();
        for (int i = 0; i < k; i++) {
            setMOP(mp,"Server " + (i + 1) + " Utilization",
                    (getSrvPhase(i) != 0) ? 1 : 0);
            if (i > 0) {
                setMOP(mp,"Number in Buffer " + (i + 1), getCostumersInStation(i));
            }
        }
    }

    private static int numStat = -1;

    /**
     * Returns the number of stations in the line
     * @return Number of stations in the line
     */
    public int getNumStations() {
        if (numStat == -1) {
            numStat = (prop.length + 1) / 2;
        }
        return numStat;
    }

    /**
     * Returns the service phase of the specified station
     * @param station station index
     * @return Service phase of the specified station
     */
    public int getSrvPhase(int station) {
        return this.prop[station];
    }

    /**
     * Returns the service phase of all the stations
     * @return Service phase of all the stations
     */
    public int[] getSrvPhase() {
        int k = getNumStations();
        int phases[] = new int[k];
        // first k positions
        System.arraycopy(prop, 0, phases, 0, k);
        return phases;
    }

    /**
     * Returns the number of entities in all the buffers
     * @return Number of entities in all the buffers
     */
    public int[] getBufferSize() {
        int k = getNumStations();
        // get last k-1 positions
        int buff[] = new int[k - 1];
        System.arraycopy(prop, k, buff, 0, k - 1);
        return buff;
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
     * Buffer in form of station s, for s = 1,2, etc. Not used for s==0 , since
     * that is the level. It does not include the one in service.
     * 
     * @param station
     * @return buffer size in station s.
     */
    public int getBufferSize(int station) {
        assert (station > 0);
        return this.prop[getNumStations() + station - 1];
    }
    
    /**
     * Returns number of customers in the specified station 
     * @param station index
     * @return Number of customers in the specified station
     */
    public int getCostumersInStation(int station) {
        return getBufferSize(station) + getSrvStatus(station);
    }
    
    /**
     * Returns the service status of the specified station
     * (1 = busy, 0 = free)
     * @param station station index
     * @return Service status of the specified station 
     * (1 = busy, 0 = free)
     */
    public int getSrvStatus(int station) {
        return (getSrvPhase(station) == 0) ? 0 : 1;
    }

    @Override
    public PHLineState clone() {
        return new PHLineState(getBufferSize(), getSrvPhase());
    }

    /**
     * Creates a new state where it changes station phase form the given phase j
     * 
     * @param station
     * @param j
     * @return a new state with the given change.
     */
    public PHLineState changePhase(int station, int j) {
        int[] bufferSize = getBufferSize();
        int[] srvPhase = getSrvPhase();
        srvPhase[station] = j;
        return new PHLineState(bufferSize, srvPhase);
    }

    /**
     * Creates a new state where server s is idle.
     * 
     * @param station
     * @return a new state with the given change.
     */
    public PHLineState setIdle(int station) {
        int[] bufferSize = getBufferSize();
        int[] srvPhase = getSrvPhase();
        srvPhase[station] = 0;
        return new PHLineState(bufferSize, srvPhase);
    }

    /**
     * Creates a new state where it changes station phase form the given phase j
     * 
     * @param station
     *            station > 0 (not 0)
     * @param b
     *            new buffer size
     * @return a new state with the given change.
     */
    private PHLineState changeBuffer(int station, int b) {
        int[] bufferSize = getBufferSize();
        int[] srvPhase = getSrvPhase();
        bufferSize[station - 1] = b;
        return new PHLineState(bufferSize, srvPhase);
    }
    
    /**
     * Reduce the number of entities in the buffer of the
     * specified station
     * @param station station index
     * @return New state of the line
     */
    public PHLineState reduceBuffer(int station) {
        int[] bufferSize = getBufferSize();
        int[] srvPhase = getSrvPhase();
        if (station > 0)
            bufferSize[station - 1]--;
        return new PHLineState(bufferSize, srvPhase);
    }

    /**
     * Generates a new State, where apiece moves from old station to new
     * station.
     * 
     * @param oldStation
     * @param newStation
     * @return the state
     */
    public PHLineState move(int oldStation, int newStation) {
        int[] bufferSize = getBufferSize();
        int[] srvPhase = getSrvPhase();
        if (oldStation > 0)
            bufferSize[oldStation - 1] = bufferSize[oldStation - 1] - 1;
        bufferSize[newStation - 1] = bufferSize[newStation - 1] + 1;
        return new PHLineState(bufferSize, srvPhase);
    }

    /**
     * changes phase in station s1 to p1, and in s2 to p2
     * 
     * @param s1
     * @param p1
     * @param s2
     * @param p2 *
     * @return a new state with the given change.
     */
    public PHLineState changePhase(int s1, int p1, int s2, int p2) {
        int[] bufferSize = getBufferSize();
        int[] srvPhase = getSrvPhase();
        srvPhase[s1] = p1;
        srvPhase[s2] = p2;
        return new PHLineState(bufferSize, srvPhase);

    }

    @Override
    public String label() {
        String stg = "";
        int k = getNumStations();
        for (int s = 0; s < k; s++) {
            if (s > 0)
                stg += "B" + getBufferSize(s);
            stg += "F" + getSrvPhase(s);
        }
        return stg;
    }

    /*
     * (non-Javadoc)
     * 
     * @see jmarkov.State#description()
     */
    @Override
    public String description() {
        return label();
    }

}
