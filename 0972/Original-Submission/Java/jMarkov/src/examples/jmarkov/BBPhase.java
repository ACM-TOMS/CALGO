package examples.jmarkov;

//BUCKET BRIGADES - JAVA MODELLING
import java.util.ArrayList;

import jmarkov.MarkovProcess;
import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.NotUnichainException;
import jphase.ContPhaseVar;
import jphase.DenseContPhaseVar;
import no.uib.cipr.matrix.Matrix;

/**
 * This is the Java implementation for the Bucket Brigades System
 * using Phase Distribution. This code is based on the Bucket Brigades -
 * Java Modelling created by Diego Rojas
 * @author Juan Pablo Alvarado
 * @author German Riaño, PhD
 * @version 4
 */
class BBPhaseState extends PropertiesState {

    private static int M; // Number of machines in the production
                            // line
    private static int N; // Numeber of workers in the production
                            // line

    /**
     * First Contructor: This creates the inicial state i0 In this
     * State all the workers but the last one are inactive.
     * @param N
     */

    public BBPhaseState(int N) {
        super(2 * N);
        BBPhaseState.N = N;
        prop[N - 1] = 1; // The last Workers is now on Phase 1

        for (int j = N; j < 2 * N; j++) {
            prop[j] = 1; // All workers are on the first machine
        }
        // The first state is has been created
    }

    /**
     * Second Constructor: This creates all other states
     * @param phases
     * @param machines
     */
    BBPhaseState(int[] phases, int[] machines) {
        super(2 * phases.length);
        System.arraycopy(phases, 0, prop, 0, N);
        System.arraycopy(machines, 0, prop, N, N);
    }

    // All other state have been created

    /**
     * The following method returns 1 if the worker "wk" is stationed
     * at the machine "k", 0 otherwise.
     * @param wk
     * @param k
     * @return 1 or 0
     */
    public int askWorkerMachine(int wk, int k) {
        int result = 0;
        if (prop[N + wk] == k && prop[wk] > 0) {
            result = 1;
        }
        return result;
    }

    /**
     * Measures of permormances
     * @author Diego Rojas
     * @version Juan Pablo
     */

    @Override
    public void computeMOPs(MarkovProcess mp) {

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                setMOP(mp,
                        "Machine " + (j + 1) + " usage by worker " + (i + 1),
                        askWorkerMachine(i, j + 1));
            }
        }

        // The number of active workers is calculated
        int sumOp = 0;
        for (int i = 0; i < N; i++) {
            sumOp += getWkrStatus(i);
            setMOP(mp, "Status Worker " + (i + 1), getWkrStatus(i));
        }
        setMOP(mp, "Number  of active workers in the production line", sumOp);
        // The average machine utilization is now calculated
        for (int i = 0; i < M; i++) {
            setMOP(mp, "Status Machine " + (i + 1), getMachStatus(i + 1));
        }

        for (int i = 0; i < N; i++) {
            setMOP(mp, "Blocked worker " + (i + 1), getWksBlocked(i));
        }
    }

    /**
     * This method returns 1 if the worker "wk" is blocked
     * @param wk
     * @return 1 or 0
     */

    public int getWksBlocked(int wk) {
        if (prop[wk] == 0) {
            return 1;
        } else
            return 0;
    }

    /**
     * This Method returns the machine number where the worker "wk" is
     * stationed at
     * @param wk
     * @return prop[wk]
     */
    public int getWrkMach(int wk) {
        return prop[N + wk];
    }

    /**
     * Returns an integer vector with the machines where all workers
     * are
     * @return machines array
     */
    public int[] getWrkMach() {
        int[] opeMach = new int[N];
        System.arraycopy(prop, N, opeMach, 0, N);
        return opeMach;
    }

    /**
     * This Method returns the phase of the worker "wk"
     * @param wk worker
     * @return phase of wk
     */
    public int getWkrPhase(int wk) {
        return prop[wk];
    }

    /**
     * This Method returns a vector with the phases of all workers
     * @return phases array
     */
    public int[] getWkrPhase() {
        int[] opePhase = new int[N];
        System.arraycopy(prop, 0, opePhase, 0, N);
        return opePhase;
    }

    @Override
    public BBPhaseState clone() {
        return new BBPhaseState(getWkrPhase(), getWrkMach());
    }

    /**
     * This Method returns 1 if the k-th worker is active, 0 otherwise
     * @param wk worker
     * @return status
     */
    public int getWkrStatus(int wk) {
        return ((getWkrPhase(wk) > 0) ? 1 : 0);
    }

    /**
     * This method returns 1 if the k-th machine is active, 0
     * otherwise
     * @param k machine number
     * @return 1 if busy
     */
    public int getMachStatus(int k) {
        int result = 0;
        for (int i = N; i < 2 * N && result == 0; i++) {
            if (prop[i] == k)
                result = 1;
        }
        return result;
    }

    /**
     * Creates a new state where it changes worker phase
     * @param curOpe
     * @param newPhase
     * @param curMac
     * @return a new state with the given change
     */

    public BBPhaseState newPhase(int curOpe, int curMac, int newPhase) {
        int[] machOpe = getWrkMach();
        int[] phaseOpe = getWkrPhase();
        machOpe[curOpe] = curMac;
        phaseOpe[curOpe] = newPhase;
        return new BBPhaseState(phaseOpe, machOpe);
    }

    /**
     * Creates a new state where it changes worker machine
     * @param curWkr
     * @param newMach
     * @param newPhase
     * @return A new state with the given changes
     */

    public BBPhaseState newMach(int curWkr, int newMach, int newPhase) {
        int[] machWkr = getWrkMach();
        int[] phaseWkr = getWkrPhase();
        machWkr[curWkr] = newMach;
        phaseWkr[curWkr] = newPhase;
        return new BBPhaseState(phaseWkr, machWkr);
    }

    /**
     * Creates a the state generated after a reset
     * @param newPhases
     * @param newMacs
     * @return BBPhaseState with the given changes
     */

    public BBPhaseState reset(int[] newPhases, int[] newMacs) {
        int[] machines = newMacs;
        int[] phase = newPhases;
        return new BBPhaseState(phase, machines);
    }

    /**
     * Creates the state generated when a worker begins a new process
     * @param curOpe
     * @param fase1
     * @param fase2
     * @param Mac1
     * @param Mac2
     * @return BBPhaseState with the given changes
     */
    public BBPhaseState newMach(int curOpe, int fase1, int fase2, int Mac1,
            int Mac2) {
        int[] machOpe = getWrkMach();
        int[] phaseOpe = getWkrPhase();
        machOpe[curOpe] = Mac1;
        machOpe[curOpe - 1] = Mac2;
        phaseOpe[curOpe] = fase1;
        phaseOpe[curOpe - 1] = fase2;
        return new BBPhaseState(phaseOpe, machOpe);
    }

    /**
     * This function creates a label for the current state in the
     * state matrix
     */
    @Override
    public String label() {
        String stg = "";
        for (int i = 0; i < N; i++)
            stg += prop[i];
        stg += "/";
        for (int i = N; i < 2 * N; i++)
            stg += prop[i];
        return stg;
    }

    @Override
    public String description() {
        String stg = "PhOpe: {";
        for (int k = 0; k < N - 1; k++)
            stg += (getWkrPhase(k)) + ",";
        stg += (getWkrPhase(N - 1));
        stg += "} MachOpe: {";
        for (int j = 0; j < N - 1; j++)
            stg += (getWrkMach(j)) + ",";
        stg += (getWrkMach(N - 1));
        stg += "} Blocked. Wk: {";
        for (int j = 0; j < N; j++)
            stg += (getWkrPhase(j) == 0) ? (j + 1) + "" : "";
        stg += "} ";
        return stg;
    }

    /**
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        // TODO Complete
        return true;
    }
}

// CLASS END

/**
 * An event occurs when a worker finishes her work on its current
 * machine's phase
 * @author Juan Pablo
 */
class BBPhaseEvent extends Event {
    /**
     * Event types
     */
    public enum Type {
        /**
         * Change of phase
         */
        Change_Phase,
        /**
         * worker finishes and switch machines
         */
        Change_Machine
    };

    private int worker;
    private Type type;

    /**
     * Constructor
     * @param curWkr worker
     * @param type type
     */
    public BBPhaseEvent(int curWkr, Type type) {
        this.worker = curWkr;
        this.type = type;
    }

    /**
     * @return Type of Event
     */
    public Type getType() {
        return type;
    }

    /**
     * @return Type of Event.
     */
    public int getWorker() {
        return worker;
    }

    /**
     * A set with all posible events is created
     */
    static EventsSet<BBPhaseEvent> getAllEvents(int numWork) {
        EventsSet<BBPhaseEvent> E = new EventsSet<BBPhaseEvent>();
        for (int i = 0; i < numWork; i++) {
            for (Type type : Type.values()) {
                E.add(new BBPhaseEvent(i, type));
            }
        }
        return E;
    }

    @Override
    public String label() {
        String stg = "";
        switch (getType()) {
        case Change_Phase: {
            stg = "Change Phase worker " + (worker + 1);
            break;
        }
        case Change_Machine: {
            stg = "Change Machine worker " + (worker + 1);
            break;
        }
        }
        return stg;
    }
}

/**
 * This is the Java implementation for the Bucket Brigades System
 * using Phase Distribution. This code is based on the Bucket Brigades -
 * Java Modelling created by Diego Rojas
 * @author Juan Pablo Alvarado
 * @author Germán Riaño, PhD
 * @version 4
 */
public class BBPhase extends SimpleMarkovProcess<BBPhaseState, BBPhaseEvent> {
    private int N; // Number of Workers
    private int M; // Number of machines

    private ContPhaseVar[][] processTimes; // Ph distribution of each
                                            // worker @

    // Constructors
    /**
     * this is a general constructor, but should not be used
     */
    private BBPhase(int Op, int M, ContPhaseVar[][] procTimes) {
        super(new BBPhaseState(Op), BBPhaseEvent.getAllEvents(Op));
        this.N = Op;
        this.M = M;
        this.processTimes = procTimes;

    }

    /**
     * Builds a Bucket Brigade system.
     * @param N Number of worker
     * @param M Number of machines
     * @param machineTimes PH distributions
     * @param vels Vels[i][i] is the relative velocities for workeri
     *        at machine j
     */
    public BBPhase(int N, int M, ContPhaseVar[] machineTimes, double[][] vels) {
        this(N, M, setPhasesMatrix(machineTimes, vels));
    }

    /**
     * This utility method creates the ph var Matrix. This Matrix
     * contains the PH distribution of each worker at each station
     * @return Vars Matrix
     */
    private static ContPhaseVar[][] setPhasesMatrix(
            ContPhaseVar[] machineTimes, double vels[][]) {
        int Op = vels.length;
        int M = machineTimes.length;
        ContPhaseVar[][] processTimes = new ContPhaseVar[Op][M];
        for (int i = 0; i < Op; i++) {
            for (int j = 0; j < M; j++) {
                processTimes[i][j] = machineTimes[j].copy();
                processTimes[i][j] = machineTimes[j].times(1.0 / vels[i][j]);
            }
        }
        return processTimes;
    }

    @Override
    public boolean active(BBPhaseState i, BBPhaseEvent e) {
        boolean result = false;
        int curWorker = e.getWorker();
        int[] curMachines = new int[N];
        curMachines = i.getWrkMach();
        int[] curPhases = new int[N];
        curPhases = i.getWkrPhase();
        if (i.getWkrStatus(curWorker) == 1) {
            switch (e.getType()) {
            case Change_Phase: {
                int m = processTimes[curWorker][curMachines[curWorker] - 1]
                        .getNumPhases();
                result = (m > 1 && curPhases[curWorker] != 0);
                break;
            }
            case Change_Machine: {
                double a[] = processTimes[curWorker][curMachines[curWorker] - 1]
                        .getMat0Array();
                //TODO: es esto un error ??
                result = (curMachines[curWorker] != 0
                        && curPhases[curWorker] != 0 && a[curPhases[curWorker] - 1] > 0);
                break;
            }
            }
        }
        return result;
    } // end active

    /**
     * These are the SimpleMarkovProcess methods.
     */

    @Override
    public double rate(BBPhaseState i, BBPhaseState j, BBPhaseEvent e) {
        double rate = -1;
        int actWkr = e.getWorker();
        int[] curPhs = i.getWkrPhase();
        int[] newPhs = j.getWkrPhase();
        int[] curMacs = i.getWrkMach();
        int[] newMacs = j.getWrkMach();
        int curPh = curPhs[actWkr];
        int newPh = newPhs[actWkr];
        int curMac = curMacs[actWkr];
        int newMac = newMacs[actWkr];
        switch (e.getType()) {
        case Change_Phase: {
            Matrix A = processTimes[actWkr][curMac - 1].getMatrix();
            rate = A.get(curPh - 1, newPh - 1);
            break;
        }
        case Change_Machine: {
            double alfa[] = processTimes[actWkr][newMac - 1].getVectorArray();
            double a[] = processTimes[actWkr][curMac - 1].getMat0Array();
            if (newPh != 0) {
                if (alfa[newPh - 1] == 0) {
                    rate = a[curPh - 1];
                } else {
                    rate = a[curPh - 1] * alfa[newPh - 1];
                }
            } else
                rate = a[curPh - 1];
            break;
        }
        }
        return rate;
    }// end rate

    @Override
    public States<BBPhaseState> dests(BBPhaseState i, BBPhaseEvent e) {
        StatesSet<BBPhaseState> destState = new StatesSet<BBPhaseState>();
        int[] curPhs = i.getWkrPhase();
        int[] curMachines = i.getWrkMach();
        int actOpe = e.getWorker();
        int curMac = curMachines[actOpe];
        switch (e.getType()) {
        case Change_Phase: {
            int curPh = curPhs[actOpe];
            ArrayList<Integer> phases = getNewPhases(actOpe, curMac, curPh,
                    processTimes);
            for (int newPhase : phases) {
                destState.add(i.newPhase(actOpe, curMac, newPhase));
            }
            break;
        }
        case Change_Machine: {
            if (curMac == M) {// the event is a reset
                int[] newPhs = new int[N];
                int[] newMacs = new int[N];
                newMacs = getResetMacs(curMachines);
                newPhs = getResetPhases(curPhs, newMacs);
                destState.add(i.reset(newPhs, newMacs));
            } else {
                if (N > 1) {
                    // case1: actOpe is the last one in the line, and
                    // actOpe is
                    // blocking actOpe-1
                    if (actOpe == N - 1 && curPhs[actOpe - 1] == 0) {
                        int newMac1 = curMac + 1;
                        int newMac2 = curMachines[actOpe - 1];
                        ArrayList<Integer> fases1 = getNewMachPhases(actOpe,
                                newMac1, processTimes);
                        ArrayList<Integer> fases2 = getNewMachPhases(
                                actOpe - 1, newMac2, processTimes);
                        for (Integer newPhs1 : fases1) {
                            for (Integer newPhs2 : fases2) {
                                destState.add(i.newMach(actOpe, newPhs1,
                                        newPhs2, newMac1, newMac2));
                            }
                        }
                    }
                    // case2: The last worker isn't blocking actOpe-1,
                    // he moves
                    // to the next machine
                    else if (actOpe == N - 1 && curPhs[actOpe - 1] > 0) {
                        int newMac = curMac + 1;
                        ArrayList<Integer> fases = getNewMachPhases(actOpe,
                                newMac, processTimes);
                        for (Integer newPhs : fases) {
                            destState.add(i.newMach(actOpe, newMac, newPhs));
                        }
                    }
                    // case3:actOpe is the first one in the line
                    else if (actOpe == 0) {
                        int newMac = curMac + 1;
                        // case 3.1: New machine is idle, first worker
                        // moves to
                        // this station
                        if (i.getMachStatus(newMac) == 0) {
                            ArrayList<Integer> fases = getNewMachPhases(actOpe,
                                    newMac, processTimes);
                            for (Integer newPhs : fases) {
                                destState
                                        .add(i.newMach(actOpe, newMac, newPhs));
                            }
                        }
                        // case 3.2: New machine is being used, first
                        // worker is
                        // now blocked
                        else {
                            int newPhs = 0;
                            destState.add(i.newMach(actOpe, newMac, newPhs));
                        }
                    }
                    // case4: ActOpe is neither the last one not the
                    // first one
                    else {
                        int newMac = curMac + 1;
                        // case4.1: newMac is idle, actOpe moves to
                        // the next
                        // machine
                        if (i.getMachStatus(newMac) == 0) {
                            if (curPhs[actOpe - 1] == 0) {
                                int newMac2 = curMachines[actOpe - 1];
                                ArrayList<Integer> fases = getNewMachPhases(
                                        actOpe, newMac, processTimes);
                                ArrayList<Integer> fases2 = getNewMachPhases(
                                        actOpe - 1, newMac2, processTimes);
                                for (Integer newPhs1 : fases) {
                                    for (Integer newPhs2 : fases2) {
                                        destState.add(i.newMach(actOpe,
                                                newPhs1, newPhs2, newMac,
                                                newMac2));
                                    }
                                }
                            } else {
                                ArrayList<Integer> fases = getNewMachPhases(
                                        actOpe, newMac, processTimes);
                                for (Integer newPhs : fases) {
                                    destState.add(i.newMach(actOpe, newMac,
                                            newPhs));
                                }
                            }
                        }
                        // case5: ActOpe is now blocked
                        else {
                            int newPhs = 0;
                            if (curPhs[actOpe - 1] == 1)// actOpe-1
                                                        // isn't
                                // blocked.
                                destState
                                        .add(i.newMach(actOpe, newMac, newPhs));
                            else {// actOpe is now blocked but
                                    // actOpe-1 is now
                                // active.
                                int newMac2 = curMachines[actOpe - 1];
                                ArrayList<Integer> fases2 = getNewMachPhases(
                                        actOpe - 1, newMac2, processTimes);
                                for (Integer newPhs2 : fases2) {
                                    destState.add(i.newMach(actOpe, newPhs,
                                            newPhs2, newMac, newMac2));
                                }
                            }

                        }
                    }
                } else {// there is only one worker along the line
                    int newMac = curMac + 1;
                    ArrayList<Integer> fases = getNewMachPhases(actOpe, newMac,
                            processTimes);
                    for (Integer newPhs : fases) {
                        destState.add(i.newMach(actOpe, newMac, newPhs));
                    }
                }
            }
        }
        }
        return destState;
    } // end dest

    /**
     * This Method finds the new phase for a given worker
     * @param actOpe
     * @param curMach
     * @param curPhase
     * @param phases
     * @return newPhase
     */
    public ArrayList<Integer> getNewPhases(int actOpe, int curMach,
            int curPhase, ContPhaseVar[][] phases) {
        ArrayList<Integer> newPhases = new ArrayList<Integer>();
        int m = phases[actOpe][curMach - 1].getNumPhases();
        Matrix A = phases[actOpe][curMach - 1].getMatrix();
        for (int i = 1; i <= m; i++) {
            if (A.get(curPhase - 1, i - 1) > 0)
                newPhases.add(i);
        }
        return newPhases;
    }

    /*****************************************************************
     * This Method returns the new phase when a worker starts a work
     * in a new machine
     * @param actOpe
     * @param newMach
     * @param phases
     * @return A list of Phases
     */

    public ArrayList<Integer> getNewMachPhases(int actOpe, int newMach,
            ContPhaseVar[][] phases) {
        ArrayList<Integer> newMachPhs = new ArrayList<Integer>();
        double alpha[] = phases[actOpe][newMach - 1].getVectorArray();
        int m = phases[actOpe][newMach - 1].getNumPhases();
        for (int i = 1; i <= m; i++) {
            if (alpha[i - 1] > 0)
                newMachPhs.add(i);
        }
        return newMachPhs;
    }

    /**
     * This method creates a vector with the new phases after a reset
     * @param curPhs
     * @param curMacs
     * @return newPhases[]
     */
    public int[] getResetPhases(int[] curPhs, int[] curMacs) {
        int newPhs[] = new int[N];
        if (N > 1) {
            for (int i = N - 1; i > 0; i--) {
                if (curPhs[i - 1] > 0)
                    newPhs[i] = curPhs[i - 1];
                else if (curPhs[i - 1] == 0) {
                    ArrayList<Integer> fases = getNewMachPhases(i, curMacs[i],
                            processTimes);
                    for (Integer newPh : fases) {
                        newPhs[i] = newPh;
                    }
                }
            }
            for (int i = 0; i < N - 1; i++) {
                if (curMacs[i] == curMacs[i + 1])
                    newPhs[i] = 0;
                else {
                    ArrayList<Integer> fases = getNewMachPhases(i, curMacs[i],
                            processTimes);
                    for (Integer newPh : fases) {
                        newPhs[i] = newPh;
                    }
                }
            }
        } else {
            ArrayList<Integer> fases = getNewMachPhases(0, curMacs[0],
                    processTimes);
            for (Integer newPhases : fases) {
                newPhs[0] = newPhases;
            }
        }
        return newPhs;
    }

    /**
     * This method creates a vector with the new machines after a
     * reset
     * @param curMacs
     * @return newMacs [];
     */
    public int[] getResetMacs(int[] curMacs) {
        int newMacs[] = new int[N];
        for (int i = N - 1; i > 0; i--) {
            newMacs[i] = curMacs[i - 1];
        }
        newMacs[0] = 1;
        return newMacs;
    }

    /**
     * Returns the steadystate rate of occurrance of a Reset
     * resetRate;
     * @return The long-run reste rate i.e., the throughput of the
     *         system.
     * @throws NotUnichainException
     */

    public double getResetRate() throws NotUnichainException {
        BBPhaseEvent e;
        States<BBPhaseState> dsst;
        States<BBPhaseState> theStates = getStates();
        BBPhaseEvent[] theEvents = getEvents();
        double resetRate = 0.0;
        double pi[] = getSteadyState();
        int E = theEvents.length;
        for (int l = 0; l < E; l++) {
            if (theEvents[l].getWorker() == N - 1) {
                e = theEvents[l];
                int s = 0;
                for (BBPhaseState i : theStates) {
                    if (i.getWrkMach(e.getWorker()) == M) {
                        if (active(i, e)) {
                            dsst = dests(i, e);
                            if (dsst != null) {
                                for (BBPhaseState j : dsst) {
                                    switch (e.getType()) {
                                    case Change_Machine: {
                                        resetRate += pi[s++] * rate(i, j, e);
                                        break;
                                    }
                                    case Change_Phase: {
                                        resetRate += 0.0;
                                        break;
                                    }
                                    }
                                }
                            }
                        }
                    }
                }
            }

        }
        return resetRate;
    }

    @Override
    public String description() {
        String stg = "  ";
        stg += "Bucket Brigades Production System with " + this.M
                + " Machines and " + this.N + " Workers\n";
        stg += "\nMatrix Phases: Process rate of worker w. on machine m.\n ";
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                stg += "{ Worker (" + (i + 1) + "):Station (" + (j + 1)
                        + ") }.\n";
                stg += processTimes[i][j].toString() + "";
            }
        }
        return stg;
    }// end of description.

    /**
     * Runs a stupid trial
     * @param a
     */
    public static void main(String[] a) {
        int Op = 2;
        int M = 2;
        DenseContPhaseVar v = DenseContPhaseVar.Erlang(1.0, 2);
        DenseContPhaseVar[] distMaq = { v, v };
        double[][] vels = { { 1, 1 }, { 2, 3 }, };
        BBPhase model = new BBPhase(Op, M, distMaq, vels);
        // model.showGUI();
        model.generate();
        model.setDebugLevel(0);
        model.printAll();
    }

}
