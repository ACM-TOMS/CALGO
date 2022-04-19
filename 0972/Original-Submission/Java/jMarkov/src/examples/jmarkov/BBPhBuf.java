package examples.jmarkov;

//BUCKET BRIGADES - JAVA MODELING 
//Bucket Brigades system using phase type distributions 

import static examples.jmarkov.BBPhBufEv.Type.RESET;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;

import jmarkov.MarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.Transition;
import jmarkov.basic.Transitions;
import jmarkov.basic.TransitionsSet;
import jmarkov.basic.exceptions.NotUnichainException;
import jphase.ContPhaseVar;
import jphase.DenseContPhaseVar;
import no.uib.cipr.matrix.Matrix;

/**
* This class is a java implementation of Bucket Brigades System
* using Phase Distribution and buffers among the stations. 
* This model is based in the source provided by Juan Pablo Alvarado
* and Germán Riaño. 
* @author Julio Goez. Universidad de los Andes. (C) 2006
*/
class BBPhBufSt extends PropertiesState {
    private static int M; // Number of machines in the production
    // line
    private static int N; // Numeber of workers in the production
    // line
    private static BBPhBuf model = null;

    /**
     * Second Constructor: This creates all other states
     * @param machines
     * @param phases
     * @param machBuffers
     */

    BBPhBufSt(int[] machines, int[] phases, int[] machBuffers) {
        super(2 * phases.length + machBuffers.length);
        N = phases.length;
        assert (N == machines.length);
        M = machBuffers.length + 1;
        System.arraycopy(machines, 0, prop, 0, N);
        System.arraycopy(phases, 0, prop, N, N);
        System.arraycopy(machBuffers, 0, prop, 2 * N, M - 1);
    }// All the states have been created

    /**
     * First Contructor: This creates the inicial state i0 In this
     * State all the works but the last one are inactive All the
     * Machines but the last one are inactive
     * @param numWorkers
     * @param numMachines
     */
    public BBPhBufSt(int numWorkers, int numMachines) {
        this(initStateMachines(numWorkers), initStatePhases(numWorkers),
                new int[numMachines - 1]);
    }

    private static int[] initStateMachines(int N) {
        int[] machines = new int[N];
        // Al start in machine first machine (0)
        return machines;
    }

    private static int[] initStatePhases(int N) {
        int[] phases = new int[N];
        // Only last worker is active: (phase 1)
        phases[N - 1] = 1;
        return phases;
    }

    /**
     * This Method returns a vector with the phases of all workers
     * @return A vector with the phases of all workers
     */
    public int[] getPhases() {
        int[] phases = new int[N];
        System.arraycopy(prop, N, phases, 0, N);
        return phases;
    }

    /**
     * This Method returns the phase of the worker "i"
     * @param i 0-based index for worker
     * @return the phase of the worker "i".
     */
    public int getPhase(int i) {
        assert (0 <= i && i < N);
        return prop[i + N];
    }

    /**
     * This Method returns the machine number of the worker i.
     * @param i index of worker( 0-based )
     * @return current Machine of worker i
     */
    public int getMachine(int i) {
        assert (0 <= i && i < N);
        return prop[i];
    }

    /**
     * Returns an integer vector with the machines where all workers
     * are.
     * @return vector with the machines.
     */
    public int[] getMachines() {
        int[] machines = new int[N];
        System.arraycopy(prop, 0, machines, 0, N);
        return machines;
    }

    /**
     * This method returns the numbers of items in the j-th machine's
     * buffer This is n_j in the paper.
     * @param j index of machine (0-based)
     * @return the numbers of items in the j-th machine's buffer
     */

    public int getBuffer(int j) {
        assert (j >= -1 && j < M - 1);
        if (j == -1)
            return Integer.MAX_VALUE; // by convention.
        return prop[2 * N + j];
    }

    /**
     * The number in buffers from machines s to j, including s but not
     * j.
     * @param s
     * @param j
     * @return teh number in buffers.
     */
    public int getBuffer(int s, int j) {
        assert (s >= -1 && s < j);
        int sum = 0;
        for (int t = s; t < j; t++) {
            sum += getBuffer(t);
        }
        return sum;
    }

    /**
     * This method returns a vector with the number of units in each
     * machine's buffer
     * @return MachBuffer
     */
    public int[] getBuffers() {
        int[] machBuffer = new int[M - 1];
        System.arraycopy(prop, 2 * N, machBuffer, 0, M - 1);
        return machBuffer;
    }

    /**
     * Checks that if more than one worker is at a machine, only one
     * is busy
     * @return true if consistent
     */
    @Override
    public boolean isConsistent() {
        boolean isOK = true;
        for (int i = 1; i < N && isOK; i++) {
            if (getMachine(i) == getMachine(i - 1)) {
                isOK &= !isBusy(i - 1);
                assert (isOK);
            } else {
                assert (getMachine(i - 1) < getMachine(i));
                isOK &= isBusy(i - 1);
                assert (isOK);
            }
        }
        isOK &= isBusy(N - 1);// last worker MUST be busy
        assert (isOK);
        for (int j = 0; j < M; j++) {
            int numAtMachine = getNumAtAmachine(j);
            isOK &= (getMachStatus(j) <= numAtMachine);
            assert (isOK);
            if (numAtMachine >= 1)
                isOK &= (getMachStatus(j) == 1);
            assert (isOK);
        }
        return isOK;
    }

    /**
     * We are going to calculate the measures of permormances
     * @author Diego Rojas
     * @version Juan Pablo
     */

    @Override
    public void computeMOPs(MarkovProcess mp) {

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                setMOP(mp,
                        "Machine " + (j + 1) + " usage by worker " + (i + 1),
                        getWorkAciveAtMachine(i, j + 1));
            }
        }

        int sumOp = 0;
        for (int i = 0; i < N; i++) {
            setMOP(mp, "Worker " + (i + 1) + " Status  ", isBusy(i) ? 1.0 : 0.0);
            sumOp += isBusy(i) ? 1 : 0;
        }
        setMOP(mp, "Number  of active workers in the production line", sumOp);

        // The average machine utilisation is now calculated
        for (int i = 0; i < M; i++) {
            setMOP(mp, "Status Machine " + (i + 1), getMachStatus(i + 1));
        }
        for (int j = 0; j < M - 1; j++) {
            setMOP(mp,
                    "Number of units in buffer after the machine " + (j + 1),
                    getBuffer(j));
        }
        for (int i = 0; i < N; i++) {
            setMOP(mp, "Blocked worker " + (i + 1), getBlocked(i));
        }
    }

    /**
     * The following method returns 1 if the worker "i" is stationed
     * at the machine "j", 0 otherwise.
     * @param i
     * @param j
     * @return 1 or 0
     */
    public int getWorkAciveAtMachine(int i, int j) {
        return (getMachine(i) == j && getPhase(i) > 0) ? 1 : 0;
    }

    /**
     * This method returns 1 if the worker i is blocked
     * @param i worker (number 0- based)
     * @return 1 or 0
     */

    public int getBlocked(int i) {
        return (getPhase(i) == 0) ? 1 : 0;
    }

    /**
     * This Method returns 1 if the i-th worker is active, 0 otherwise
     * @param i
     * @return 1 if the i-th worker is active, 0 otherwise
     */
    public boolean isBusy(int i) {
        return getPhase(i) > 0;
    }

    /**
     * Returns the worker active at machine j. Formally, it returns
     * max{i:m_i<= j}.
     * @param j Machine number
     * @return The worker index.
     */
    public int getWorker(int j) {
        int i = N - 1;
        while (getMachine(i) > j) {
            i--;
        }
        return i;
    }

    /**
     * This method returns 1 if the j-th machine is active, 0
     * otherwise
     * @param j
     * @return 1 if the j-th machine is active, 0 otherwise
     */
    public int getMachStatus(int j) {
        int result = 0;
        for (int i = 0; i < N && result == 0; i++) {
            if (getWorkAciveAtMachine(i, j) == 1)
                result = 1;
        }
        return result;
    }

    /**
     * Count the number of worker at the machine j
     * @param j Machine number
     * @return Number of Workers.
     */
    public int getNumAtAmachine(int j) {
        int result = 0;
        for (int i = 0; i < N; i++) {
            result += (getMachine(i) == j) ? 1 : 0;
        }
        return result;
    }

    /*****************************************************************
     * TRANSFORMATIONS
     * *******************************************************************
     */

    /**
     * Creates a new state where curOpe has changed his phase. This is
     * T^F_{ik}(x) in the paper.
     * @param curOpe
     * @param newPhase
     * @return a new state with the given change
     */

    public BBPhBufSt phaseChange(int curOpe, int newPhase) {
        int[] newMachines = getMachines();
        int[] newPhases = getPhases();
        int[] newBuffers = getBuffers();
        newPhases[curOpe] = newPhase;
        return new BBPhBufSt(newMachines, newPhases, newBuffers);
    }

    /**
     * Creates a new state where worker i has changed his machine, and
     * starts in the next machine in phase k2. This is T^M_{ik}(x) in
     * the paper.
     * @param i
     * @param k
     * @return a new state with the given change
     */

    public BBPhBufSt moveToNext(int i, int k) {
        int[] newMachines = getMachines();
        int[] newPhases = getPhases();
        int[] newBuffers = getBuffers();
        newMachines[i] = newMachines[i] + 1;
        newPhases[i] = k;
        return new BBPhBufSt(newMachines, newPhases, newBuffers);
    }

    /**
     * Creates a new state where i has changed his machine and gets
     * blocked. This is T^B_{ik}(x) in the paper.
     * @param curOpe current worker
     * @return a new state with the given change
     */

    public BBPhBufSt block(int curOpe) {
        return moveToNext(curOpe, 0);
    }

    /**
     * Creates a new state where worker i walks back, and worker w
     * starts in phase k. It is T^R_{isk}(x) in the paper.
     * @param i worker that initiates the reset
     * @param w worker that
     * @param s
     * @return the new state
     */
    public BBPhBufSt reset(int i, int w, int s) {
        int[] newBuffers = getBuffers();
        int j = getMachine(i);
        // drop item:
        if (j < M - 1)
            newBuffers[j]++;
        // shift everyone:
        int[] newMachines = shift(w, i, getMachines(), s);
        int[] newPhases = shift(w, i, getPhases(), 0);
        // Pick up new item
        if (s > 0)
            newBuffers[s - 1]--;
        return new BBPhBufSt(newMachines, newPhases, newBuffers);
    }

    /**
     * This is T_{wi}(ve, val) in the paper
     * @param low lower index worker
     * @param high high index worker
     * @param vec array
     * @param value new value
     * @return an array resulting from deleting vec[i], shifting from
     *         w to i, puttinf val in vec[w].
     */
    private int[] shift(int low, int high, int[] vec, int value) {
        assert (vec.length == N);
        int[] result = new int[N];
        System.arraycopy(vec, 0, result, 0, N);
        for (int k = low; k < high; k++) {
            result[k + 1] = vec[k];
        }
        result[low] = value;
        return result;
    }

    /**
     * Return the worker w that satisfies w = min {i : machine[i]>=s }
     * @param s
     * @return the worker w.
     */
    public int getResetWorker(int s) {
        int i = 0;
        while (i < N && getMachine(i) < s)
            i++;
        return i;
    }

    /**
     * Returns the first empty machine with material for a reset
     * originated at machine j.
     * @param j The machine
     * @return The last machine it propagates to.
     */
    public int getResetMachine(int j) {
        int s = j;
        while (s > 0
                //&& (getMachStatus(s) == 1 || getBuffer(s-1) == 0 || getNumAtAmachine(s) > 0)) {
                && (getMachStatus(s) == 1 || getBuffer(s - 1) == 0 || getNumAtAmachine(s) > 1)) {
            s--;
        }
        return s;
    }

    @Override
    public BBPhBufSt clone() {
        return new BBPhBufSt(getMachines(), getPhases(), getBuffers());
    }

    /** Separator */
    String sp(int i, int N) {
        return (N < 10) ? "" : (i < N) ? "," : "";
    }

    /**
     * This function creates a label for the current state in the
     * state matrix
     */
    @Override
    public String label() {
        String stg = "";
        for (int i = 0; i < N; i++)
            stg += getMachine(i) + 1 + sp(i, N);
        stg += "/";
        for (int i = 0; i < N; i++)
            stg += getPhase(i) + sp(i, N);
        if (M > 1)
            stg += "/";
        for (int j = 0; j < M - 1; j++)
            stg += getBuffer(j) + sp(j, M - 1);
        return stg;
    }

    @Override
    public String description() {
        String stg = "Machines: (";
        for (int i = 0; i < N; i++)
            stg += (getMachine(i) + 1) + ((i < N - 1) ? "," : "");
        stg += ") Phases: (";
        for (int i = 0; i < N; i++)
            stg += (getPhase(i)) + ((i < N - 1) ? "," : "");
        stg += ") Buffer Units: (";
        for (int j = 0; j < M - 1; j++)
            stg += getBuffer(j) + ((j < M - 2) ? "," : "");
        stg += "), Blocked Wks: (";
        for (int i = 0; i < N; i++)
            stg += (getPhase(i) == 0) ? (i + 1) + "" : "";
        stg += ")";
        return stg;
    }

    /**
     * @return Returns the model.
     */
    public static final BBPhBuf getModel() {
        return model;
    }

    /**
     * @param newModel The model to set.
     */
    static final void setModel(BBPhBuf newModel) {
        model = newModel;
    }
}

/**
 * An event occurs when a worker finishes its work on its current
 * machine's phase
 * @author German Riano, Juan Pablo Alvarado.
 */
class BBPhBufEv extends jmarkov.basic.Event {
    /**
     * Event types.
     */
    enum Type {
        /** Change of phase */
        CHANGE_PHASE,
        /** Worker moves to the next machine */
        MOVE,
        /** Worker moves to the next machine and a new one starts */
        MOVE_START,
        /** Worker moves and gets blocked */
        BLOCKING,
        /**
         * A worker finishes and gets blocked, the previous one starts
         * a job
         */
        BLOCKING_START,
        /** A reset or wlk-back is initiated */
        RESET
    };

    private int curWorker;// current worker
    private Type type; // Event type

    /**
     * Constructor for CHANGE_PHASE, MOVE, BLOCKING
     * @param curOp Current oparator.
     * @param type Type of Event.
     */
    public BBPhBufEv(int curOp, Type type) {
        this.curWorker = curOp;
        this.type = type;
    }

    /**
     * @return the type
     */
    public Type getType() {
        return type;
    }

    /**
     * @return The worker who produces this event.
     */
    public int getWorker() {
        return curWorker;
    }

    /**
     * A set with all posible events is created
     */
    static EventsSet<BBPhBufEv> getAllEvents(int numWorkers) {
        EventsSet<BBPhBufEv> E = new EventsSet<BBPhBufEv>();
        for (int i = 0; i < numWorkers; i++) {
            for (Type type : Type.values()) {
                E.add(new BBPhBufEv(i, type));
            }
        }
        return E;
    }

    @Override
    public String label() {
        String stg = "";
        switch (getType()) {
        case CHANGE_PHASE:
            stg = "PhChge(" + (curWorker + 1) + ")";
            break;
        case MOVE:
            stg = "Move(" + (curWorker + 1) + ")";
            break;
        case MOVE_START:
            stg = "MoveStart(" + (curWorker + 1) + ")";
            break;
        case BLOCKING:
            stg = "Block(" + (curWorker + 1) + ")";
            break;
        case BLOCKING_START:
            stg = "BlockStart(" + (curWorker + 1) + ")";
            break;
        case RESET:
            stg = "Reset(" + (curWorker + 1) + ")";
            break;
        }
        return stg;
    }
}

/**
* This class is a java implementation of Bucket Brigades System
* using Phase Distribution and buffers among the stations. 
* This model is based in the source provided by Juan Pablo Alvarado
* and Germán Riaño. 
* @author Julio Goez. Universidad de los Andes. (C) 2006
*/
public class BBPhBuf extends MarkovProcess<BBPhBufSt, BBPhBufEv> {
    private int N; // Number of Workers
    private int M; // Number of machines
    // PH distribution of each worker @ each statation:
    private ContPhaseVar[] machineTimes;
    private double vels[][];
    private int[] capBuffer;

    /**
     * Builds a Bucket brigades system.
     * @param numWorkers Number of workers
     * @param numMachines Number of machines
     * @param processTimes M- dimensional vector of distributions.
     * @param vels array of velocities (N x M)
     * @param capBuffer Buffer capacity. (size M-1)
     */
    public BBPhBuf(int numWorkers, int numMachines,
            ContPhaseVar[] processTimes, double[][] vels, int[] capBuffer) {
        super(new BBPhBufSt(numWorkers, numMachines), BBPhBufEv
                .getAllEvents(numWorkers));
        this.N = numWorkers;
        this.M = numMachines;
        this.capBuffer = capBuffer;
        this.machineTimes = processTimes;
        this.vels = vels;
        assert (M == processTimes.length);
        assert (N == vels.length);
        assert (M - 1 == capBuffer.length);
        buildReachablePhases();
        BBPhBufSt.setModel(this);
    }

    /**
     * Used by GUI.
     */
    public BBPhBuf() {
        this(2, 2,//
                new ContPhaseVar[] //
                { DenseContPhaseVar.Erlang(1.0, 2),
                        DenseContPhaseVar.Erlang(1.0, 2) }, //
                new double[][] { { 1.0, 1.0 }, { 2.0, 3.0 } }, //
                new int[] { 0 });
        setDebugLevel(0);
        setMaxStates(4000);
    }

    // See below...
    private int reachablePhases[][][];
    private int startPhases[][];
    private boolean canFinish[][];
    private boolean canChange[][];

    /**
     * This method builds the structure of reachablePhases, canStart,
     * canFinish and canChange. reachablePhases[j][k] is an array of
     * all rechable phases from phase k for the distribution that
     * characterizes machine j. canStart[j][k] is true if machine j
     * can start processing on phase k. canFinish[j][k] is true if
     * machine j can finish processing on phase k. canChange[j][k] is
     * tue if machine j can change phase on phase k.
     */
    private void buildReachablePhases() {
        reachablePhases = new int[M][][];
        startPhases = new int[M][];
        canFinish = new boolean[M][];
        canChange = new boolean[M][];
        for (int j = 0; j < M; j++) {
            ContPhaseVar var = machineTimes[j];
            int m = var.getNumPhases();
            double alpha[] = var.getVectorArray();
            double a[] = var.getMat0Array();
            double A[][] = var.getMatrixArray();
            // this are m+1 to count the zero phase case
            // which is initialized to false by default.
            reachablePhases[j] = new int[m + 1][];
            canFinish[j] = new boolean[m + 1];
            canChange[j] = new boolean[m + 1];
            ArrayList<Integer> listStart = new ArrayList<Integer>();
            for (int k = 1; k <= m; k++) {
                if (alpha[k - 1] > 0)
                    listStart.add(k);
                canFinish[j][k] = (a[k - 1] > 0);
                ArrayList<Integer> listReachable = new ArrayList<Integer>();
                for (int k2 = 1; k2 <= m; k2++) {
                    if (A[k - 1][k2 - 1] > 0) {
                        canChange[j][k] |= true;
                        listReachable.add(k2);
                    }
                }

                int n = listReachable.size();
                reachablePhases[j][k] = new int[n];
                for (int k2 = 0; k2 < n; k2++) {
                    reachablePhases[j][k][k2] = listReachable.get(k2);
                }
            }
            int numStart = listStart.size();
            startPhases[j] = new int[numStart];
            for (int k = 0; k < numStart; k++) {
                startPhases[j][k] = listStart.get(k);
            }
        }
    }

    /**
     * Returns a list of starting phases for machine j
     * @param j machine number.
     * @return A list of possible phases.
     */
    public int[] getStartingPhases(int j) {
        assert (0 <= j && j < M);
        if (startPhases == null || startPhases[j] == null)
            buildReachablePhases();
        return startPhases[j];
    }

    /**
     * Returns the reachable phases when the machine j is processing
     * on phase k.
     * @param j Machine number
     * @param k Current Phase
     * @return List of reachble phases for machine j when in phase k.
     */
    public int[] getReachablePhases(int j, int k) {
        assert (0 <= j && j < M);
        assert (1 <= k && k <= machineTimes[j].getNumPhases());
        if (reachablePhases == null || reachablePhases[j] == null
                || reachablePhases[j][k] == null)
            buildReachablePhases();
        return reachablePhases[j][k];

    }

    /**
     * Gets the starting probability at machine j for this phase.
     * @param j Machine
     * @param k Phase
     * @return starting probability: alpha_j(k)
     */
    public double getStartingProb(int j, int k) {
        assert (0 <= j && j < M);
        assert (1 <= k && k <= machineTimes[j].getNumPhases());
        return machineTimes[j].getVector().get(k - 1);
    }

    /**
     * Gets the transition rate for machine j, from phases f1 to f2.
     * @param j Machine.
     * @param k1 origin phase
     * @param k2 destination phase
     * @return Rate. A_j(k1,k2).
     */
    public double getTransitionRate(int j, int k1, int k2) {
        assert (0 <= j && j < M);
        assert (1 <= k1 && k1 <= machineTimes[j].getNumPhases());
        assert (1 <= k2 && k2 <= machineTimes[j].getNumPhases());
        return machineTimes[j].getMatrix().get(k1 - 1, k2 - 1);
    }

    /**
     * The absortion rate from phase f when machine j is working.
     * @param j machine
     * @param f Phase
     * @return Rate a_j(f), where a_j = -A1.
     */
    public double getAbsorptionRate(int j, int f) {
        assert (0 <= j && j < M);
        assert (f <= machineTimes[j].getNumPhases());
        return machineTimes[j].getMat0().get(f - 1);
    }

    /**
     * This Boolean function returns true if one events e is active
     * when the system is in state i, and false otherwise
     * @param x Current state
     * @param e Current Event
     * @return True if the event can occur.
     */

    public boolean active(BBPhBufSt x, BBPhBufEv e) {
        boolean result = false;
        int i = e.getWorker();
        int k = x.getPhase(i);
        int j = x.getMachine(i);
        boolean isLastMach = (j == M - 1);
        boolean nextIsBusy = !isLastMach && (i < N - 1)
                && x.getMachine(i + 1) == j + 1;
        boolean isFirstWorker = (i == 0);
        boolean prevIsBlocked = (!isFirstWorker && (x.getPhase(i - 1) == 0));
        boolean thereIsRoom = (j < M - 1 && x.getBuffer(j) < capBuffer[j]);
        switch (e.getType()) {
        case CHANGE_PHASE:
            result = canChange[j][k];
            break;
        case MOVE:
            result = !isLastMach && !prevIsBlocked && canFinish[j][k]
                    && !nextIsBusy;
            break;
        case MOVE_START:
            result = !isLastMach && prevIsBlocked && canFinish[j][k]
                    && !nextIsBusy;
            break;
        case BLOCKING:
            result = !prevIsBlocked && canFinish[j][k] && nextIsBusy
                    && !thereIsRoom;
            break;
        case BLOCKING_START:
            result = prevIsBlocked && canFinish[j][k] && nextIsBusy
                    && !thereIsRoom;
            ;
            break;
        case RESET:
            result = canFinish[j][k];
            if (!isLastMach) {
                result &= nextIsBusy && thereIsRoom;
            }
            break;
        }
        return result;
    } // end active

    /**
     * @see jmarkov.MarkovProcess#activeTransitions(State, Event)
     */
    @Override
    public Transitions<BBPhBufSt> activeTransitions(BBPhBufSt x, BBPhBufEv e) {
        TransitionsSet<BBPhBufSt> trans = new TransitionsSet<BBPhBufSt>();
        if (!active(x, e))
            return trans;
        int i = e.getWorker();
        int j = x.getMachine(i);
        int k = x.getPhase(i);
        double rate = -1;
        switch (e.getType()) {
        case CHANGE_PHASE:
            for (int k2 : getReachablePhases(j, k)) {
                rate = vels[i][j] * getTransitionRate(j, k, k2);
                trans.add(x.phaseChange(i, k2), rate);
            }
            break;
        case MOVE:
            for (int k1 : getStartingPhases(j + 1)) {
                rate = vels[i][j] * getAbsorptionRate(j, k)
                        * getStartingProb(j + 1, k1);
                trans.add(x.moveToNext(i, k1), rate);
            }
            break;
        case MOVE_START:
            for (int k1 : startPhases[j + 1]) {
                for (int k2 : startPhases[j]) {
                    rate = vels[i][j] * getAbsorptionRate(j, k)
                            * getStartingProb(j + 1, k1)
                            * getStartingProb(j, k2);
                    trans.add(x.moveToNext(i, k1).phaseChange(i - 1, k2), rate);
                }
            }
            break;
        case BLOCKING:
            rate = vels[i][j] * getAbsorptionRate(j, k);
            trans.add(x.block(i), rate);
            break;
        case BLOCKING_START:
            for (int k2 : startPhases[j]) {
                rate = vels[i][j] * getAbsorptionRate(j, k)
                        * getStartingProb(j, k2);
                trans.add(x.block(i).phaseChange(i - 1, k2), rate);
            }
            break;
        case RESET:
            BBPhBufSt y = x.phaseChange(i, 0);// idle the worker
            int s = y.getResetMachine(j);// find the machine
            int w = y.getResetWorker(s); // and worker
            BBPhBufSt z = x.reset(i, w, s);// do reset
            rate = vels[i][j] * getAbsorptionRate(j, k);
            Transition<BBPhBufSt> tr = new Transition<BBPhBufSt>(z, rate);
            trans.add(revive(tr, i, w));
            break;
        }
        return trans;
    }

    /**
     * This method calculates the rate of transition from i to j when
     * occurs the event ev.
     * @param x initial state.
     * @param y final state.
     * @param e event.
     * @return the rate of a transition from i to j when ocurr
     */

    public double rate(BBPhBufSt x, BBPhBufSt y, BBPhBufEv e) {
        int i = e.getWorker();
        int k = x.getPhase(i);
        int k2 = y.getPhase(i);// k'on the paper
        int j = x.getMachine(i);
        double a[] = null;
        double alpha[] = null;
        double rate = -1;
        switch (e.getType()) {
        case CHANGE_PHASE:
            Matrix A = machineTimes[j].getMatrix();
            rate = vels[i][j] * A.get(k - 1, k2 - 1);
            break;
        case MOVE:
            a = machineTimes[j].getMat0Array();
            alpha = machineTimes[j + 1].getVectorArray();
            rate = vels[i][j] * a[k - 1] * alpha[k2 - 1];
            break;
        case MOVE_START:
            a = machineTimes[j].getMat0Array();
            alpha = machineTimes[j + 1].getVectorArray();
            double[] alpha1 = machineTimes[j].getVectorArray();
            int k1 = y.getPhase(i - 1);
            rate = vels[i][j] * a[k - 1] * alpha1[k1 - 1] * alpha[k2 - 1];
            break;
        case BLOCKING:
            a = machineTimes[j].getMat0Array();
            rate = vels[i][j] * a[k - 1];
            break;
        case BLOCKING_START:
            a = machineTimes[j].getMat0Array();
            alpha1 = machineTimes[j].getVectorArray();
            k1 = y.getPhase(i - 1);
            rate = vels[i][j] * a[k - 1] * alpha1[k1 - 1];
            break;
        case RESET:
            int s = x.getResetMachine(j);
            int w = x.getResetWorker(s);
            a = machineTimes[j].getMat0Array();
            alpha = machineTimes[w].getVectorArray();
            rate = vels[i][j] * a[k - 1];
            if (w == 0 && x.getMachine(0) == 0 && x.isBusy(0)) {
                rate = vels[i][j] * a[k - 1];
            } else {
                k2 = y.getPhase(w);// k'on the paper
                rate = vels[i][j] * a[k - 1] * alpha[k2 - 1];
            }
            break;
        }
        return rate;
    }// end rate

    private Transitions<BBPhBufSt> revive(Transitions<BBPhBufSt> trs, int i) {
        TransitionsSet<BBPhBufSt> theSet = new TransitionsSet<BBPhBufSt>();
        for (Transition<BBPhBufSt> tr : trs) {
            theSet.add(revive(tr, i));
        }
        return theSet;
    }

    /**
     * If workers i and w can start working this returns all posiible
     * new states.
     * @param tr Transition.
     * @param i worker
     * @return Al possible states generated when worker i starts
     *         working
     */
    private Transitions<BBPhBufSt> revive(Transition<BBPhBufSt> tr, int i, int w) {
        return revive(revive(tr, i), w);
    }

    /**
     * If worker i can start working this returns all posiible new
     * states.
     * @param tr Transition.
     * @param i worker
     * @return Al possible states generated when worker i starts
     *         working
     */
    public Transitions<BBPhBufSt> revive(Transition<BBPhBufSt> tr, int i) {
        assert (0 <= i && i < N);
        BBPhBufSt x = tr.getState();
        double baseRate = tr.getRate();
        Transitions<BBPhBufSt> transSet = new TransitionsSet<BBPhBufSt>();
        int j = x.getMachine(i);
        boolean isFirstAtMachine = (i < N - 1 && x.getMachine((i + 1)) > j)
                || (i == N - 1);
        if (isFirstAtMachine && !x.isBusy(i)) {
            for (int f : getStartingPhases(j)) {
                transSet.add(x.phaseChange(i, f), baseRate
                        * getStartingProb(j, f));
            }
        } else {
            transSet.add(x, baseRate);
        }
        return transSet;
    }

    /**
     * Returns the steady-state rate of occurrance of a Reset (i.e.,
     * throughput);
     * @return The rate.
     * @throws NotUnichainException
     */

    public double getResetRate() throws NotUnichainException {
        Transitions<BBPhBufSt> trans;
        States<BBPhBufSt> theStates = getStates();
        BBPhBufEv[] theEvents = getEvents();
        double resetRate = 0.0;
        double pi[] = getSteadyState();
        for (BBPhBufEv e : theEvents) {
            if (e.getWorker() == N - 1 && e.getType() == RESET) {
                int s = 0;
                for (BBPhBufSt i : theStates) {
                    if (i.getMachine(e.getWorker()) == M - 1) {
                        if (active(i, e)) {
                            trans = activeTransitions(i, e);
                            if (trans != null) {
                                for (Transition<BBPhBufSt> t : trans) {
                                    resetRate += pi[s] * t.getRate();
                                }
                            }
                        }
                    }
                    s++;
                }
            }

        }
        return resetRate;
    }

    @Override
    public String description() {
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw);
        pw.print("Bucket Brigades Production System with " + M
                + " Machines and " + N + " Workers.\n");
        pw.print("\nBuffers Capacity (");
        for (int j = 0; j < M - 1; j++)
            pw.print(capBuffer[j] + ((j < M - 2) ? "," : ""));
        pw.print(")\n");
        pw.print("\nDistributions: Process Time on machine m.\n");
        for (int j = 0; j < M; j++) {
            pw.println("Machine " + (j + 1) + ": ");
            pw.println(machineTimes[j].description());
        }
        pw.print("\nRelitive Velocities:\n");
        Jama.Matrix mat = new Jama.Matrix(vels);
        mat.print(pw, 8, 2);
        return sw.toString();
    }// end of description.

    @Override
    public int printMOPs(PrintWriter out, int width, int decimals) {
        int namesWidth = super.printMOPs(out, width, decimals);
        try {
            out.printf(pad("Throughput Rate", namesWidth, false)
                    + pad(getResetRate(), width, decimals));
        } catch (NotUnichainException e) {
            out.println(e);
        }
        return namesWidth;
    }

    /**
     * @see jmarkov.MarkovProcess#finalize()
     */
    @Override
    protected void finalize() throws Throwable {
        super.finalize();
        BBPhBufSt.setModel(null);
    }

    /**
     * Test Program.
     * @param a
     */
    public static void main(String[] a) {
        int N = 3;
        int M = 2;
        int[] bufferCaps = { 1 };
        //ContPhaseVar v = DenseContPhaseVar.Erlang(1.0, 2);
        ContPhaseVar v = DenseContPhaseVar.expo(1.0);
        ContPhaseVar v2 = DenseContPhaseVar.expo(1.0);
        ContPhaseVar[] machineTimes = { v, v, v2 };
        double[][] velocities = { { 1, 1 }, { 1, 1}, { 1, 1}};
        BBPhBuf model = new BBPhBuf(N, M, machineTimes, velocities, bufferCaps);
        model.setDebugLevel(4);
        // model.setSteadyStateSolver(new JamaSolver(model));
        //model.showGUI();
        model.generate();
        model.setDebugLevel(4);
        //model.printAll();
    }

}
