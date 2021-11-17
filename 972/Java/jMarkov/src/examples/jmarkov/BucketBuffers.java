package examples.jmarkov;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import jmarkov.MarkovProcess;
import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import Jama.Matrix;

// CLASS END


/**
 * Model of a Bucket Brigades production system with buffers. 
 * @author Diego Rojas, Germán Riaño.
 */
public class BucketBuffers extends
        SimpleMarkovProcess<BucketStateBuf, BucketEventBuf> {

    private int N;
    private int M;
    private double[][] mu;
    private int[] capBuffers;

    /**
     * Creates a Bucket Brigades system with the given number of
     * workers, machines, process rates and buffer capacities.
     * @param N Number of workers
     * @param M Number of machines
     * @param mu Processing rates for each worker at each machine.
     * @param capBuffers Buffers capacities.
     */
    private BucketBuffers(int N, int M, double[][] mu, int[] capBuffers) {
        super(new BucketStateBuf(N, M), BucketEventBuf.getAllEvents(M));
        this.M = M;
        this.N = N;
        this.mu = mu;
        this.capBuffers = capBuffers;
        assert (mu.length == N);
        assert (mu[0].length == M);
        assert (capBuffers.length == M);
    }

    /**
     * Creates a Bucket Brigades system with the given process rates
     * and buffer capacities.
     * @param mu Processing rates for each worker at each machine.
     * @param capBuffers Buffers capacities.
     */

    public BucketBuffers(double[][] mu, int[] capBuffers) {
        this(mu.length, capBuffers.length, mu, capBuffers);
    }

    /**
     * Creates a Bucket Brigades system with the given process rates
     * and buffer capacities read form the files..
     * @param ratesFile Rates file
     * @param buffersFile Capacity file
     * @throws IOException
     */
    public BucketBuffers(String ratesFile, String buffersFile)
            throws IOException {
        this(loadRateMatrix(ratesFile), loadBuffers(buffersFile));
    }
    
    /**
     * Used by GUI
     * @throws IOException
     */
    public BucketBuffers() throws IOException {
        this("BucketFiles/Mu2.txt","BucketFiles/Buf2.txt");
    }


    /**
     * This boolean function returns 1 if one event e is active (may
     * happen) when the system is in a state i, and 0 otherwise.
     */
    @Override
    public boolean active(BucketStateBuf i, BucketEventBuf e) {
        return (i.getMachineStatus(e.getMachine()) == 1);
    }

    /**
     * This function returns an array of states that are reached from
     * a state i after the event e has ocurred. An auxiliar function
     * getOpEvent(state,machine) returns the number of the worker who
     * finished his work on the machine m. Since this is a pure bucket
     * brigade system, this worker can only be the one with the
     * highest ordinal number stationed at machine m.
     */
    @Override
    public States<BucketStateBuf> dests(BucketStateBuf i, BucketEventBuf e) {
        int[] newMachines = i.getMachines();
        int[] newBuffers = i.getBuffers();
        int evtWorker = getWorkerEvt(i, e.getMachine());

        // First case. Machine m+1 is idle. Worker j moves to machine
        // m+1 and
        // begins work
        if (e.getMachine() < M - 1
                && i.getMachineStatus(e.getMachine() + 1) == 0) {
            newMachines[evtWorker] = e.getMachine() + 1;
        }
        // Second case. Machine m+1 is busy and there is no space in
        // the buffer
        // in front of it.
        // Worker waits on machine m+1
        else if (e.getMachine() < M - 1
                && i.getMachineStatus(e.getMachine() + 1) == 1
                && !spaceInBuffer(e.getMachine() + 1, i)) {
            newMachines[evtWorker] = e.getMachine() + 1;
        }
        // Third case. Machine m+1 is busy and there is space in the
        // buffer in
        // front of it.
        // Worker drops his piece into the buffer of machine m+1 and
        // moves
        // backwards.
        else if (e.getMachine() < M - 1
                && i.getMachineStatus(e.getMachine() + 1) == 1
                && spaceInBuffer(e.getMachine() + 1, i)) {
            newBuffers[e.getMachine() + 1]++;
            return new StatesSet<BucketStateBuf>(i.moveBackwards(evtWorker,
                    newMachines, newBuffers));
        }
        // Fourth case. The last worker has finished on the last
        // machine.
        // This worker goes backwards
        else {
            assert (e.getMachine() == M - 1);
            return new StatesSet<BucketStateBuf>(i.moveBackwards(evtWorker,
                    newMachines, newBuffers));
        }
        return new StatesSet<BucketStateBuf>(new BucketStateBuf(newMachines,
                newBuffers));
    }

    /**
     * Now the rate at which an event occurs, given a state, is
     * defined.
     */
    @Override
    public double rate(BucketStateBuf i, BucketStateBuf j, BucketEventBuf e) {
        return mu[getWorkerEvt(i, e.getMachine())][e.getMachine()];
    }

    /**
     * SpaceInBuffer returns 1 if the ith buffer has space available,
     * and 0 otherwise
     * @param i
     * @return 1 if the ith buffer has space available, and 0
     *         otherwise
     */
    private boolean spaceInBuffer(int i, BucketStateBuf s) {
        return (s.getBuffer(i) < capBuffers[i]);
    }

    /**
     * @see jmarkov.SimpleMarkovProcess#description()
     */
    @Override
    public String description() {
        String stg = "  ";
        stg += "Bucket Brigades Production System with " + this.M
                + " Machines and " + this.N + " Workers\n";
        stg += "\n"
                + "Matrix Mu: Process rate of worker w. on machine m. \n      ";
        for (int i = 0; i < M; i++)
            stg += "M" + i + "  ";

        stg += "\n";
        for (int i = 0; i < N; i++) {
            stg += "N" + i + "  ";
            for (int j = 0; j < M; j++) {
                stg += (mu[i][j]) + " ";
            }
            stg += "\n";
        }
        stg += "\nB.C. ";
        for (int i = 1; i < M; i++) {
            stg += capBuffers[i] + "   ";
        }
        return stg;
    }

    /**
     * getOpEvent returns the number of the worker that has finished
     * on machine ev.machine. Since this is a common BucketBrigades
     * system, the worker will be the one with highest number that was
     * active on ev.machine
     */
    int getWorkerEvt(BucketStateBuf s, int m) {
        int result = 0;
        for (int i = 0; i < N; i++) {
            if (s.getMachine(i) == m)
                result = i;
        }
        return result;
    }

    /**
     * Reads the maximum buffer capacity for each buffer between
     * machines
     * @param fileBuf
     * @return maximum buffer capacity
     * @throws IOException
     */
    private static int[] loadBuffers(String fileBuf) throws IOException {
        Matrix Buf;
        double[] buf = null;

        try {
            Buf = Matrix.read(new BufferedReader(new FileReader(fileBuf)));
            buf = Buf.getRowPackedCopy();
        } catch (IOException e) {
            throw new IOException("Unable to read " + fileBuf);
        }
        int m = buf.length + 1;
        int[] buffy = new int[m];
        for (int i = 1; i < m; i++) {
            buffy[i] = (int) buf[i - 1];
        }

        return buffy;
    }

    /**
     * Reads the process rates for every worker on every machine.
     * @param fileMu
     * @throws IOException
     * @returns the read process rate matrix.
     */
    private static double[][] loadRateMatrix(String fileMu) throws IOException {
        double[][] mu = null;
        Matrix muMatrix;
        try {
            muMatrix = Matrix.read(new BufferedReader(new FileReader(fileMu)));
            mu = muMatrix.getArrayCopy();
        } catch (IOException e) {
            throw new IOException("Unable to read " + fileMu);
        }
        return mu;
    }

    /**
     * Test program.
     * @param args
     * @throws IOException
     */
    public static void main(String args[]) throws IOException {
        String fileMu = "BucketFiles/Mu2.txt";
        String fileBuf = "BucketFiles/Buf2.txt";
        BucketBuffers theModel = new BucketBuffers(fileMu, fileBuf);
   //     theModel.showGUI();
        theModel.setDebugLevel(4);
        theModel.generate();

    }

}

/**
 * The state is described by 4 vectors. The first one, is the machine
 * number every worker is stationed at. The second one, is the status
 * of each worker: 1 if busy, 0 if idle The third one, is the status
 * of each machine: 1 if busy, 0 if idle The last one is the number of
 * units there are in the M-1 possible buffers between machines The
 * second and third vectors can be deduced from the first vector.
 */

class BucketStateBuf extends PropertiesState {
    // Number of machines
    private static int M;
    // Number of workers
    private static int N;

    /**
     * First constructor. Creates the first state i0
     * @param numWorkers
     * @param numMachines
     */

    public BucketStateBuf(int numWorkers, int numMachines) {
        super(2 * numWorkers + 2 * numMachines);
        M = numMachines;
        N = numWorkers;

        for (int i = 0; i < numWorkers; i++) {
            prop[i] = 0; // All workers begin at the first machine
            prop[i + numWorkers] = 0; // (almost) all workers are on
            // hold
        }
        prop[2 * numWorkers - 1] = 1; // Only the last worker is
        // active at the
        // beginning
        for (int i = 1; i < numMachines; i++) {
            prop[2 * numWorkers + i] = 0; // All machines but the
            // first one
            // are
            // inactive
        }
        prop[2 * numWorkers] = 1; // Only the first machine is
        // active
        // Initially, all buffers contain no units.
        // Now the Measures of Performance are set.
    }

    /**
     * Second constructor. Creates all other states
     * @param machines
     * @param workerStatus
     * @param machineStatus
     * @param numInBuffers
     */
    public BucketStateBuf(int[] machines, int[] workerStatus,
            int[] machineStatus, int numInBuffers[]) {
        super(2 * (machines.length) + 2 * machineStatus.length);
        System.arraycopy(machines, 0, prop, 0, N);
        System.arraycopy(workerStatus, 0, prop, N, N);
        System.arraycopy(machineStatus, 0, prop, 2 * N, M);
        System.arraycopy(numInBuffers, 0, prop, 2 * N + M, M);
    }

    /**
     * Third constructor. Creates all other states
     * @param machines
     * @param numInBuffers
     */
    public BucketStateBuf(int[] machines, int numInBuffers[]) {
        // TODO: kill superfluos info. Or may be not?
        this(machines, getWorkerStatus(machines), getMachineStatus(machines),
                numInBuffers);
    }

    
    /**
     * Sets all measures of performance.
     * @see jmarkov.basic.State#computeMOPs(MarkovProcess)
     */
    @Override
    public void computeMOPs(MarkovProcess mp) {
        int sumWrks = 0;
        int sumBuf = 0;

        /*
         * We want to know how much time every worker spends on every
         * machine For this purpose, the function int askMaqOp is
         * created. It returns 1 if the worker i is currently
         * stationed at the machine j, and 0 if it is not
         */

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++)
                setMOP(mp, "Machine " + (j+1) + " used by worker " + (i+1),
                        isWkrAtMachine(i, j) ? 1 : 0);
        }
        // Now, the number of active workers in the system is
        // calculated.
        for (int i = 0; i < N; i++) {
            sumWrks += getWorkerStatus(i);
            setMOP(mp, "Status Worker " + (i+1), getWorkerStatus(i));
        }
        // The average machine utilisation is also calculated
        for (int i = 0; i < M; i++)
            setMOP(mp, "Status Machine" + (i+1), getMachineStatus(i));
        // The average units in buffers are caculated
        for (int i = 1; i < M; i++) {
            setMOP(mp, "Units in Buffer" + (i+1), getBuffer(i));
            sumBuf += getBuffer(i);
        }
        setMOP(mp, "Number of active workers in the system", sumWrks);
        setMOP(mp, "Number of total units in buffer", sumBuf);
    }

    /**
     * The following function returns true if the "worker" is
     * stationed at "machine" and false otherwise
     * @param worker
     * @param machine
     * @return true or false.
     */
    public final boolean isWkrAtMachine(int worker, int machine) {
        return (prop[worker] == machine);
    }

    /**
     * Returns the machine number the k-th worker is stationed at
     * @param k
     * @return the machine number
     */
    public final int getMachine(int k) {
        return prop[k];
    }

    /**
     * Returns an integer vector with the machines where all workers
     * are
     * @return an arrray with the machine numbers.
     */
    public final int[] getMachines() {
        int[] MaqOp = new int[N];
        System.arraycopy(prop, 0, MaqOp, 0, N);
        return MaqOp;
    }

    /**
     * Returns 1 if the k-th worker is active and 0 otherwise
     * @param i worker index.
     * @return Worker status
     */
    public final int getWorkerStatus(int i) {
        return prop[i + N];
    }

    /**
     * Returns 1 if the j-th machine is active and 0 otherwise
     * @param j
     * @return The machine status
     */
    public final int getMachineStatus(int j) {
        return prop[j + 2 * N];
    }

    /**
     * Returns the number of units in each buffer
     * @param j Buffer index
     * @return Number currently in buffer.
     */
    public final int getBuffer(int j) {
        assert (0 <= j && j < M);
        return prop[j + 2 * N + M];
    }

    /**
     * Returns an integer vector with the number of units in each
     * buffer
     * @return Array with number in buffers.
     */
    public final int[] getBuffers() {
        int buffers[] = new int[M];
        System.arraycopy(prop, 2 * N + M, buffers, 0, M);
        return buffers;
    }

    /**
     * Backwards returns the new state after one woker has moved
     * backwards searching for an unit.
     * @param j worker
     * @param machines
     * @param buffers
     * @return The state resulting from a backward movement
     */
    public BucketStateBuf moveBackwards(int j, int[] machines, int[] buffers) {
        int[] newMachines = new int[N];
        int[] newBuffers = new int[M];

        System.arraycopy(machines, 0, newMachines, 0, N);
        System.arraycopy(buffers, 0, newBuffers, 0, M);

        // First case. Worker j-1 is waiting at machine m.
        // Worker j takes the piece that worker j-1 holds. Worker j-1
        // moves then backwards.
        if (j > 0 && machines[j - 1] == machines[j]) {
            return moveBackwards(j - 1, newMachines, newBuffers);
        }
        // Second case. Worker j finds a piece in the buffer for the
        // machine where j is standing and the machine is idle (or the
        // first one).
        // Worker j takes one piece from the buffer and stays on the
        // same machine.
        else if (getBuffer(machines[j]) > 0
                && (getMachineStatus(machines[j]) == 0 || machines[j] == 0)) {
            newBuffers[machines[j]]--;
        }
        /*
         * Third case. Worker j-1 is not waiting on the same machine
         * as j, there is no piece in the buffer at machine m and
         * machine m-1 is idle. Worker j moves back to machine m-1 and
         * then continues backwards.
         */
        else if (j > 0 && machines[j - 1] != machines[j]
                && getMachineStatus(machines[j] - 1) == 0
                && getBuffer(machines[j]) == 0) {
            newMachines[j]--;
            return moveBackwards(j, newMachines, newBuffers);
        }
        /*
         * Fourth case. No piece is in the buffer at machine m and
         * worker j-1 is busy at machine m-1. Worker j uses the
         * preemption rule and preempts worker j-1. Worker j-1 then
         * moves backward.
         */
        else if (j > 0 && getBuffer(machines[j]) == 0
                && getMachineStatus(machines[j] - 1) == 1) {
            newMachines[j]--;
            return moveBackwards(j - 1, newMachines, newBuffers);
        }
        // Fifth case. If the first worker is stationed at the first
        // machine,
        // he can not move backward.
        else if (j == 0 & machines[j] == 0)
            ;

        // Sixth case. If the first worker does not find any units in
        // the buffer
        // at the machine he is on, it not being the first one,
        // then he moves to the previous machine
        else {
            // assert (j == 0 && machines[j] > 0 &&
            // getBuffer(machines[j]) == 0);
            newMachines[j]--;
            return moveBackwards(j, newMachines, newBuffers);
        }

        // Copy of the result vector
        return new BucketStateBuf(newMachines, newBuffers);
    }

    /**
     * This function creates a label for the current state in the
     * state matrix
     */
    @Override
    public String label() {
        String stg = "";
        for (int i = 0; i < N; i++)
            stg += (getMachine(i) + 1);
        /*
         * stg += "/"; for (int i = 0; i < N; i++) stg +=
         * getWorkerStatus(i); stg += "/"; for (int j = 0; j < M; j++)
         * stg += getMachine(j);
         */
        stg += "/";
        for (int j = 1; j < M; j++)
            stg += getBuffer(j);

        return stg;
    }

    /** Provides a description of each state */
    @Override
    public String description() {
        String stg = "Machines (";
        for (int i = 0; i < N; i++)
            stg += (getMachine(i)+1) + "";
        stg += ") Buf. (";
        for (int i = 1; i < M; i++)
            stg += (getBuffer(i)) + "";
        stg += "). Act.Wkrs. = {";
        for (int i = 0; i < N; i++)
            stg += (getWorkerStatus(i) == 1) ? (i+1) + " " : "";
        stg += "} Act. Mach. {";
        for (int j = 0; j < M; j++)
            stg += (getMachineStatus(j) == 1) ? (j+1) + " " : "";
        stg += "}";
        return stg;
    }

    /**
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        // TODO Complete
        boolean result = true;
        return result;
    }

    /**
     * Returns an integer vector where the status of each worker is
     * registered, 1 meaning busy and 0 meaning idle.
     * @param machines
     * @return the workers status.
     */
    private static int[] getWorkerStatus(int[] machines) {
        int[] workersStatus = new int[N];
        workersStatus[N - 1] = 1; // The last worker is always
        // active
        for (int i = N - 2; i >= 0; i--) {
            if (machines[i] != machines[i + 1])
                workersStatus[i] = 1;
        }
        return workersStatus;
    }

    /**
     * Returns a vector with the status of all machines for a new
     * state, given an event e and a state s
     * @param machines
     * @return the machines status
     */
    private static int[] getMachineStatus(int[] machines) {
        int[] StatusMaq = new int[M];
        for (int i = 0; i < N; i++)
            StatusMaq[machines[i]] = 1;
        return StatusMaq;
    }

}

// CLASS END
/**
 * An event occurs when one worker finishes its work on its current
 * machine. In this way, the events can be described by one integer:
 * the machine number where the work has finished.
 */

class BucketEventBuf extends Event {
    private static int M;
    private int machine;

    /**
     * Constructor
     * @param machine The machine that finished processing.
     */
    public BucketEventBuf(int machine) {
        this.machine = machine;
    }

    /** A set with all posible events is created */
    static EventsSet<BucketEventBuf> getAllEvents(int M) {
        EventsSet<BucketEventBuf> E = new EventsSet<BucketEventBuf>();
        for (int i = 0; i < M; i++)
            E.add(new BucketEventBuf(i));
        BucketEventBuf.M = M;
        return E;
    }

    /**
     * @return Returns the machine.
     */
    public final int getMachine() {
        return machine;
    }

    /** This function returns a description of each event */
    @Override
    public String description() {
        return "Process end on machine " + (machine + 1);
    }

    /** This function returns a description of each event */
    @Override
    public String label() {
        return "End(" + (machine + 1) + ")";
    }

}

// CLASS END

