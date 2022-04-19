package examples.jmarkov;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import jmarkov.MarkovProcess;
import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import Jama.Matrix;

/**
 * This class describes a Bucket Brigades system with no intermediate 
 * buffers, N workers and M machines with possibly uneven working velocities.
 * @author Diego Rojas, Germán Riaño.
 */
public class BucketBrigades extends
        SimpleMarkovProcess<BucketState, BucketEvent> {
    private int N;// Number of workers
    private int M; // Number of machines

    private double[][] mu;

    // Constructor

    /**
     * @param N Number of workers
     * @param M Number of machines
     * @param mu Processing rates for each worker at each machine.
     */
    public BucketBrigades(int N, int M, double[][] mu) {
        super(new BucketState(N, M), BucketEvent.getAllEvents(M));
        this.M = M;
        this.N = N;
        this.mu = mu;
        assert (mu.length == N);
        assert (mu[0].length == M);
    }

    /**
     * This method construct a BB where workers consistently dominate
     * each other.
     * @param N Number of workers
     * @param M Number of machines
     * @param processRates Processing rates for each machine.
     * @param velocities Relative worker velocities.
     */
    public BucketBrigades(int N, int M, double[] processRates,
            double[] velocities) {
        super(new BucketState(N, M), BucketEvent.getAllEvents(M));
        this.M = M;
        this.N = N;
        mu = new double[N][M];
        for (int i = 0; i < N; i++)
            for (int j = 0; j < M; j++)
                mu[i][j] = processRates[j] * velocities[i];
        assert (mu.length == N);
        assert (mu[0].length == M);
    }

    /**
     * This boolean function returns true if one event e is active
     * (may occur) when the system is in a state i, and 0 otherwise
     */
    @Override
    public boolean active(BucketState i, BucketEvent e) {
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
    public States<BucketState> dests(BucketState i, BucketEvent e) {
        int[] newMachines = i.getMachines(); // The initial MaqOp
        // vector is
        // copied.
        int[] newWorkersStatus = new int[N];
        int[] newMachineStatus = new int[M];
        int eventWorker = getWorkerEvent(i, e.getMachine());
        // If the event is not a reset:Worker goes to the next machine
        if (e.getMachine() < (M - 1)) {
            newMachines[eventWorker] = e.getMachine() + 1;
        }

        // If the event is a reset
        if (e.getMachine() == (M - 1)) {
            // Every worker, but the first one, takes the work of the
            // previous
            // one
            for (int k = N - 1; k > 0; k--) {
                newMachines[k] = newMachines[k - 1];
            }
            // The first worker starts on the first machine
            newMachines[0] = 0;
        }
        // The vectors StatusOp and StatusMaq for the new state are
        // calculated
        newWorkersStatus = getWorkersStatus(newMachines);
        newMachineStatus = gsetMachineStatus(newMachines);
        return new StatesSet<BucketState>(new BucketState(newMachines,
                newWorkersStatus, newMachineStatus));
    }

    /**
     * Now the rate at which an event occurs, given a state, is
     * defined.
     */
    @Override
    public double rate(BucketState i, BucketState j, BucketEvent e) {
        return mu[getWorkerEvent(i, e.getMachine())][e.getMachine()];
    }

    /**
     * int[] SetStatusOp returns a vector with the status of all
     * workers for a new state, given an event e and a state s
     * @param machines
     * @return The operators status, consistent with the machines they
     *         are at.
     */

    private int[] getWorkersStatus(int[] machines) {
        int[] newWorkerStatus = new int[N];
        newWorkerStatus[N - 1] = 1; // The last worker is always
        // active
        for (int i = N - 2; i >= 0; i--) {
            if (machines[i] != machines[i + 1]) {
                newWorkerStatus[i] = 1;
            }
        }

        return newWorkerStatus;
    }

    /**
     * int[] setMachineStatus returns a vector with the status of all
     * machines for a new state, given an event e and a state s
     * @param machines
     * @return a vector with the status of all machines.
     */

    private int[] gsetMachineStatus(int[] machines) {
        int[] newMachinesStatus = new int[M];
        for (int i = 0; i < N; i++) {
            newMachinesStatus[machines[i]] = 1;
        }

        return newMachinesStatus;
    }

    /**
     * getWorkerEvent returns the index of the worker that has
     * finished on machine ev.machine. Since this is a common
     * BucketBrigades system, the worker will be the one with highest
     * number that was active on ev.machine
     */

    int getWorkerEvent(BucketState s, int m) {
        int result = 0;
        for (int i = 0; i < N; i++) {
            if (s.getMachine(i) == m)
                result = i;
        }
        return result;
    }

    /**
     * This is just s test program.
     * @param args Not used.
     */
    public static void main(String args[]) {
        String fileMu = "./BucketFiles/Mu.txt";

        double mu[][] = loadRateMatrix(fileMu);
        Matrix x = new Matrix(mu);
        final int Op = x.getRowDimension();
        final int M = x.getColumnDimension();
        // La siguiente instrucción crea el estado i0
        BucketBrigades theModel = new BucketBrigades(Op, M, mu);
        theModel.setDebugLevel(0);
        theModel.setMaxStates(2000);
        theModel.showGUI();
        theModel.generate();

        theModel.printAll();
        /*
         * System.out.println(theModel.MOPsToString());
         * System.out.println(theModel.statesToString());
         * System.out.println(theModel.eventsRatesToString());
         */
    }

    /**
     * Reads the process rates for every worker on every machine.
     * @param fileName
     * @returns the read process rate matrix.
     */
    private static double[][] loadRateMatrix(String fileName) {
        double[][] mu = null;
        Matrix Mu;
        File file = new File(fileName);
        String fullName = fileName;
        try {
            fullName = file.getCanonicalPath();
            Mu = Matrix.read(new BufferedReader(new FileReader(file)));
            mu = Mu.getArrayCopy();
        } catch (Exception e) {
            System.out.println("Unable to read " + fullName);
            System.exit(0);
        }
        return mu;
    }

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

        return stg;
    }

    /**
     * @see jmarkov.MarkovProcess#label()
     */
    @Override
    public String label() {
        return super.label() + " BB Expo (" + N + " workers and " + M
                + " machines).";
    }
}

/**
 * The state is characterized by the machine for each worker, the
 * status for each worker and machine. (Actually the first array
 * contains all needed information).
 * @author Diego Rojas, German Riano. Universidad de los Andes. (C)
 *         2006
 */
class BucketState extends PropertiesState {
    private static int M; // Number of machines
    private static int N; // Number of workers

    /**
     * First constructor. Creates the first state i0
     * @param numWorkers
     * @param numMachines
     */

    public BucketState(int numWorkers, int numMachines) {
        super(2 * numWorkers + numMachines);
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
            // All machines but the first one are inactive
            prop[2 * numWorkers + i] = 0;
        }
        prop[2 * numWorkers] = 1; // Only the first machine is
        // active
        // Now the Measures of Performance are set.
    }

    /**
     * Second constructor. Creates all other states
     * @param machines
     * @param workerStatus
     * @param machineStatus
     */

    BucketState(int[] machines, int[] workerStatus, int[] machineStatus) {
        super(2 * (machines.length) + machineStatus.length);

        System.arraycopy(machines, 0, prop, 0, N);
        System.arraycopy(workerStatus, 0, prop, N, N);
        System.arraycopy(machineStatus, 0, prop, 2 * N, M);
    }

    /** Sets all measures of performance */

    @Override
    public void computeMOPs(MarkovProcess mp) {
        int sumWk = 0;

        /*
         * We want to know how much time every worker spends on every
         * machine For this purpose, the function int askMaqOp is
         * created. It returns 1 if the worker i is currently
         * stationed at the machine j, and 0 if it is not
         */

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++)
                setMOP(mp, "Machine " + (j + 1) + " used by worker " + (i + 1),
                        askMachWorker(i, j));
        }
        // Now, the number of active workers in the system is
        // calculated.
        // Those workers which are not active, are assumed to be on
        // hold.
        for (int i = 0; i < N; i++) {
            sumWk += getWorkerStatus(i);
            setMOP(mp, "Status Worker " + (i + 1), getWorkerStatus(i));
        }
        // The average machine utilisation is also calculated
        for (int i = 0; i < M; i++)
            setMOP(mp, "Status Machine" + (i + 1), getMachineStatus(i));

        setMOP(mp, "Number of active workers in the system", sumWk);
    }

    /**
     * The following function returns 1 if the "worker" is stationed
     * at "machine" and 0 otherwise
     * @return 1 or 0
     */
    private int askMachWorker(int worker, int machine) {
        int result = 0;
        if (prop[worker] == machine) {
            result = 1;
        }
        return result;
    }

    /**
     * Returns the machine number the j-th worker is stationed at
     * @param j
     * @return the machine number
     */

    public int getMachine(int j) {
        return prop[j];
    }

    /**
     * Returns an integer-valued array with the machines where all
     * workers are.
     * @return an array with the machines.
     */
    public int[] getMachines() {
        int[] machines = new int[N];
        System.arraycopy(prop, 0, machines, 0, N);
        return machines;
    }

    /**
     * Returns 1 if the k-th worker is active and 0 otherwise
     * @param i
     * @return 1 if worker i is busy.
     */
    public int getWorkerStatus(int i) {
        return prop[i + N];
    }

    /**
     * Returns 1 if the k-th machine is active and 0 otherwise
     * @param j Machine number
     * @return 1 if busy
     */
    public int getMachineStatus(int j) {
        return prop[j + 2 * N];
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
        return stg;
    }

    /** Provides a description of each state */
    @Override
    public String description() {
        String stg = "Macines: ";
        for (int i = 0; i < N - 1; i++)
            stg += (getMachine(i)) + ",";
        stg += (getMachine(N - 1)) + ". Act.Wkrs. (";
        for (int i = 0; i < N; i++)
            stg += (getWorkerStatus(i) == 1) ? (i) + "" : "";
        stg += ") Wkrs hold. (";
        for (int i = 0; i < N; i++)
            stg += (getWorkerStatus(i) == 0) ? (i) + "" : "";
        stg += ") Act. M. (";
        for (int j = 0; j < M; j++)
            stg += (getMachineStatus(j) == 1) ? (j) + "" : "";
        stg += "} M hold {";
        for (int j = 0; j < M; j++)
            stg += (getMachineStatus(j) == 0) ? (j) + "" : "";
        stg += "}";
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
 * An event occurs when one worker finishes its work on its current
 * machine. In this way, the events can be described by one integer:
 * the machine number where the work has finished. If it is the last
 * one, this is a reset.
 */

class BucketEvent extends Event {
    // Revisar que el valor del número de máquinas esté correcto.
    static int M;// Number of machines.
    private int machine;

    /**
     * Constructor
     * @param machine The machine whwere work is completed.
     */
    public BucketEvent(int machine) {
        this.machine = machine;
    }

    /** A set with all posible events is created */
    static EventsSet<BucketEvent> getAllEvents(int M) {
        EventsSet<BucketEvent> E = new EventsSet<BucketEvent>();
        for (int i = 0; i < M; i++) {
            E.add(new BucketEvent(i));
        }
        BucketEvent.M = M;
        return E;
    }

    /** This function returns a description of each event */
    @Override
    public String label() {
        return (machine == (M - 1)) ? "Reset" : "Process end on machine "
                + (machine + 1);
    }

    /**
     * @return Returns the machine.
     */
    public final int getMachine() {
        return machine;
    }
}

// FALTAN LOS printMeasures() DEFINIDOS EN JACKSON.JAVA??
