package examples.jmarkov;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.URL;

import jmarkov.MarkovProcess;
import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.NotUnichainException;
import Jama.Matrix;

/**
 * This class represents a "Closed Jackson Network" with numStations multi-server
 * workstations and N entities in the system. Every entity goes from workstation
 * "i" to workstation "j" with probability "Pij".
 * 
 * @author Rodrigo Alberto Cáliz, Luis Felipe Parra (2003) Mod: Germán Riaño (2004-2006)
 */

public class Jackson extends SimpleMarkovProcess<JacksonState, JacksonEvent> {
    private double[][] probs;
    private double[] mus;
    private int[] numServers;
    private int[] wip;
    private int numStations;

    /**
     * Constructor using file names
     * 
     * @param filewip
     *            File with Initial WIP
     * @param fileservers
     *            File with number of servers
     * @param fileservicesrates
     *            File with service rates
     * @param fileprob
     *            File with transfer probability rates
     */
    public Jackson(String filewip, String fileservers,
            String fileservicesrates, String fileprob) {
        this(loadInt(filewip), loadInt(fileservers),
                loadDouble(fileservicesrates), loadMatrix(fileprob));
    }

    /**
     * Default Constructor used by GUI
     */
    public Jackson() {
        this("JacksonFiles/Wiplevel.txt", "JacksonFiles/Numservers.txt",
                "JacksonFiles/Servicesrates.txt",
                "JacksonFiles/Probabilities.txt");
    }

    /**
     * General constructor.
     * 
     * @param wip
     *            Initial wip array.
     * @param servers
     *            Number of servers per station
     * @param srvRates
     *            Service rates
     * @param probabilities
     *            Transfer probabilities
     */
    public Jackson(int[] wip, int[] servers, double[] srvRates,
            double[][] probabilities) {
        super(new JacksonState(wip), JacksonEvent.getAllEvents(wip.length,
                probabilities));
        numStations = wip.length;
        this.wip = wip;
        numServers = servers;
        mus = srvRates;
        probs = probabilities;
    }

    /**
     * Constructs a flow CONWIP line with M stations, and the given snumber of
     * servers and service rates.
     * 
     * @param wip
     *            Total wip
     * @param M
     *            Number of stations
     * @param servers
     *            Number of servers at each station.
     * @param srvRates
     *            M-dimensional array with the rates.
     */
    public Jackson(int wip, int M, int servers[], double[] srvRates) {
        this(initWip(wip, M), servers, srvRates, flowProb(M));
    }

    /**
     * Constructs a flow CONWIP line with M stations, single servers and the
     * given service rates.
     * 
     * @param wip
     *            Total wip
     * @param M
     *            Number of stations
     * @param srvRates
     *            M-dimensional array with the rates.
     */
    public Jackson(int wip, int M, double[] srvRates) {
        this(initWip(wip, M), initSingleServers(M), srvRates, flowProb(M));
    }

    /**
     * Constructs a flow CONWIP line with M stations, single servers and the
     * given service rates.
     * 
     * @param wip
     *            Total wip
     * @param M
     *            Number of stations
     * @param srvRates
     *            M-dimensional array with the rates.
     * @param prob
     *            Transfer probabilities.
     */
    public Jackson(int wip, int M, double[] srvRates, double[][] prob) {
        this(initWip(wip, M), initSingleServers(M), srvRates, prob);
    }

    /**
     * Construct an initial wip array with all wip on first station.
     * 
     * @param wip
     *            Total wip
     * @param M
     *            Number of stations.
     * @return Aforementioned array.
     */
    private static int[] initWip(int wip, int M) {
        int wipp[] = new int[M];
        wipp[0] = wip;
        return wipp;
    }

    /**
     * Used to build a single server array.
     * 
     * @param M
     *            Number of stations
     * @return (1,1,1...,1)
     */
    private static int[] initSingleServers(int M) {
        int srv[] = new int[M];
        for (int i = 0; i < M; i++) {
            srv[i] = 1;
        }
        return srv;
    }

    /**
     * Creates a Transfer probability corresponfing to a conwip flow line. i.e.
     * the items go to station i+1 with probability 1, and from the last station
     * a signal authorizes the new material in the first.
     * 
     * @param M
     *            Stations.
     * @return Aforementioned matrix.
     */
    private static double[][] flowProb(int M) {
        double[][] prob = new double[M][M];
        for (int i = 0; i < M - 1; i++) {
            prob[i][i + 1] = 1.0;
        }
        prob[M - 1][0] = 1.0;
        return prob;
    }

    /**
     * Loads a JAMA matrix from a file. The format is a text file, entries
     * separated by spaces, and each row of the matrix i a row of the file.
     * 
     * @param fileName
     *            The name fo the file.
     * @return jama.Matrix object read from file.
     */
    public static Matrix loadJamaMatrix(String fileName) {
        Matrix Mat = null;
        File fil = null;
        try {
            URL url = Jackson.class.getResource(fileName);
            fil = new File(url.toURI());
            System.out.println("Reading file: " + fil.getAbsolutePath());
            Mat = Matrix.read(new BufferedReader(new FileReader(fil)));
        } catch (Exception e) {
            System.out.println("Unable to read the file "
                    + fil.getAbsolutePath() + "\n");
            e.printStackTrace();
            System.exit(0);
        }
        return Mat;

    }

    /**
     * Load an array of int from a file.
     * 
     * @param fileName
     *            The name fo the file.
     * @return the read array.
     */
    public static int[] loadInt(String fileName) {
        int[] temp2 = null;
        Matrix C = loadJamaMatrix(fileName);
        double[] temp1 = C.getRowPackedCopy();
        int N = temp1.length;
        temp2 = new int[N];
        for (int i = 0; i < temp1.length; i++)
            temp2[i] = (int) temp1[i];
        return temp2;
    }

    /**
     * Reads an array of double form a file.
     * 
     * @param fileName
     *            The name fo the file.
     * @return the read array.
     */
    public static double[] loadDouble(String fileName) {
        Matrix C = loadJamaMatrix(fileName);
        return C.getRowPackedCopy();
    }

    /**
     * Load a matrix as an aary.
     * 
     * @param fileName
     *            The name fo the file.
     * @return the 2-D array read.
     */
    public static double[][] loadMatrix(String fileName) {
        Matrix C = loadJamaMatrix(fileName);
        return C.getArrayCopy();
    }

    /**
     * Determine the active events.
     */

    @Override
    public boolean active(JacksonState i, JacksonEvent e) {
        boolean result = false;
        int origin = e.getOrigin();
        int dest = e.getDest();
        result = (i.getWip(origin) > 0 && probs[origin][dest] > 0);
        return result;
    }

    @Override
    public States<JacksonState> dests(JacksonState i, JacksonEvent e) {
        int origin = e.getOrigin();
        int dest = e.getDest();
        return new StatesSet<JacksonState>(i.Trans(origin, dest));
    }

    /**
     * Returns the transition rate from State i to State j.
     */

    @Override
    public double rate(JacksonState i, JacksonState j, JacksonEvent e) {

        int origin = e.getOrigin();
        int dest = e.getDest();
        return mus[origin] * Math.min(i.getWip(origin), numServers[origin])
                * probs[origin][dest];
    }

    @Override
    public String description() {
        int totWIP = 0;
        for (int i = 0; i < numStations; i++)
            totWIP += wip[i];
        String stg = "CLOSED JACKSON NETWORK \n\nThere are " + totWIP
                + " units in the system\n\n";
        for (int i = 0; i < numStations; i++) {
            stg += "WorkStation " + (i + 1) + ":\n " + numServers[i]
                    + " Exponential servers with service rate " + mus[i] + "\n";
        }
        stg += "Transition Probabilities:\n";
        Matrix probMat = new Matrix(probs);
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw, true);
        probMat.print(pw, 7, 4);// 4 decimals, 7 width
        stg += sw.toString();
        return stg;
    }

    /**
     * Returns the throughput rate for every Station
     * 
     * @return An array with the rates.
     * @throws NotUnichainException
     */

    public double[] effLambdas() throws NotUnichainException {

        double sums[] = new double[numStations];
        Event[] events = getEvents();
        for (Event ev : events) {
            JacksonEvent e = (JacksonEvent) ev;
            sums[e.getOrigin()] += getEventRate(e.getIndex());
        }
        return sums;
    }

    /**
     * Replaces the method to printo MOPs. This method is called in
     * <code>toString()</code>, which is also used in the MOP's tab in the
     * GUI.
     */
    @Override
    public int printMOPs(PrintWriter out, int width, int decimals) {
        // super.printMOPs(out); // /optional: call default method
        double lambda[];
        try {
            lambda = effLambdas();
            for (int i = 0; i < numStations; i++) {
                out.println("STATION " + (i + 1) + ":");
                double wip = getMOPsAvg("Average WIP at the buffer " + (i + 1));
                double busyServers = getMOPsAvg("Busy Servers at Station "
                        + (i + 1));
                double numInQ = wip - busyServers;
                double queueTime = numInQ / lambda[i];
                double flowTime = wip / lambda[i];
                double rho = busyServers / this.numServers[i];
                out.printf("    Throughput        = %6.4f ", lambda[i]);
                out.printf("\n    Avg. WIP          = %6.4f ", wip);
                out.printf("\n    Avg. Flow Time    = %6.4f ", flowTime);
                out.printf("\n    Avg. Num in Queue = %6.4f ", numInQ);
                out.printf("\n    Avg. Queue Time   = %6.4f ", queueTime);
                out.printf("\n    Avg. Busy Servers = %6.4f ", busyServers);
                out.printf("\n    Avg. Utilization  = %6.4f ", rho);
                out.println();
            }
        } catch (NotUnichainException e) {
            out.println(e);
        }

        return 0;
    } // Main

    /**
     * This is just a test program.
     * 
     * @param a
     *            If given these are the files to read the data from.
     */
    public static void main(String a[]) {
        String file1 = "";
        String file2 = "";
        String file3 = "";
        String file4 = "";
        if (a.length < 4) {
            file1 = "JacksonFiles/Wiplevel.txt";
            file2 = "JacksonFiles/Numservers.txt";
            file3 = "JacksonFiles/Servicesrates.txt";
            file4 = "JacksonFiles/Probabilities.txt";
        } else {
            file1 = a[0];
            file2 = a[1];
            file3 = a[2];
            file4 = a[3];
        }
        Jackson theNet = new Jackson(file1, file2, file3, file4);
        theNet.setDebugLevel(1);
        theNet.showGUI();
        theNet.printAll();
    }

    /**
     * @return Returns the numServers.
     */
    public int[] getNumServers() {
        return numServers;
    }

    /**
     * @return Returns the numStations.
     */
    public int getNumStations() {
        return numStations;
    }

} // End of Jackson class

/**
 * A State in a Closed Jackson Network. This class represents a State in a
 * "Closed Jackson Network" with numStations multiserver workstations and N
 * entities in the system. Every entity goes from workstation "i" to workstation
 * "j" with probability "Pij". The State consist of a numStations dimensional
 * vector. Every position represents the number of entities at every
 * workstation.
 */

class JacksonState extends PropertiesState {

    /**
     * Constructor for the JacksonState class
     * 
     * @param wip
     *            number of entities at each workstation
     */
    public JacksonState(int[] wip) {
        super(wip);
    }

    /**
     * @see State#computeMOPs(MarkovProcess)
     */
    @Override
    public void computeMOPs(MarkovProcess mp) {
        int M = getNumStations();
        Jackson jacky = (Jackson) mp;
        int servers[] = jacky.getNumServers();

        for (int i = 0; i < M; i++) {
            setMOP(mp, "Busy Servers at Station " + (i + 1), Math.min(prop[i],
                    servers[i]));
        }
        for (int i = 0; i < M; i++) {
            setMOP(mp, "Average WIP at the buffer " + (i + 1), prop[i]);
        }
    }

    /**
     * Returns the number of entities in a given workstation.
     * 
     * @param s
     *            Workstation index
     * @return Number of entities at the specified workstation
     */
    public int getWip(int s) {
        return prop[s];
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
     * @return The number of workstations.
     */
    public int getNumStations() {
        return prop.length;
    }

    /**
     * Creates a new state obtained by moving an entity from workstation i to j.
     * 
     * @param i
     *            Origin station
     * @param j
     *            Destination station
     * @return The new state
     */
    public JacksonState Trans(int i, int j) {
        int[] newWip = getProperties();
        if (i != j) {
            newWip[i] = newWip[i] - 1;
            newWip[j] = newWip[j] + 1;
        }
        return new JacksonState(newWip);
    }

    /**
     * Describes the States
     * 
     * @return a String description of the State.
     */
    @Override
    public String description() {
        int M = getNumStations();

        String stg = "";
        for (int i = 0; i < M; i++) {
            stg += prop[i] + " at Station " + (i + 1)
                    + ((i < M - 1) ? ", " : "");
        }
        return stg;
    }

} // End of JacksonState class.

/**
 * Each Event characterizes the movement from i to j.
 * 
 * @author Rodrigo Alberto Caliz, Luis Felipe Parra, Germán Riaño
 */
class JacksonEvent extends Event {

    private int origin;

    private int dest;

    /**
     * Constructor of JacksonEvent
     * 
     * @param origin
     *            Origin state index
     * @param dest
     *            Destination State index
     */
    public JacksonEvent(int origin, int dest) {
        this.origin = origin;
        this.dest = dest;
    }

    /**
     * Returns the destination state of this event
     * 
     * @return Destination state index of this event
     */
    public int getDest() {
        return dest;
    }

    /**
     * Returns the origin state of this event
     * 
     * @return Origin state index of this event
     */
    public int getOrigin() {
        return origin;
    }

    /**
     * Returns all the possible events in a given state
     * 
     * @param M
     *            State index
     * @param prob
     *            transition probability between workstations
     * @return All the possible events in a given state
     */
    public static EventsSet<JacksonEvent> getAllEvents(int M, double[][] prob) {
        EventsSet<JacksonEvent> theSet = new EventsSet<JacksonEvent>();
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                if (prob[i][j] > 0)
                    theSet.add(new JacksonEvent(i, j));
            }
        }
        return theSet;
    }

    /**
     * @see java.lang.Object#toString()
     */
    @Override
    public String label() {
        return "From " + (origin + 1) + " to " + (dest + 1);
    }

}