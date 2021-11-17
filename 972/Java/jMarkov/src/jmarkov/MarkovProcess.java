package jmarkov;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.JMarkovElement;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.Transition;
import jmarkov.basic.Transitions;
import jmarkov.basic.TransitionsSet;
import jmarkov.basic.exceptions.NotUnichainException;
import jmarkov.gui.MarkovGUI;
import jmarkov.solvers.JamaTransientSolver;
import jmarkov.solvers.MtjSolver;
import jmarkov.solvers.SteadyStateSolver;
import jmarkov.solvers.TransientSolver;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;

/**
 * The abstract class SimpleMarkovProcess represents a Continuous or
 * Discrete Time Markov Chain. In order to model a particular problem
 * the user has to extend this class. The class can generate the model
 * through the buildRS algorithm. This enables it to generate all
 * states and the transition matrix, from behavior rule given by the
 * user. These rules are determined by implementing the methods,
 * <code>active</code>,<code>dests</code> and <code>rate</code>.
 * The user should also determine how to code the space state. This is
 * accomplished by implementing the State class. A particular
 * implementation of State is provided where each state is coded with
 * k integer properties. Examples are included in this release.
 * @see SimpleMarkovProcess#dests(State, Event)
 * @see SimpleMarkovProcess#active(State, Event)
 * @see SimpleMarkovProcess#rate(State, State, Event)
 * @see State
 * @see Event
 * @see PropertiesState
 * @author German Riaño. Universidad de los Andes.
 * @version 1.0a
 * @param <S> The State class for the model.
 * @param <E> The Event class for the model.
 */
public abstract class MarkovProcess<S extends State, E extends Event>
        implements JMarkovElement {

    /** Initial State */
    private S i0; // Initial state;
    /**
     * Number of completed states. Use this counter so that GUI
     * updates correctly.
     */
    protected int cnt = 0;

    /** Set of fully analyzed States */
    protected StatesSet<S> theStates = new StatesSet<S>();
    SortedSet<S> uncheckedStates = new TreeSet<S>(); // Unchecked
    /** states found so far */
    private SortedMap<S, S> foundStates = new TreeMap<S, S>();
    /** Number of events */
    private SortedMap<S, Transitions<S>> rates = new TreeMap<S, Transitions<S>>();

    private int numE; // Number of Events

    boolean generated = false; // tells whether the model has been
    // generated.

    /**
     * Maximum number of States to Generate. It may be limited by your
     * license
     */
    private long maxStates = 1000;
    private double[][] ratesMatrix = null; // Array to hold dense
    // matrix.
    private double[][] genMatrix = null; // Array to hold dense
    // matrix.
    private double[] thePi = null; // SS probabilities

    /**
     * The name of the model. For a long description override
     * <code>description()<code>.
     * @see #description()
     */
    protected String name = "";
    // private static ArrayList mopsNames = new ArrayList();

    // StatesSet<S> theStates = null; // Array to hold the list of
    // States.
    E[] theEvents = null; // array to hold the Events.

    DebugReporter reporter = new DebugReporter(1);

    // SteadyStateSolver steadyStateSolver = new JamaSolver(this);
    SteadyStateSolver steadyStateSolver = null;
    TransientSolver transientSolver = null;

    /** Default Transient Solver */
    protected SteadyStateSolver defaultSteadyStateSolver = null;
    /** Deafualt Transient solver */
    protected TransientSolver defaultTransientSolver = null;

    /** Status variables */
    public enum Status {
        /** Model has not been generated */
        IDLE,
        /** Model is being generated */
        RUNNING,
        /** Model has been generated */
        GENERATED,
        /** Model execution has been suspended */
        SUSPENDED,
        /** Model generated, writing matrix */
        WRITING,
        /** Model execution generated an error */
        ERROR,
        /** No Model loaded. (used by GUI) */
        NoModel
    };

    // model status
    private Status _status = Status.IDLE;

    // This code romoved july 01 2005 by GRM. With generics it makes
    // no sense
    // since the class does need an Event class E.
    // /**
    // * Builds a SimpleMarkovProcess that contains all states
    // reachable from
    // i0, and
    // * with numE possible events. The Events are numbered from 0 to
    // numE-1.
    // *
    // * @param i0
    // * The initial State.
    // */
    // public SimpleMarkovProcess(S i0, int numE) {
    // this(i0, EventsSet.createEventSet(numE));
    // }

    /**
     * Builds a SimpleMarkovProcess that contains all states reachable
     * from i0, and with E being the set of all possible events.
     * @param i0 The initial State.
     * @param eSet The set of all Events.
     * @param name The name of the Model.
     */
    public MarkovProcess(S i0, EventsSet<E> eSet, String name) {
        // this.i0 = i0;
        // uncheckedStates.add(i0);
        // theEvents = eSet.toEventArray();
        clearMOPs();
        setInitialState(i0);
        setEventSet(eSet);
        this.name = name;
    }

    /**
     * Builds a SimpleMarkovProcess that contains all states reachable
     * from i0, and with E being the set of all possible events.
     * @param i0 The initial State.
     * @param eSet The set of all Events.
     */
    public MarkovProcess(S i0, EventsSet<E> eSet) {
        this(i0, eSet, "");
    }

    /**
     * If a constructor calls this constructor then it MUST call
     * setEvents and setInitialState afterwards.
     */
    protected MarkovProcess() {
        reset(false);

    };

    /**
     * Sets the Events set. It causes the model to be reset. This
     * method should be called only in a constructor.
     * @param eSet the Events set.
     */
    protected void setEventSet(EventsSet<E> eSet) {
        theEvents = eSet.toEventArray();
        reset(false);
    }

    /**
     * Sets the initial state. It casues the model to be reset. This
     * method should be called only in a constructor.
     * @param i0 The initial state.
     */
    protected void setInitialState(S i0) {
        this.i0 = i0;
        uncheckedStates = new TreeSet<S>(); // Estados no revisados
        clearMOPs();
        if (i0 != null)
            uncheckedStates.add(i0);
        reset(false);
    }

    /**
     * returns the initial state.
     * @return The initial state.
     */
    protected S getInitialState() {
        return i0;
    }

    /**
     * Resets the Model. It erases all found states and transition
     * rates. Keeps the initial state and Events set.
     */
    public synchronized void reset() {
        reset(false);
    }

    /**
     * Resets the Model. It erases all found states and transition
     * rates.
     * @param resetEvents whether the Events are deleted. WARNNING! if
     *        this is true you must call setEventSet.
     * @see #setInitialState(State)
     * @see #setEventSet(EventsSet)
     */
    protected synchronized void reset(boolean resetEvents) {

        uncheckedStates = new TreeSet<S>(); // Estados no revisados
        clearMOPs();
        if (i0 != null)
            uncheckedStates.add(i0);
        if (resetEvents)
            theEvents = null;
        theStates = null;
        generated = false;
        ratesMatrix = null;
        genMatrix = null;
        foundStates = new TreeMap<S, S>(); // All the states
        // Encontrados
        theStates = new StatesSet<S>(); // Explored states
        rates = new TreeMap<S, Transitions<S>>();
        resetResults();
        setStatus(Status.IDLE);
    }

    /**
     * Resets the result of the model. If it has been generated this
     * method keeps the Graph, but erases steady state, and transient
     * probabilities.
     */
    public synchronized void resetResults() {
        thePi = null;
        cnt = 0;
        setStatus(Status.IDLE);
    }

    // ******************************************
    // ABSTRACT FUNCTIONS
    // ******************************************

    /**
     * The user MUST implement this Function in order to describe the
     * dynamics of the model. For the current state i, and on action
     * e, the user has to describe the transtions that can occur. This
     * implies finding all destination states and the rate at which
     * the transtions occur. There is no guarantee that the event is
     * active, so the user should check for this. If the event is not
     * active an empty Transition element should be returned. A
     * typical code for a queuing system should look like this: <code>
     * public abstract Transitions<MyState> activeTransitions(MyState i, MyEvent e){
     *      TransitionsSet<MyState> trans = new TransitionsSet<MyState>();
     *      case (ARRIVAL)
     *          if (i.size() < capacity)
     *              trans.add(i.doArrival(), arrRate);
     *          break;
     *      case(DEPARTURE)
     *          if (i.size() >=1)
     *              trans.add(i.doDeparture, serviceRate);
     *          break;
     *      }
     *      return trans;
     * }
     * </code>
     * @param i The current State.
     * @param e The ocurring event.
     * @return The transitions that occur at this state when (and if)
     *         this events occurs.
     * @see Transitions
     * @see TransitionsSet
     * @see Transitions
     */
    public abstract Transitions<S> activeTransitions(S i, E e);

    /**
     * generate() builds the space state and rate matrix using the
     * algorithm BuildsSR. (See Ciardo, G. "Tools for formulating
     * Markov Processes", chapter 2 in Grassman W. "Computational
     * Probability". Kluwer). The states can be collected later with
     * getStates and the rates can be accessed in disperse form with
     * the method <code>getRate(j)</code> of every state.
     * Alternatively the method <code>getGenerator()</code> and
     * <code>getRates()</code> access the generator matriz or rate
     * matrix in compact form.
     */

    public void generate() {
        if (getStatus() == Status.IDLE) {
            debug(1, "MODEL: " + label());
            if (getDebugLevel() >= 2) {
                EventsSet<E> eSet = new EventsSet<E>(theEvents);
                debug(2, "Events Set: " + eSet);
            }
            debug(1, "Generating model ...");
        } else
            debug(1, "Continuing model generation ...");
        setStatus(Status.RUNNING);
        // / main cycle of SR Algorithm:
        while (!uncheckedStates.isEmpty() && (cnt < maxStates) && canGo()) {
            generateStep();
        }
        setEndStatus();
    }

    /**
     * Sets the status at the end of generate() or generateStep().
     */
    private void setEndStatus() {
        // generated = canGo() && (cnt < maxStates);
        generated = uncheckedStates.isEmpty();
        debug(1, ""); // newline
        if (generated) {
            debug(1, "DONE. " + cnt + " states found.");
            theStates.numerateStates();
            // new line if debug >=1
            setStatus(Status.GENERATED);
        } else if (cnt == maxStates) {
            setStatus(Status.ERROR);
            debug(0, "Maximum Number of States Reached! (" + cnt + ")");
            throw new RuntimeException("Maximum Number of States Reached! ("
                    + cnt + ")");
        } else {
            debug(0, "PROCESS INTERRUPTED");
            setStatus(Status.SUSPENDED);
        }
    }

    /**
     * Updates the current rate thru R(i,j) := R(i,j) + rate(e,i)
     * @param i current State.
     * @param j destination State.
     * @param e event.
     */
    private void updateRates(S i, Transitions<S> currentRates, S j, double rate) {
        double curVal, newVal;
        curVal = currentRates.addRate(j, rate);
        newVal = curVal + rate;
        debug(5, "Rate(" + i + "," + j + ") = " + newVal + " (It was " + curVal
                + ")");
    }

    /**
     * Generates a single state step. (i.e. completely explores an
     * state). It assumes that there are indeed states needing to be
     * explored
     */
    private void generateStep() {
        S i;
        S curj;
        i = uncheckedStates.first();
        if (reporter.getDebugLevel() == 1) {
            if (cnt % 100 == 0) {
                debug(1, "", true); // newline
                debug(1, "", false, true); // indent
            }
            debug(1, ".", false);
        }
        cnt++;
        debug(2, "Building State " + i);
        debug(5, "S = " + theStates);
        debug(5, "U = " + uncheckedStates);
        uncheckedStates.remove(i);
        theStates.add(i);
        foundStates.put(i, i);
        i.computeMOPs(this);
        Transitions<S> totalRates = new TransitionsSet<S>();
        rates.put(i, totalRates);
        for (E e : theEvents) {
            Transitions<S> trans = activeTransitions(i, e);
            if (trans == null)
                break;
            for (Transition<S> tr : trans) {
                Thread.yield();
                S j = tr.getState();
                curj = foundStates.get(j);
                if (curj == null) { // not in the set
                    debug(5, "New state found j = " + j);
                    uncheckedStates.add(j);
                    foundStates.put(j, j);
                    // j.computeMOPs();
                    assert (j.isConsistent());
                } else { // it is already in the set
                    j = curj; // we ensure it is the same
                }
                double rate = tr.getRate();
                debug(4, "Event [" + e + "]: Transition " + i + " ---> " + j
                        + " with rate " + rate);
                updateRates(i, totalRates, j, rate);
            } // end for
        } // end for
        Thread.yield();
    }

    /**
     * Return the number of States in the model. If the model has not
     * been generated, then it will be automatically generated.
     * @return the number of States in the model.
     */

    public int getNumStates() {
        States<S> stts = getStates();
        return stts.size();
    }

    /**
     * Returns an array with all the States in the model. It generates
     * the model if it has not been generated.
     * @return The States
     */

    public StatesSet<S> getStates() {
        return getStates(true);
    }

    /**
     * Returns an array with the States in the model that have been
     * checked so far. If <code>generate</code> is true it generates
     * the model if it has not been generated. If no states have been
     * generates it returns null.
     * @param causesGeneration whether the model should be generated.
     * @return ann array with the states found and checked so far.
     */
    public StatesSet<S> getStates(boolean causesGeneration) {

        if (!generated && causesGeneration)
            generate();
        if (theStates != null)
            return theStates;
        debug(1, "Building States...");
        theStates.numerateStates();
        return (theStates);
    }

    /**
     * Returns all The events defined in the model.
     * @return Events array
     */
    public E[] getEvents() {
        return theEvents;
    }

    /**
     * This method returns a dynamic data structure with the rate from
     * State i to all reachable states.
     * @param i State
     * @return Rates to reachable states
     * @see Transitions
     */
    public synchronized Transitions<S> getRates(S i) {
        return rates.get(i);
    }

    /**
     * Gets the current total rate form i to j. <b>Warnning:</B> If
     * the model has not been generated id returns the current total
     * rate, it does not causes the generation of the model. Run
     * <code>generate()</code> first.
     * @see #generate()
     * @param i origin State
     * @param j destination state
     * @return total rate
     */
    public synchronized double getRate(S i, S j) {
        if (!generated)
            generate();
        return rates.get(i).getRate(j);
    }

    /**
     * Gets the total rate form State number i to j. The number is
     * relative to the ordered set of all states, therefore this
     * method causes the model to be generated if it has not been
     * generated.
     * @see #getRate(State, State)
     * @param i origin State
     * @param j destination state
     * @return total rate
     */
    public synchronized double getFinalRate(S i, S j) {
        if (!generated)
            generate();
        return rates.get(i).getRate(j);
    }

    /**
     * Returns the transition rates matrix <b>R</b> in dense format.
     * It generates the model if it has not been generated.
     * @return an array with the matrix.
     * @see SimpleMarkovProcess#getRate(State, State)
     * @see SimpleMarkovProcess#getRates(State)
     * @see SimpleMarkovProcess#getMtjRates()
     */
    public synchronized double[][] getRates() {
        if (!generated)
            generate();
        if (ratesMatrix != null)
            return ratesMatrix;
        setStatus(Status.RUNNING);
        debug(1, "Building dense Rates Matrix ..");
        int n = getNumStates();
        ratesMatrix = new double[n][n];
        cnt = 0;
        for (Map.Entry<S, Transitions<S>> entry : rates.entrySet()) {
            debug(3, "Row " + cnt + ", State: " + entry.getKey());
            Transitions<S> trans = entry.getValue();
            for (Transition<S> tr : trans) {
                double valObj = tr.getRate();
                S s2 = tr.getState();
                int j = s2.getIndex();
                debug(4, "Col " + j + ", State: " + s2);
                ratesMatrix[cnt][j] = valObj;
                assert (cnt == entry.getKey().getIndex());
            }
            cnt++;
        }
        setStatus(Status.GENERATED);
        return ratesMatrix;
    }

    /** The Rates matrix in MTJ */
    private Matrix mtRatesMatrix = null;
    private PrintWriter out;

    /**
     * Returns the transition rates matrix <b>R</b> in MTJ format. It
     * generates the model if it has not been generated.
     * @return an mt.Matrix with the matrix.
     */
    public synchronized Matrix getMtjRates() {
        if (!generated)
            generate();
        if (mtRatesMatrix != null)
            return mtRatesMatrix;
        setStatus(Status.RUNNING);
        debug(1, "Building Rates Matrix ..");
        int n = getNumStates();
        int[] nz = new int[n];
        int i = 0;
        for (Transitions<S> trans : rates.values()) {
            nz[i] = trans.size();
            i++;
        }

        // mtRatesMatrix = new CompRowMatrix(n, n, nz);
        // TODO: es esto lo mejor??
        mtRatesMatrix = new FlexCompRowMatrix(n, n);

        i = 0;
        for (Map.Entry<S, Transitions<S>> entry : rates.entrySet()) {
            cnt = i;
            Transitions<S> trans = entry.getValue();
            debug(3, "Row " + i + ", State: " + entry.getKey());
            for (Transition<S> tr : trans) {
                S s = tr.getState();
                double rate = tr.getRate();
                int j = s.getIndex();
                debug(4, "Col " + j + ", State: " + s);
                mtRatesMatrix.set(i, j, rate);
            }
            i++;
        }
        debug(3, "The Rates matrix is: \n" + mtRatesMatrix);
        setStatus(Status.GENERATED);
        return mtRatesMatrix;
    }

    /**
     * @return true if the model has been completely generated
     */
    public boolean isGenerated() {
        return generated;
    }

    /**
     * Returns the infinitesimal generator matrix <b>Q </b>, in dense
     * format. It generates the model if it has not been generated.
     * @return The generator Matrix.
     */

    public double[][] getGenerator() {
        if (genMatrix != null)
            return genMatrix;
        debug(1, "Building Generator Matrix ...");
        double sum;
        double[][] R = getRates();
        int n = R.length;
        double[][] Q = new double[n][n];
        for (int i = 0; i < n; i++) {
            sum = 0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    Q[i][j] = R[i][j];
                    sum += R[i][j];
                }
            }
            Q[i][i] = -sum;
        }
        return (genMatrix = Q);
    }

    /**
     * The generator <B>Q</B> as an MTJ Matrix
     * @return The matrix Q.
     */
    public Matrix getMtjGenerator() {
        MtjSolver sol;
        if (steadyStateSolver instanceof MtjSolver) {
            sol = (MtjSolver) steadyStateSolver;
        } else {
            sol = new MtjSolver(this);
        }
        return sol.getGenerator();
    }

    /**
     * Returns the steady state probabilities for this model. Theat
     * is, it solves the balance equations. It returns an array of
     * zeros if there is no unique solution.
     * @return An array with the steady-state probabilities.
     * @throws NotUnichainException
     */
    public double[] getSteadyState() throws NotUnichainException {
        if (!generated)
            generate();
        if (thePi != null)
            return thePi;
        thePi = getSteadyStateSolver().getSteadyState();
        return thePi;
    }

    /**
     * The Class for the states in this model.
     * @return The class for the states
     */
    public Class getStateClass() {
        return i0.getClass();
    }

    /**
     * The Class for the Events in the system.
     * @return The event Class in the model
     */
    public Class getEventClass() {
        return (theEvents.length > 0) ? theEvents[0].getClass() : Event.class;
    }

    // ***************************************************************************
    // MOPS STUFF
    // ***************************************************************************

    private List<String> mopsNames = new ArrayList<String>();

    /**
     * This method declares the existance of a measure of performance
     * (MOP). The MOP for every state is calculated in the class that
     * extends the State class.
     * @param mopName The name of the new MOP.
     * @return true if the name already existed.
     */
    public boolean addMOP(String mopName) {
        return mopsNames.add(mopName);
    }

    /**
     * Clear all MOPs defined in the system.
     */
    public void clearMOPs() {
        mopsNames = new ArrayList<String>();
    }

    /**
     * Sets the names of all MOPs (measures of performance).
     * @param mopNames
     */
    public void setMOPs(String[] mopNames) {
        for (int i = 0; i < mopNames.length; i++)
            addMOP(mopNames[i]);
    }

    /**
     * Return all the names of defined MOPs.
     * @return an array with all the MOP's defined.
     */
    public String[] getMOPNames() {
        return mopsNames.toArray(new String[1]);
    }

    /**
     * Returns an array with the average of all the steady state
     * measures of performance. The order is the same as in
     * getMOPNames.
     * @return An array containing the values of all MOPs averages.
     * @throws NotUnichainException
     * @see #getMOPsMoment(int)
     */
    public double[] getMOPsAvg() throws NotUnichainException {
        return getMOPsMoment(1);
    }

    /**
     * Returns an array with the m-th moment of all the steady state
     * measures of performance. The order is the same as in
     * getMOPNames.
     * @param m the order of the moment desired. m=1 is the expected
     *        value.
     * @return An array containing the values of all MOPs m-th
     *         moments.
     * @throws NotUnichainException
     */
    public double[] getMOPsMoment(int m) throws NotUnichainException {
        debug(1, "Computing MOPs moment " + m);
        if (!generated)
            generate();
        String[] mopNames = getMOPNames();
        int M = mopNames.length;
        double[] result = new double[M];
        for (int idx = 0; idx < M; idx++) {
            result[idx] = getMOPsMoment(idx, m);
        }
        return result;
    }

    /**
     * Gets the index that correspond to this MOP.
     * @param name MOP name
     * @return The index of the MOP with this name
     */
    public int getMOPIndex(String name) {
        return mopsNames.indexOf(name);
    }

    /**
     * Returns the number of defined Measures of performance (MOPs).
     * @return the number of MOPs defined so far.
     */
    public int numMOPs() {
        return mopsNames.size();
    }

    /**
     * Returns the steady state measures average of the MOP numbre
     * mopNum.
     * @param mopNum The Number of the MOP of the which the average is
     *        to be computed.
     * @return The long run averagefor this MOP.
     * @throws NotUnichainException
     */
    public double getMOPsAvg(int mopNum) throws NotUnichainException {
        return getMOPsMoment(mopNum, 1);
    }

    /**
     * Returns the steady state measures average of the MOP with name
     * mopName.
     * @param mopName The name whose Averga is to be computed.
     * @return The long run averagefor this MOP.
     * @throws NotUnichainException
     */
    public double getMOPsAvg(String mopName) throws NotUnichainException {
        return getMOPsMoment(mopName, 1);
    }

    /**
     * Returns the steady state measures m-th moment of the MOP number
     * mopNum. m=1 is the long-run expected value, m=2 expected value
     * of the square, etc.
     * @param mopNum The number for the MOP
     * @param m The value of m.
     * @return The m-th moments for this MOP.
     * @throws NotUnichainException
     */
    public double getMOPsMoment(int mopNum, int m) throws NotUnichainException {
        thePi = getSteadyState();
        theStates = getStates();
        String mopName = getMOPNames(mopNum);
        debug(1, "Computing moment " + m + " for MOP " + mopName);
        double sum = 0.0;
        int i = 0;
        for (S s : theStates) {
            double mopVal = s.getMOP(mopNum);
            if (m == 1)
                sum += thePi[i++] * mopVal;
            else if (m == 2)
                sum += thePi[i++] * mopVal * mopVal;
            else
                sum += thePi[i++] * Math.pow(mopVal, m);
        }
        return sum;
    }

    /**
     * Returns the steady state measures m-th moment of the MOP with
     * name mopName. m=1 is the long-run expected value.
     * @param mopName The name of the MOP that is to be computed
     * @param m Valu of the moment
     * @return The m-th moments for this MOP.
     * @throws NotUnichainException
     */
    public double getMOPsMoment(String mopName, int m)
            throws NotUnichainException {
        int idx = getMOPIndex(mopName);
        return getMOPsMoment(idx, m);
    }

    /**
     * Return the names of the i-th MOP.
     * @param mopNum The number i of the MOP
     * @return The name of the i-th MOP.
     */
    public String getMOPNames(int mopNum) {
        return (mopsNames.get(mopNum).toString());
    }

    /**
     * Return a String description of all MOPs in steady state (it
     * reports mean and standard deviation).
     * @return A String description of all MOPs.
     */
    public String MOPsToString() {
        StringWriter stw = new StringWriter();
        printMOPs(new PrintWriter(stw));
        return stw.toString();
    }

    /**
     * Return a String description of all MOPs in steady state (it
     * reports mean and standard deviation).
     * @param width the columns width
     * @param decimals the number of decimals to use.
     * @return String with a table representing MOPs names, means and
     *         standard deviations.
     */
    public String MOPsToString(int width, int decimals) {
        StringWriter stw = new StringWriter();
        printMOPs(new PrintWriter(stw), width, decimals);
        return stw.toString();
    }

    /**
     * Return a string as eventsRatesToString, with width 8 and 4
     * decimals
     * @return A string with the information.
     * @see #eventRatesToString(int,int)
     */
    public String eventsRatesToString() {
        return eventRatesToString(10, 5);
    }

    /**
     * Return a String as printed by printEventsrates
     * @param width Maximum width for each number
     * @param decimals Number of decimals
     * @return A string with the valus of all Rates
     * @see #printEventsRates(PrintWriter,int,int)
     */
    public String eventRatesToString(int width, int decimals) {
        StringWriter stw = new StringWriter();
        printEventsRates(new PrintWriter(stw), width, decimals);
        return stw.toString();
    }

    /**
     * Prints a table reporting the steadystate occurrance of all
     * events.
     * @param out where the table will be printed.
     */
    public void printEventsRates(PrintWriter out) {
        printEventsRates(out, 10, 5);
    }

    /**
     * Prints a table reporting the steadystate occurrance of all
     * events.
     * @param out where the table will be printed.
     * @param width The column width
     * @param decimals The number of decimals to use.
     */
    public void printEventsRates(PrintWriter out, int width, int decimals) {
        debug(1, "Printing Rates...");
        out.println("EVENTS OCCURANCE RATES");
        out.println();
        double[] means;
        try {
            means = getEventsRates();
            int M = means.length;
            String names[] = new String[M];
            int nameWidth = 0;
            for (int k = 0; k < M; k++) {
                names[k] = theEvents[k].toString();
                nameWidth = Math.max(nameWidth, names[k].length());
            }
            nameWidth = Math.max(nameWidth + 1, 7);
            out
                    .println(pad("NAME", nameWidth, false)
                            + pad("MEAN RATE", width));
            for (int k = 0; k < M; k++) {
                out.println(pad(names[k], nameWidth, false)
                        + pad(means[k], width, decimals));
            }
        } catch (NotUnichainException e) {
            out.println(e);
        }
        out.flush();
    }

    /**
     * Returns the defined events
     * @return aan array with the names of all Events
     */
    public String[] getEventNames() {
        int numE = theEvents.length;
        String[] names = new String[numE];
        for (int i = 1; i < numE; i++) {
            names[i] = theEvents[i].toString();
        }
        return names;
    }

    /**
     * Return the steadystate rate of occurrance of the Events number
     * eNum.
     * @param eNum The even number in the event set.
     * @return an array for all the steady state rates. The order is
     *         that of the events set.
     * @throws NotUnichainException
     */
    public double getEventRate(int eNum) throws NotUnichainException {
        E e;
        Transitions<S> trans;
        e = theEvents[eNum];
        debug(1, "Computing Events Rate for Event " + e);
        double sumRate = 0.0;
        double pi[] = getSteadyState();
        cnt = 0;
        for (S s : theStates) {
            trans = activeTransitions(s, theEvents[eNum]);
            for (Transition<S> tr : trans) {
                sumRate += pi[cnt] * tr.getRate();
            }
            cnt++;
        }
        return sumRate;
    }

    /**
     * Return an array with the steadystate rate of occurrance of all
     * the Events.
     * @return an array for all the steady state rates. The order is
     *         that of the events set.
     * @throws NotUnichainException
     */
    public double[] getEventsRates() throws NotUnichainException {
        debug(1, "Computing Events Rates");
        int M = theEvents.length;
        double[] result = new double[M];
        for (int k = 0; k < theEvents.length; k++) {
            result[k] = getEventRate(k);
        }
        return result;
    }

    /**
     * Prints the Measures of performance (MOPS) on standard output.
     * @see #printMOPs(PrintWriter)
     * @see #printMOPs(PrintWriter, int, int)
     */
    public final void printMOPs() {
        printMOPs(new PrintWriter(System.out, true));
    }

    /**
     * Prints a String description of all MOPs in steady state (it
     * reports mean and standard deviation), with a width of 10 and 5
     * decimal figures.
     * @param out The printer where the MOPS will be printed.
     * @see #printMOPs()
     * @see #printMOPs(PrintWriter, int, int)
     */

    public final void printMOPs(PrintWriter out) {
        printMOPs(out, 10, 5);
    }

    /**
     * Prints a String description of all MOPs in steady state (it
     * reports mean and standard deviation). You can override this
     * method to print your own MOPs. You can call it in the first
     * line like this <code>
     * public void printMOPs(PrintWriter out, int width, int decimals) {
     *      int namesWidth = super.printMOPs(out,width, decimals);
     *      // your oun code here:
     *      out.println(pad("Another MOP", namesWidth, false)
     *                  + pad(Value, width, decimals);
     *      }
     * </code>
     * @param out The printer where the MOPS will be printed.
     * @param width the columns width
     * @param decimals the number of decimals to use.
     * @return The max width among the declared MOPs. You can use this
     *         to align nicely your own MOPs.
     * @see #printMOPs(PrintWriter)
     * @see #printMOPs(PrintWriter, int, int)
     */

    public int printMOPs(PrintWriter out, int width, int decimals) {
        if (numMOPs() == 0)
            return 0;
        out.println("MEASURES OF PERFORMANCE");
        out.println();
        int namesWidth = 0;
        double[] means;
        try {
            means = getMOPsMoment(1);
            double[] mom2 = getMOPsMoment(2);
            String[] names = getMOPNames();
            int M = means.length;
            for (int k = 0; k < M; k++)
                namesWidth = Math.max(namesWidth, names[k].length());
            namesWidth++;
            out.println(pad("NAME", namesWidth, false) + pad("MEAN", width)
                    + pad("SDEV", width) + "\n");
            for (int k = 0; k < M; k++) {
                double mean = means[k];
                double dif = mom2[k] - mean * mean;
                if (dif < 0.0 && dif >= -1.0e-4)
                    dif = 0.0;// some solvers might give an error.
                assert (dif >= 0.0);
                double sdev = Math.sqrt(dif);
                out.println(pad(names[k], namesWidth, false)
                        + pad(means[k], width, decimals)
                        + pad(sdev, width, decimals));
            }
        } catch (NotUnichainException e) {
            out.print(e);
        }
        out.flush();
        return namesWidth;
    }

    @Override
    public String toString() {
        return label();
    }

    /**
     * Retuns a String description of the model and solution.
     * @return a String wit the information of printAll.
     * @see #printAll()
     */
    public String allToString() {
        StringWriter stw = new StringWriter();
        printAll(new PrintWriter(stw));
        return stw.toString();
    }

    /**
     * Prints a description of the Model: the States and the
     * Transition Matrix. If the model has less than 100 states it
     * shows all states and transition matrix and steady state
     * probabilities. Otherwise only the description, measures of
     * performance and events rates are shown.
     * @see #allToString()
     */
    public void printAll() {
        // creates a PrintWriter with auto-flush.
        PrintWriter pw = new PrintWriter(System.out, true);
        printAll(pw);
    }

    /**
     * Prints to the given PrintWriter a summary of the information
     * related to this MarkovChain. The information is the same as as
     * in the method <code>printAll()</code>.
     * @see #toString()
     * @see #printAll()
     * @param out
     */
    public void printAll(PrintWriter out) {
        int N = getNumStates();
        out.println(description() + "\n");
        out.println("System has " + N + " States.\n");
        if (N < 100) {
            printStates(out);
            printDenseMatrix(out);
        }
        out.println();
        printMOPs(out);
        out.println();
        printEventsRates(out);
        out.println();
    }

    /**
     * Returns a String with a description of the Model: the States
     * and the Transition Matrix. Its use is not recommended for large
     * models.
     * @return A string with the Matrix
     */
    public String denseMatrixToString() {
        StringWriter stw = new StringWriter();
        printDenseMatrix(new PrintWriter(stw));
        return stw.toString();
    }

    /**
     * Returns the Transition Matrix as a String. Its use is not
     * recommended for large models.
     * @param width The width of each column.
     * @param rateDecimals The number of decimals for the rates.
     * @param printZeros Whether zeros or blanks should be printed.
     * @param useGenerator whether the generator matrix <b>Q</>,
     *        rather than the rates matrix should be printed.
     * @return A String description of the rates or generator matrix.
     */

    public String denseMatrixToString(int width, int rateDecimals,
            boolean printZeros, boolean useGenerator) {
        StringWriter stw = new StringWriter();
        printDenseMatrix(new PrintWriter(stw), width, rateDecimals, printZeros,
                useGenerator);
        return stw.toString();
    }

    /**
     * Prints a the Transition Matrix. It will use default values for
     * width and decimals. Its use is not recommended for large
     * models.
     * @param out The writer to write to.
     */
    public void printDenseMatrix(PrintWriter out) {
        printDenseMatrix(out, 10, 4, false, false);
    }

    /**
     * Prints a description of the Model using the given PrintWriter:
     * the States and the Transition Matrix. Its use is not
     * recommended for large models.
     * @param out The writer to write to.
     * @param width The width of each column.
     * @param rateDecimals The number of decimals for the rates.
     * @param printZeros Whether zeros or blanks should be printed.
     * @param useGenerator whether the generator matrix <b>Q</>,
     *        rather than the rates matrix should be printed.
     */
    public void printDenseMatrix(PrintWriter out, int width, int rateDecimals,
            boolean printZeros, boolean useGenerator) {
        printDenseMatrix(out, width, rateDecimals, printZeros, useGenerator,
                null);
    }

    /**
     * Prints a description of the Model using the given PrintWriter:
     * the States and the Transition Matrix. Its use is not
     * recommended for large models.
     * @param width The width of each column.
     * @param rateDecimals The number of decimals for the rates.
     * @param printZeros Whether zeros or blanks should be printed.
     * @param useGenerator whether the generator matrix <b>Q</>,
     *        rather than the rates matrix should be printed.
     * @param idx the indices of separators.
     */
    protected void printDenseMatrix(PrintWriter out, int width,
            int rateDecimals, boolean printZeros, boolean useGenerator,
            int idx[]) {
        setStatus(Status.WRITING);
        double[][] theMat;
        if (useGenerator) {
            out.println("GENERATOR MATRIX: ");
            theMat = getGenerator();
        } else {
            out.println("RATES MATRIX: ");
            theMat = getRates();
        }
        States<S> stts = getStates();
        int n = getNumStates();
        if (idx == null)
            (idx = new int[1])[0] = n + 1;// No separators!
        int w = statesLableMaxWidth(width);
        out.print(pad(" ", w) + vLine());
        Iterator<S> itr = stts.iterator();
        // print states label.
        for (int i = 0, k1 = 0; i < n; i++) {
            out.print(pad(itr.next().toString(), w));
            if (i == idx[k1]) {
                out.print(vLine());
                k1++;
            }
        }
        // print the matrix
        String hLine = "\n" + hLine(w * (n + 1) + idx.length);
        out.print(hLine);
        itr = stts.iterator();
        for (int k1 = 0, i = 0; i < n; i++) {
            cnt = i;
            out.println();
            out.print(pad(itr.next().toString(), w) + vLine());
            for (int j = 0, k2 = 0; j < n; j++) {
                if (theMat[i][j] == 0 && !printZeros) {
                    out.print(pad("", w));
                } else {
                    out.print(pad(theMat[i][j], w, rateDecimals));
                }
                if (j == idx[k2]) {
                    out.print(vLine());
                    k2++;
                }
            }
            if (i == idx[k1]) {
                out.print(hLine);
                k1++;

            }
        }
        out.flush();
        setStatus(Status.GENERATED);
    }

    /**
     * Prints a description of the States and the Equilibrium
     * Probabilities.
     * @return A String description of the States.
     */

    public String statesToString() {
        StringWriter stw = new StringWriter();
        printStates(new PrintWriter(stw), 10, 5);
        return stw.toString();
    }

    /**
     * Prints a description of the States and the Equilibrium
     * Probabilities.
     * @param out The writer to write to.
     */

    public void printStates(PrintWriter out) {
        printStates(out, 10, 5);
    }

    /**
     * Prints a description of the States and the Equilibrium
     * Probabilities.
     * @param out The writer to write to.
     * @param width The width of each column.
     * @param probDecimals The number of decimals for the
     *        probabilities.
     */

    public void printStates(PrintWriter out, int width, int probDecimals) {
        this.out = out;
        States<S> stts = getStates();
        debug(1, "Printing States...");
        int maxL = 9; // Lets find max label width
        double[] pi;
        try {
            pi = getSteadyState();
            for (S s : stts) {
                maxL = Math.max(maxL, s.label().length());
            }
            int w = maxL + 3;
            int w2 = width + 3;
            // Print headers
            out.println(pad("", w) + pad("EQUILIBRUM", w2, false));
            out.println(pad("STATE", w, false) + pad("PROBAB.", w2, false)
                    + "DESCRIPTION");
            out.println();
            int i = 0;
            for (S s : stts) {
                out
                        .println(pad(s.label(), w, false)
                                + pad(pi[i++], w, probDecimals, false)
                                + ((s.description() != "") ? ""
                                        + s.description() : ""));
            }
        } catch (NotUnichainException e) {
            out.println(e);
        }
        out.println();
    }

    /**
     * This method should be implemented by the subclass to give word
     * description of the model. For example it should say: Queueing
     * system with 2 servers and exponential arrivals with rate 50 per
     * hour.
     * @return A description of the Model
     */
	public abstract String description();

    /**
     * pad fills with blanks up to width w
     * @param s The String to print.
     * @param w The width to pad up to.
     */
    protected String pad(String s, int w) {
        return pad(s, w, true);
    }

    /**
     * pad fills with blanks up to width w
     * @param s The String to print.
     * @param w The width to pad up to.
     * @param right whether the string should be eligned to the right.
     */
    protected String pad(String s, int w, boolean right) {
        String stg = "";
        int padding = Math.max(1, w - s.length()); // At _least_ 1
        // space
        for (int k = 0; k < padding; k++)
            stg += ' ';
        return right ? (stg + s) : (s + stg);
    }

    /**
     * pad generates a string representing the double v, padded with
     * spaces up to width w. *
     * @param v The double to print.
     * @param w The width to pad up to.
     */

    protected String pad(double v, int w) {
        return pad(v, w, true);
    }

    /**
     * pad generates a string representing the double v, padded with
     * spaces up to width w.
     * @param v The double to print to print.
     * @param w The width to pad up to.
     * @param right whether the string should be eligned to the right.
     */
    protected String pad(double v, int w, boolean right) {
        int d = 3;
        return pad(v, w, d, right);
    }

    /**
     * pad fills with blanks up to width w. Alignment is to the right.
     * @param v The number to print.
     * @param w The width to pad up to.
     * @param d number of decimals.
     */
    protected String pad(double v, int w, int d) {
        return pad(v, w, d, true);
    }

    /**
     * pad fills with blanks up to width w
     * @param v The number to print.
     * @param w The width to pad up to.
     * @param d number of decimals.
     * @param right whether the string should be eligned to the right.
     */
    protected String pad(double v, int w, int d, boolean right) {
        DecimalFormat format = new DecimalFormat();
        format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
        format.setMinimumIntegerDigits(1);
        format.setMaximumFractionDigits(d);
        format.setMinimumFractionDigits(d);
        format.setGroupingUsed(false);
        String s = format.format(v); // format the number
        return pad(s, w, right);
    }

    /**
     * Computes the maximum used by the state's labels.
     * @param width minimum width acceptable.
     * @return Max(width, max label).
     */
    public int statesLableMaxWidth(int width) {
        States<S> stts = getStates();
        int maxL = 0; // Lets find max label width
        for (S s : stts) {
            maxL = Math.max(maxL, s.label().length());
        }
        return Math.max(maxL + 1, width);
    }

    /**
     * Retuns an horizontal text line of the given length.
     * @param length
     * @return a horizontal line.
     */
    protected String hLine(int length) {
        StringBuilder bl = new StringBuilder();
        for (int j = 0; j < length; j++)
            bl.append("\u2014");// u2015 '?'
        return bl.toString();
    }

    /**
     * Returns a text vertical line.
     * @return vertical line.
     */
    protected String vLine() {
        return "|";// u200C "?"
    }

    /*
     * The next statements are used for debug control.
     */
    /**
     * Sets the DebugReporter to use.
     * @param reporter The reporter tah will capture the debug
     *        information.
     * @see DebugReporter
     */
    public void setDebugReporter(DebugReporter reporter) {
        this.reporter = reporter;
    }

    /**
     * Gets the DebugReporter currently in use.
     * @return a DebugReporter where debug information is sent.
     * @see DebugReporter
     */
    public DebugReporter getDebugReporter() {
        return reporter;
    }

    /**
     * Prints debug information with this importance level
     * @param level
     * @param s The message to send.
     */
    public void debug(int level, String s) {
        reporter.debug(level, s, true);
    }

    /**
     * Prints debug information with this importance level
     * @param level The level of importance (0=show always, 5= show on
     *        debug level is 5).
     * @param s The message
     * @param newline Whether a new line should be written.
     */
    public void debug(int level, String s, boolean newline) {
        reporter.debug(level, s, newline);
    }

    /**
     * Prints debug information with this importance level
     * @param level The level of importance (0=show always, 5= show on
     *        debug level is 5).
     * @param s The string to write.
     * @param newline Whether to use a new line.
     * @param indent Whether it should indent according to level.
     */
    public void debug(int level, String s, boolean newline, boolean indent) {
        reporter.debug(level, s, newline, indent);
    }

    /**
     * @return current debug level, where level=0 means no debug info
     *         and level = 4 verbose info.
     */
    public int getDebugLevel() {
        return reporter.getCurLevel();
    }

    /**
     * Sets the debug level, where level=0 means no debug info, level =
     * 5 verbose info.
     * @param level New debug level
     */
    public void setDebugLevel(int level) {
        reporter.setCurLevel(level);
    }

    // / UI stuff

    private MarkovGUI theGUI = null;

    /**
     * Loads the Graphic User Interface (GUI) that represent this
     * Markov Chain.
     */
    public void loadGUI() {
        if (theGUI == null) {
            debug(1, "Loading Graphic User Interface...");
            theGUI = new MarkovGUI(this);
        }
    }

    /**
     * Shows the Graphic User Interface (GUI) that represent this
     * Markov Chain. It loads one if none has been defined.
     */
    public void showGUI() {
        if (theGUI == null) {
            loadGUI();
        }
        if (theGUI != null) {
            theGUI.setVisible(true);
        }
    }

    /**
     * Hides the Graphic User Interface (GUI) that represent this
     * Markov Chain if one is defined.
     */
    public void hideGUI() {
        if (theGUI != null) {
            theGUI.setVisible(false);
        }
    }

    /**
     * Destroys the Graphic User Interface (GUI) that represent this
     * Markov Chain if one is defined.
     */
    public void killGUI() {
        if (theGUI != null) {
            theGUI.setVisible(false);
            theGUI = null;
        }
    }

    /**
     * Allos to stop model execution by graphica user interface.
     * @return Used to check with the GUI if the user has requested to
     *         stop.
     */
    public boolean canGo() {
        return (theGUI == null) ? true : theGUI.isAllowedToRun();
    }

    /**
     * Pauses the current execution of the model. This is called by
     * the GUI.
     */
    public void pause() {
        // canGo = false;
        setStatus(Status.SUSPENDED);
    }

    /**
     * Runs the model, or resumes execution if it had been suspended.
     * Its use is intended for Graphical User Interface(GUI). A
     * standard user should use <code>generate()</code> instead.
     * @see SimpleMarkovProcess#generate()
     */

    public synchronized void go() {
        // canGo = true;
        setStatus(Status.RUNNING);
        generate();
    }

    /**
     * Runs the model for a single step. Its use is intended for
     * Graphical User Interface(GUI). A standard user should use
     * <code>generate()</code> instead.
     * @see SimpleMarkovProcess#generate()
     * @see SimpleMarkovProcess#go()
     */

    public synchronized void goStep() {
        // canGo = true;
        if (!uncheckedStates.isEmpty() && cnt < maxStates) {
            setStatus(Status.RUNNING);
            generateStep();
        }
        setEndStatus();
    }

    /**
     * Returns the current status of the model.
     * @return One of the constants IDLE, RUNNING, GENERATED,
     *         SUSPENDED, ERROR
     */

    public Status getStatus() {
        return _status;
    }

    /**
     * Return the number of states processed so far in the current
     * process.
     * @return the number of states processed so far.
     */
    public long getProgress() {
        return cnt;
    }

    /**
     * Returns a String describing the current status of the model.
     * @return a String describing the current status of the model.
     */
    public String getStatusMsg() {
        String stg = "";
        switch (_status) {
        case IDLE:
            stg = "Model has not been generated. ";
            break;
        case RUNNING:
            stg = "Model is running. ";
            break;
        case SUSPENDED:
            stg = "Model is suspended. ";
            break;
        case ERROR:
            stg = "Model generation caused an ERROR. ";
            break;
        case GENERATED:
            stg = "Model was succesfully generated. It has " + getNumStates()
                    + " states.";
            break;
        case WRITING:
            stg += "Writing information";
            break;
        default: // do nothing
        }
        return stg;
    }

    /**
     * Sets the value for the current status of the model.
     * @param status
     */
    private void setStatus(Status status) {
        _status = status;
        if (theGUI != null)
            theGUI.updateStatus();
    }

    /**
     * @return the Maximum number of states to generate.
     */
    public long getMaxStates() {
        return maxStates;
    }

    /**
     * Sets the maximum number of states to generate. Increase this
     * only if you are sure that your model has a big number od
     * states. Your current License may prevent you from setting this
     * number.
     * @param num Maximum Number of States to generate.
     */
    public void setMaxStates(long num) {
        maxStates = num;
    }

    /*
     * SOLVERS MANAGEMENT
     */

    /**
     * Returns the default SteadyStateSolver.
     * @return the default SteadyStateSolver.
     */
    protected final SteadyStateSolver getDefaultSteadyStateSolver() {
        if (defaultSteadyStateSolver == null) {
            defaultSteadyStateSolver = new MtjSolver(this);
            // defaultSteadyStateSolver = new JamaSolver(this);
        }
        return defaultSteadyStateSolver;
    }

    /**
     * The currently defined solver.
     * @see jmarkov.solvers.SteadyStateSolver
     * @return Returns the steadyStateSolver.
     */
    public SteadyStateSolver getSteadyStateSolver() {
        if (steadyStateSolver == null)
            steadyStateSolver = getDefaultSteadyStateSolver();
        return steadyStateSolver;
    }

    /**
     * Allows the user to set an alternate solver.
     * @see jmarkov.solvers.SteadyStateSolver
     * @param steadyStateSolver The steadyStateSolver to set.
     */
    public void setSteadyStateSolver(SteadyStateSolver steadyStateSolver) {
        this.steadyStateSolver = steadyStateSolver;
        resetResults();
        debug(1, "New steady state solver set: " + steadyStateSolver);
    }

    /**
     * The default solver for transient state.
     * @see jmarkov.solvers.TransientSolver
     * @return Returns the default transientSolver.
     */
    protected TransientSolver getDefaultTransientSolver() {
        if (defaultTransientSolver == null) {
            defaultTransientSolver = new JamaTransientSolver(this);
        }
        return transientSolver;
    }

    /**
     * The currently defined solver for transient state.
     * @see jmarkov.solvers.TransientSolver
     * @return Returns the transientSolver.
     */
    public TransientSolver getTransientSolver() {
        if (transientSolver == null) {
            transientSolver = getDefaultTransientSolver();
        }
        return transientSolver;
    }

    /**
     * Allows the user to set an alternate solver.
     * @see jmarkov.solvers.TransientSolver
     * @param transientSolver The transientSolver to set.
     */
    public void setTransientSolver(TransientSolver transientSolver) {
        this.transientSolver = transientSolver;
        resetResults();
        debug(1, "New transient solver set: " + transientSolver);
    }

    /**
     * Returns the name of the model.
     * @return The Name of this model.
     */
	public String label() {
        if (name == "") {
            name = this.getClass().getSimpleName();
        }
        return name;
    }

    // / FINALIZE

    /*
     * @see java.lang.Object#finalize()
     */
    @Override
    protected void finalize() throws Throwable {
        super.finalize();
        debug(1, "Model " + label() + " unloaded from memory.");
    }

} // end class
