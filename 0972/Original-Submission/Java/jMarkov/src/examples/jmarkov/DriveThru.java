package examples.jmarkov;

import static examples.jmarkov.DriveThruEvent.Type.ARRIVAL;
import static examples.jmarkov.DriveThruEvent.Type.MIC_COMPLETION;
import static examples.jmarkov.DriveThruEvent.Type.SERVICE_COMPLETION;
import static examples.jmarkov.DriveThruState.CustStatus.BLOCKED_DONE;
import static examples.jmarkov.DriveThruState.CustStatus.COOKING;
import static examples.jmarkov.DriveThruState.CustStatus.EMPTY;
import static examples.jmarkov.DriveThruState.CustStatus.ORDERING;
import static examples.jmarkov.DriveThruState.CustStatus.WAIT_MIC;

import java.io.PrintWriter;

import jmarkov.MarkovProcess;
import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.exceptions.NotUnichainException;
import examples.jmarkov.DriveThruState.CustStatus;

/**
 * This class implements a Drive Thru. Extends 
 * SimpleMarkovProcess.
 * 
 * @author Margarita Arana y Gloria Díaz.  Universidad de los Andes.
 * Mod: Germán Riaño (2004)
 * @version 1.0a
 */
public class DriveThru extends
        SimpleMarkovProcess<DriveThruState, DriveThruEvent> {

    double lambda; // arrival rate
    double mu1; // Service rate for server 1
    double mu2; // Service rate for server 2
    int M; // Maximum number of clients in the system
    int S; // Number of servers
    int N; // Number of places between the window and the microphone

    /**
     * Constructor de un DriveThru.
     * 
     * @param lambda
     *            Tasa de arribos
     * @param mu1
     *            Tasa de servicios del micrï¿½fono
     * @param mu2
     *            Tasa de servicios de la ventana
     * @param M
     *            Nï¿½mero mï¿½ximo de entidades en el sistema
     * @param S
     *            Nï¿½mero de servidores
     * @param N
     *            Nï¿½mero de puestos entre la ventana y el micrï¿½fono
     */
    public DriveThru(double lambda, double mu1, double mu2, int M, int S, int N) {
        super((new DriveThruState(N, S)), DriveThruEvent.getAllEvents(N));
        this.lambda = lambda;
        this.mu1 = mu1;
        this.mu2 = mu2;
        this.M = M;
        this.S = S;
        this.N = N;

    }

    /**
     * Default constructor for GUI.
     */
    public DriveThru() {
        this(80.0, 12.0, 30.0, 4, 2, 1);
    }

    /**
     * Determines when the states are active for each state.
     * 
     * @see SimpleMarkovProcess#active(State, Event)
     */

    @Override
    public boolean active(DriveThruState s, DriveThruEvent ev) {
        boolean result = false;
        switch (ev.getType()) {
        case ARRIVAL:
            // un carro puede llegar si hay espacio en cola
            result = (s.getQLength() < M - N - 1);
            break;
        case MIC_COMPLETION:
            // se puede terminar de tomar la orden si una persona esta haciendo
            // el pedido
            result = (s.getMicStatus() == ORDERING);
            break;
        default:
            // se puede terminar una orden si la persona correspondiente la esta
            // esperando
            if (ev.getPos() == N) {
                result = (s.getMicStatus() == COOKING);
            } else {
                result = (s.getStatus(ev.getPos()) == COOKING);
            }
        }
        return result;
    }

    /**
     * Computes the rate: the rate is lambda if an arraival occurs, 
     * the rate is mu1 if a service type one is finished,
     * the rate is mu2 if an service type two is finished.
     * 
     * @see SimpleMarkovProcess#rate(State, State, Event)
     */
    @Override
    public double rate(DriveThruState i, DriveThruState j, DriveThruEvent e) {
        switch (e.getType()) {
        case ARRIVAL:
            return lambda;
        case MIC_COMPLETION:
            return mu1;
        default:
            return mu2;
        }
    }

    /**
     * Computes the status of the destination when an event occurs
     * 
     * @see SimpleMarkovProcess#dests(State, Event)
     */

    @Override
    public States<DriveThruState> dests(DriveThruState i, DriveThruEvent e) {
        int numServ = i.getAvlServs();
        CustStatus[] status = i.getStatus();
        CustStatus newMic = i.getMicStatus();
        int newQsize = i.getQLength();
        int numGone = 0;
        boolean micMoves = false;
        int k; // utility counter

        switch (e.getType()) {
        case ARRIVAL:
            if (i.getMicStatus() == EMPTY && numServ > 0) {
                
                newMic = ORDERING;
                numServ = numServ - 1;
            } else if (i.getMicStatus() == EMPTY && numServ == 0) {
                
                newMic = WAIT_MIC;
            } else if (i.getQLength() < M - N - 1) {
                
                newQsize = i.getQLength() + 1;
            }
            break;

        case MIC_COMPLETION:
            newMic = COOKING;
            for (k = 0; ((k < N) && (status[k] != EMPTY)); k++)
                ;
           
            if (k != N) { 
                status[k] = COOKING;
                newMic = EMPTY;
                micMoves = true;
            }
            break;

        default: 
            numServ = numServ + 1; 
            int p = e.getPos();
            if (p > 0 && p < N) {
                
                status[p] = BLOCKED_DONE;
            } else if (p == N) { 
                newMic = BLOCKED_DONE;
            } else { 
                
                status[0] = EMPTY;
                
                int pos1, pos2;
                
                for (k = 1; ((k < N) && status[k] == BLOCKED_DONE); k++)
                    ;
                numGone = k; 
                if (k != N) { 
                    pos1 = k;
                    pos2 = N - 1;
                    for (k = pos1; k <= pos2; k++) {
                        status[k - numGone] = status[k];
                    }
                }
                for (k = N - numGone; k < N; k++) {
                    status[k] = EMPTY;
                }
                if (newMic == COOKING) {
                    status[N - numGone] = newMic;
                    newMic = EMPTY;
                    micMoves = true;
                } else if (newMic == BLOCKED_DONE) {
                    newMic = EMPTY;
                    micMoves = true;
                }
            }
            break;
        } // end switch

        if (newMic == WAIT_MIC && numServ > 0) {
            newMic = ORDERING;
            numServ--;
        }
        if (micMoves) {
            if (i.getQLength() > 0 && numServ > 0) {
                newMic = ORDERING;
                numServ = numServ - 1;
                newQsize = i.getQLength() - 1;
            } else if (i.getQLength() > 0 && numServ == 0) {
                newMic = WAIT_MIC;
                newQsize = i.getQLength() - 1;
            }
        }
        StatesSet<DriveThruState> set = new StatesSet<DriveThruState>();
        set.add(new DriveThruState(status, newMic, newQsize, numServ));
        return set;
    } // end dests

    @Override
    public String description() {
        return "SISTEMA DRIVE THRU. " + "\nTasa de Entrada   = " + lambda
                + "\nTasa en el Mic    = " + mu1 + "\nTasa de sevicio 2 = "
                + mu2 + "\nPosiciï¿½n del mic  = " + N + "\nServidores        = "
                + S + "\nCap en el sistema = " + M;
    }

    /**
     * Print all waiting times associated with each MOP
     */
    @Override
    public int printMOPs(PrintWriter out, int width, int decimals) {
        int namesWidth = super.printMOPs(out, width, decimals);
        // this rate work for all MOPs
        double ldaEff;
        try {
            ldaEff = getEventRate(ARRIVAL.ordinal());
            String[] names = getMOPNames();
            double waitTime;
            int N = names.length;
            namesWidth += 20;
            for (int i = 0; i < N; i++) {
                waitTime = 60 * getMOPsAvg(names[i]) / ldaEff;
                String name = "Waiting time for " + names[i];
                out.println(pad(name, namesWidth, false)
                        + pad(waitTime, width, decimals) + " minutes");
            }
        } catch (NotUnichainException e) {
            out.println(e);
        }
        return namesWidth;
    }

    /**
     * Main method.
     * 
     * @param a
     *            Not used.
     */
    public static void main(String[] a) {
        // as in handout:
        DriveThru theDT = new DriveThru(80.0, 12.0, 30.0, 4, 2, 1);
        // DriveThru theDT = new DriveThru(80.0, 120.0, 30.0, 4, 2, 2);
        theDT.setDebugLevel(5);

        theDT.showGUI();
        theDT.printAll();
        theDT.printMOPs();
    }

} // class end

/**
 * This is a particular case of PropertiesState. Here, N is the position of the 
 * microphone. The first N-1 components represent the status of the first queue, the 
 * component N is the status of the microphone, the component N+1 is the number of clients in 
 * the queue, and N+2 are the available servers.
 */
class DriveThruState extends State {

    // private int micPos;
    // private CustStatus micStatus;
    private int numQ;
    private int avlServ;
    private CustStatus[] prop = null;

    /**
     * This enumeration shows the different status for a customer.
     * 
     */
    public enum CustStatus {
        /** Empty space. */
        EMPTY,
        /** In service. */
        ORDERING,
        /** A client in the microphone, but there are no servers available. */
        WAIT_MIC,
        /** The client order is being prepared. */
        COOKING,
        /** The order is ready but the client is blocked. */
        BLOCKED_DONE;
    }

    /**
     * Builds a State representing an empty system
     * 
     * @param micPos
     * @param serv
     */
    DriveThruState(int micPos, int serv) {
        this(new CustStatus[micPos], EMPTY, 0, serv);
        for (int i = 0; i < prop.length; i++) {
            prop[i] = EMPTY;
        }
    }

    /**
     * Builds a DriveThru state.
     * 
     * @param vec
     *            The states from the window until the microphone, 
     *            without including the microphone.
     * @param mic
     *            Microphone status.
     * @param numQ
     *            Number of clients in the queue.
     * @param avServs
     *            Number of servers available.
     */

    DriveThruState(CustStatus[] statusVec, CustStatus micStatus, int numQ,
            int avServs) {
        prop = new CustStatus[statusVec.length + 1];
        int micPos = statusVec.length;
        System.arraycopy(statusVec, 0, prop, 0, micPos);
        prop[micPos] = micStatus;
        this.numQ = numQ;
        this.avlServ = avServs;
    }

    /**
     * Compute all the MOPs for this state
     */
    @Override
    public void computeMOPs(MarkovProcess mp) {
        int servEtapa1 = 0;
        int servEtapa2 = 0;
        int blockedDone = 0;
        int blockedBefore = 0;
        int total = 0;
        for (CustStatus s : prop) {
            servEtapa1 += (s == ORDERING) ? 1 : 0;
            servEtapa2 += (s == COOKING) ? 1 : 0;
            blockedDone += (s == BLOCKED_DONE) ? 1 : 0;
            blockedBefore += (s == WAIT_MIC) ? 1 : 0;
            total += (s != EMPTY) ? 1 : 0;
        }
        setMOP(mp, "Tamano Cola", getQLength());
        setMOP(mp, "Serv Ocupados Microfono ", servEtapa1);
        setMOP(mp, "Serv Ocupados Cocinando", servEtapa2);
        setMOP(mp, "Serv Ocupados ", servEtapa1 + servEtapa2);
        setMOP(mp, "Clientes Bloqueados antes de ordenar", blockedBefore);
        setMOP(mp, "Clientes Bloqueados con orden lista", blockedDone);
        setMOP(mp, "Clientes Bloqueados", blockedBefore + blockedDone);
        setMOP(mp, "Total clientes en Espera", blockedBefore + blockedDone
                + getQLength());
        setMOP(mp, "Total Clientes ", total + getQLength());
    }

    /**
     * Get the number of clients in the queue.
     * 
     * @return Number of clients in the queue.
     */
    public int getQLength() {
        return numQ;
    }

    /**
     * Get the status of the of the i-th component.
     * 
     * @param i
     *         index of the component
     * 
     * @return Status of the i-th component.
     */
    public CustStatus getStatus(int i) {
        return prop[i];
    }

    /**
     * Get the vector of clients statuses. 
     * 
     * @return Status of components 0 to N-1.
     */
    public CustStatus[] getStatus() {
        int micPos = getMicPos();
        CustStatus[] status = new CustStatus[micPos];
        System.arraycopy(prop, 0, status, 0, micPos);
        return status;
    }

    /**
     * Get the status of the window.
     * 
     * @return The status of the client at the microphone.
     */
    public CustStatus getMicStatus() {
        int n = prop.length - 1;
        return prop[n];
    }

    /**
     * Return the mic position.
     * 
     * @return mic position index
     */
    public int getMicPos() {
        return prop.length - 1;
    }

    /**
     * Get the status of the window
     * 
     * @return Status of the window.
     */
    public CustStatus getVentana() {
        return prop[0];
    }

    /**
     * Computes the number of available servers.
     * 
     * @return Number of available servers.
     */
    public int getAvlServs() {
        return avlServ;
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
        String stg = "";
        for (CustStatus s : prop) {
            switch (s) {
            case EMPTY:
                stg += "0";
                break;
            case ORDERING:
                stg += "m";
                break;
            case WAIT_MIC:
                stg += "w";
                break;
            case COOKING:
                stg += "c";
                break;
            case BLOCKED_DONE:
                stg += "b";
                break;
            }
        }
        return stg + "Q" + numQ;
        // return stg + "Q" + prop[micPos + 1] + "S" + prop[micPos + 2];
    }

    String statusDesc(CustStatus stat) {
        switch (stat) {
        case EMPTY:
            return "empty";
        case ORDERING:
            return "ordering,";
        case WAIT_MIC:
            return "waiting";
        case COOKING:
            return "cooking";
        default: // DONE
            return "blocked";
        }
    }

    /**
     * Describes the State
     * 
     * @see jmarkov.basic.State#description()
     */
    @Override
    public String description() {
        String stg = "";
        int N = getMicPos();
        stg = "Queue CustStatus: (";
        for (int i = 0; i < N; i++) {
            stg += statusDesc(getStatus(i));
            stg += (i < N - 1) ? ", " : "";
        }
        stg += "). Mic status: " + statusDesc(getMicStatus());
        stg += ". Queue Size: " + getQLength();
        return stg;
    }

    /**
     * @see jmarkov.basic.State#compareTo(jmarkov.basic.State)
     */
    @Override
    public int compareTo(State j) {
        if (!(j instanceof DriveThruState))
            throw new IllegalArgumentException("Comparing wrong types!");
        DriveThruState u = (DriveThruState) j;
        int micPos = getMicPos();
        for (int k = 0; k <= micPos; k++) {
            if (getStatus(k).ordinal() > u.getStatus(k).ordinal())
                return +1;
            if (getStatus(k).ordinal() < u.getStatus(k).ordinal())
                return -1;
        }
        if (getQLength() > u.getQLength())
            return +1;
        if (getQLength() < u.getQLength())
            return -1;
        if (getAvlServs() > u.getAvlServs())
            return +1;
        if (getAvlServs() < u.getAvlServs())
            return -1;
        return 0;
    }

}

/**
 * This class implements the events in a Drive Thru.
 */
class DriveThruEvent extends jmarkov.basic.Event {
    /** Event types. */
    public static enum Type {
        /** Arrivale to the system. */
        ARRIVAL,
        /** Car at mic finishes service. */
        MIC_COMPLETION,
        /** Service completion for somebody who ordered. */
        SERVICE_COMPLETION;
    }

    private Type type; // event type
    private int position; // Position of the client whose order is complete

    /**
     * Creates an ARRIVAL or MIC_COMPLETION event.
     * 
     * @param type
     */
    public DriveThruEvent(Type type) {
        assert (type == ARRIVAL || type == MIC_COMPLETION);
        this.type = type;
    }

    /**
     * Creates a Service Completion event at he given position.
     * 
     * @param position
     *            Postion where the event occurs ( 0-based ).
     */
    public DriveThruEvent(int position) {
        this.type = SERVICE_COMPLETION;
        this.position = position;
    }

    /**
     * @return position where this event occurs. (valid only if type ==
     *         SERVICE_COMPLETION).
     */
    public int getPos() {
        assert (type == SERVICE_COMPLETION);
        return position;
    }

    /**
     * @return event type
     */
    public Type getType() {
        return type;
    }

    /**
     * @param micPos
     * @return A set with all the events in the system.
     */
    public static EventsSet<DriveThruEvent> getAllEvents(int micPos) {
        EventsSet<DriveThruEvent> eSet = new EventsSet<DriveThruEvent>();
        eSet.add(new DriveThruEvent(ARRIVAL));
        eSet.add(new DriveThruEvent(MIC_COMPLETION));
        for (int i = 0; i <= micPos; i++)
            eSet.add(new DriveThruEvent(i));
        return eSet;
    }

    @Override
    public String label() {
        String stg = "";
        switch (type) {
        case ARRIVAL:
            stg = "Arrival";
            break;
        case MIC_COMPLETION:
            stg = "MicEnd";
            break;
        default:
            stg = "SrvEnd(" + position + ")";
        }
        return stg;
    }

}
