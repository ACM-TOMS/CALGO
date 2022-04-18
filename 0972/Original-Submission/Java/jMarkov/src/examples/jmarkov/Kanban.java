package examples.jmarkov;

import jmarkov.MarkovProcess;
import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;

/**
 * This class implements a Kanban system with one card. Extends
 * SimpleMarkovProcess.
 * @author Margarita Arana, Gloria Diaz, Germán Riaño. 
 * @version 1.0a
 */
public class Kanban extends SimpleMarkovProcess<KanbanState, KanbanEvent> {
    private double[] mu;
    private int m;
    private int[] kanban;

    /**
     * Kanban constructor.
     * @param m Number of stations
     * @param mu Service rates at the stations 
     * @param kanban Number of cards in the stations
     * @param s Status of the stations
     */
    public Kanban(int m, int[] kanban, double[] mu, int[] s) {
        super((new KanbanState(s, new int[m], kanban)), KanbanEvent
                .getAllEvents(m));
        this.m = m;
        this.mu = mu;
        this.kanban = kanban;
    }

    /**
     * Default Constructor for GUI
     */
    public Kanban() {
        this(3, new int[] { 2, 2 }, new double[] { 2.0, 2.0, 2.0 }, new int[] {
                1, 0, 0 });
    }

    /**
     * Describes the system 
     * @see jmarkov.SimpleMarkovProcess#description()
     */
    @Override
    public String description() {
        String stg = "KANBAN SYSTEM \n" + m + " Stations.";
        for (int i = 0; i < m; i++) {
            stg += "\nStation "
                    + (i + 1)
                    + ": Rate = "
                    + mu[i]
                    + ((i < m - 1) ? ", Kanbans after station: " + kanban[i]
                            : "");
        }

        return stg;
    }

    /**
     * Determine the active events
     */
    @Override
    public boolean active(KanbanState i, KanbanEvent e) {
        return (i.getStatus(e.getStation()) == 1);
    }

    @Override
    public States<KanbanState> dests(KanbanState i, KanbanEvent ev) {
        int[] newBuff = new int[m];
        int e = ev.getStation();
        for (int j = 1; j < m; j++) {
            newBuff[j] = i.getBuffer(j);
        }

        int[] newKanban = new int[m - 1];
        for (int j = 0; j < m - 1; j++) {
            newKanban[j] = i.getKanban(j);
        }

        int[] newStatus = new int[m];
        for (int j = 0; j < m; j++) {
            newStatus[j] = i.getStatus(j);
        }

        newStatus[e] = 0;
        if (e== 0) { 
            newBuff[1] = newBuff[1] + 1; 
                                          
            if (m == 2 && newStatus[1] == 0) {
                newStatus[1] = 1;
                newKanban[0] = newKanban[0] + 1;
                newBuff[1] = newBuff[1] - 1;
            }
            if (m > 2 && newKanban[1] > 0 && newStatus[1] == 0) {

                newStatus[1] = 1;
                newKanban[0] = newKanban[0] + 1;
                newKanban[1] = newKanban[1] - 1;
                newBuff[1] = newBuff[1] - 1;
            }
            if (newKanban[0] > 0) { 
                newStatus[0] = 1;
                newKanban[0] = newKanban[0] - 1;
            }
            return new StatesSet<KanbanState>( new KanbanState(newStatus, newBuff, newKanban));
        } // end if(e==0)

        if (e == m - 1) { 
            if (newBuff[e] > 0) {
                newStatus[e] = 1;
                newKanban[e - 1] = newKanban[e - 1] + 1;
                newBuff[e] = newBuff[e] - 1;
            }
        } // end if(e==m-1)
        else { 
            newBuff[e + 1] = newBuff[e + 1] + 1;
            
            if ((e < m - 2) && (newKanban[e + 1] > 0)
                    && (newStatus[e + 1] == 0)) {
               
                newStatus[e + 1] = 1; 
                newKanban[e] = newKanban[e] + 1; 
                newKanban[e + 1] = newKanban[e + 1] - 1;
                
                newBuff[e + 1] = newBuff[e + 1] - 1;
                
            }
            if (e == m - 2 && newStatus[e + 1] == 0) {
                
                newStatus[e + 1] = 1;
                
                newKanban[e] = newKanban[e] + 1; 
                newBuff[e + 1] = newBuff[e + 1] - 1; 
            }
            if (newKanban[e] > 0 && newBuff[e] > 0) {
                
                newStatus[e] = 1; 
                newKanban[e - 1] = newKanban[e - 1] + 1;
                
                newKanban[e] = newKanban[e] - 1;
                
                newBuff[e] = newBuff[e] - 1; 
            }
        } // end else

        for (int j = e - 1; j >= 0; j--) {
            if (newStatus[j] == 1) { 
                return new StatesSet<KanbanState>(new KanbanState(newStatus, newBuff, newKanban));
            }

            if (j == 0) {
                if (newKanban[j] > 0) {
                    newStatus[j] = 1;
                    newKanban[j] = newKanban[j] - 1;
                }
                return new StatesSet<KanbanState>(new KanbanState(newStatus, newBuff, newKanban));
            } // end if(j==0)

            if (newBuff[j] == 0 || newKanban[j] == 0) { 
                return new StatesSet<KanbanState>(new KanbanState(newStatus, newBuff, newKanban));
            }
            newStatus[j] = 1;
            newKanban[j - 1] = newKanban[j - 1] + 1;
            newKanban[j] = newKanban[j] - 1;
            newBuff[j] = newBuff[j] - 1;
        } // end for

        return new StatesSet<KanbanState>(new KanbanState(newStatus, newBuff, newKanban));
    } // end de la funciï¿½n dest

    /**
     * Get the service rate of each station
     * @see SimpleMarkovProcess#rate(State, State, Event)
     */
    @Override
    public  double rate(KanbanState i, KanbanState j, KanbanEvent e) {
        return mu[e.getStation()];
    }

    /**
     * Main method
     * @param a Not used
     */
    public static void main(String[] a) {

        int[] tarj = { 2, 2 };
        double[] mu = { 2.0, 2.0, 2.0 };
        int[] s = { 1, 0, 0 };
        Kanban theKanban = new Kanban(s.length, tarj, mu, s);
        theKanban.showGUI();
        theKanban.generate();
        theKanban.printAll();
    }

} // class end

/**
 * This class define the states of a Kanban. Extends PropertiesState.
 */
class KanbanState extends PropertiesState {
    /**
     * Builds a new Kanban state.
     * @param s marks which stations are active.
     * @param x array with the number of pieces in each buffer
     * @param y array with the number of cards in each box
     */
    KanbanState(int[] status, int[] buffer, int[] kanbans) {
        super(status.length + (buffer.length - 1) + kanbans.length);
        int m = status.length;
        for (int i = 0; i < m; i++) {
            prop[i] = status[i];
        }
        for (int i = 0; i < m - 1; i++) {
            prop[i + m] = buffer[i + 1];
        }
        for (int i = 1; i < m; i++) {
            prop[i + (2 * (m - 1))] = kanbans[i - 1];
        }
    }

    @Override
    public void computeMOPs(MarkovProcess mp) {
        int m = getNumStations();
        for (int i = 0; i < m; i++) {
            setMOP(mp, "Server utilization " + (i + 1), getStatus(i));
        }
        for (int i = 0; i < m - 1; i++) {
            setMOP(mp, "WIP in Kanbans " + (i + 1), getKanban(i));
        }
        for (int i = 1; i < m; i++) {
            setMOP(mp, "WIP in Buffer " + i, getBuffer(i));
        }
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
     * Get the status of station i
     * @param i index of the station
     * @return 1 if the station is active
     */
    public int getStatus(int i) { 
    	// The stations indexes start from zero
        return prop[i];
    }

    /**
     * Get the quantity of material in buffer i
     * @param i buffer index
     * @return Quantity of material in buffer i.
     */
    public int getBuffer(int i) { 
    	// The buffers indexes start from one
        int m = (prop.length + 2) / 3;
        return prop[i + (m - 1)];
    }

    /**
     * Get the number of cards in box i
     * @param i box index
     * @return Number of cards in box i
     */
    public int getKanban(int i) { 
    	// The boxes indexes start from zero
        int m = getNumStations();
        return prop[(i + 1) + 2 * (m - 1)];
    }

    /**
     * Get the number of stations
     * @return Number of stations
     */
    public int getNumStations() {
        return (prop.length + 2) / 3;
    }

    @Override
    public String label() {
        String stg = "";
        int m = (prop.length + 1) / 3 + 1;
        for (int i = 0; i < m; i++) {
            stg += "S" + (i + 1) + "(" + getStatus(i) + ","
                    + ((i > 0) ? "" + getBuffer(i) : "*") + ","
                    + ((i < m - 1) ? "" + getKanban(i) : "*") + ")";
        }
        return stg;
    }

    @Override
    public String description() {
        String stg = "";
        int m = (prop.length + 1) / 3 + 1;
        for (int i = 0; i < m; i++) {
            stg += "STATION " + (i + 1) + "= (Status: " + getStatus(i)
                    + ", Buffer: " + ((i > 0) ? "" + getBuffer(i) : "*")
                    + ", Kanbans: " + ((i < m - 1) ? "" + getKanban(i) : "*")
                    + ")  ";
        }
        return stg;
    }

}

class KanbanEvent extends Event {
    private int station;

    /**
     * Creates an ending process at the given station.
     * @param station
     */
    public KanbanEvent(int station) {
        this.station = station;
    }

    /**
     * @return Returns the station.
     */
    public int getStation() {
        return station;
    }

    /**
     * @param m
     * @return a set with all Events
     */
    public static EventsSet<KanbanEvent> getAllEvents(int m) {
        KanbanEvent set[] = new KanbanEvent[m];
        for (int s = 0; s < m; s++)
            set[s] = new KanbanEvent(s);
        return new EventsSet<KanbanEvent>(set);
    }
}
