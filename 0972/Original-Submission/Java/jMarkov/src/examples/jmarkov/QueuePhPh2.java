package examples.jmarkov;

import java.util.Date;

import jmarkov.GeomProcess;
import jmarkov.GeomRelState;
import jmarkov.MarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.StatesSet;

/**
 * This class implements a system whit two tandem servers each one with a phase
 * service. The first server has an infinite capacity and the second server has
 * finite capacity.
 * 
 * @author Julio Goez. Universidad de los Andes.
 */
public class QueuePhPh2 extends GeomProcess<QueuePhPh2State, QueuePhPh2Event> {

    final int arrival = QueuePhPh2Event.arrival;
    final int st1Ph1ToSt1Ph2 = QueuePhPh2Event.st1Ph1ToSt1Ph2;
    final int st1Ph2ToSt2Ph1 = QueuePhPh2Event.st1Ph2ToSt2Ph1;
    final int st1Ph2ToSt2Ph2 = QueuePhPh2Event.st1Ph2ToSt2Ph2;
    final int st2Ph1ToSt2Ph1 = QueuePhPh2Event.st2Ph1ToSt2Ph1;
    final int st2Ph1ToSt2Ph2 = QueuePhPh2Event.st2Ph1ToSt2Ph2;
    final int st2Ph2ToSt2Ph1 = QueuePhPh2Event.st2Ph2ToSt2Ph1;
    final int st2Ph2ToSt2Ph2 = QueuePhPh2Event.st2Ph2ToSt2Ph2;

    private double lambda1;
    private double lambda2;
    private double s1Mu1;
    private double s1Mu2;
    private double s2Mu1;
    private double s2Mu2;
    private double prob1;// Probability of type 1 costumers
    private int cptyServ2;
    
    /**
     * Constructor for the PhPh2 class
     * @param lambda1 Arrival rate of the type 1 customers
     * @param lambda2 Arrival rate of the type 2 customers
     * @param s1Mu1 Service rate of type 1 customers in server 1 
     * @param s1Mu2 Service rate of type 1 customers in server 2
     * @param s2Mu1 Service rate of type 2 customers in server 1
     * @param s2Mu2 Service rate of type 2 customers in server 2
     * @param cptyServ2
     */
    public QueuePhPh2(double lambda1, double lambda2, double s1Mu1,
            double s1Mu2, double s2Mu1, double s2Mu2, int cptyServ2) {
        super(new QueuePhPh2State(0, 0, 0), QueuePhPh2Event.getAllEvents());
        this.lambda1 = lambda1;
        this.lambda2 = lambda2;
        this.s1Mu1 = s1Mu1;
        this.s1Mu2 = s1Mu2;
        this.s2Mu1 = s2Mu1;
        this.s2Mu2 = s2Mu2;
        this.prob1 = lambda1 / (lambda1 + lambda2);
        this.cptyServ2 = cptyServ2;
    }

    /**
     * Used by GUI
     */
    public QueuePhPh2() {
        this(1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2);
    }

    /*
     * (non-Javadoc)
     * 
     * @see jmarkov.GeomProcess#active(jmarkov.State, int)
     */
    @Override
    public boolean active(QueuePhPh2State i, int absLevel, QueuePhPh2Event e) {

        boolean result = false;

        switch (e.type) {
        case arrival:
            result = true;
            break;
        case st1Ph1ToSt1Ph2:
            result = (i.getServer1Status() == 1);
            break;
        case st1Ph2ToSt2Ph1:
            result = (i.getServer1Status() == 2);
            break;
        case st1Ph2ToSt2Ph2:
            result = (i.getServer1Status() == 2);
            break;
        case st2Ph1ToSt2Ph1:
            result = (i.getServer2Status() == 1);
            break;
        case st2Ph1ToSt2Ph2:
            result = (i.getServer2Status() == 1);
            break;
        case st2Ph2ToSt2Ph1:
            result = (i.getServer2Status() == 2);
            break;
        case st2Ph2ToSt2Ph2:
            result = (i.getServer2Status() == 2);
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
    public GeomRelState<QueuePhPh2State>[] dests(QueuePhPh2State i,
            int absLevel, QueuePhPh2Event e) {

        int newCostumersInServer2 = i.getCostumersInServer2();
        int newServer1Status = i.getServer1Status();
        int newServer2Status = i.getServer2Status();
        int rLevel = 0;

        if (e.type == arrival) {
            rLevel = 1;
            if (absLevel == 0)
                newServer1Status = 1;
        } else if (e.type == st1Ph1ToSt1Ph2) {
            newServer1Status++;
        } else if (e.type == st1Ph2ToSt2Ph1 || e.type == st1Ph2ToSt2Ph2) {
            if (newCostumersInServer2 < cptyServ2) {
                if (newCostumersInServer2 == 0) {
                    switch (e.type) {
                    case st1Ph2ToSt2Ph1:
                        newServer2Status = 1;
                        break;
                    case st1Ph2ToSt2Ph2:
                        newServer2Status = 2;
                    }
                }

                if (absLevel == 1) {
                    newServer1Status = 0;
                } else {
                    newServer1Status = 1;
                    rLevel = -1;
                }
                newCostumersInServer2++;
            } else {
                newServer1Status = 3;
            }
        } else if (e.type == st2Ph1ToSt2Ph1 || e.type == st2Ph1ToSt2Ph2 || //
                e.type == st2Ph2ToSt2Ph1 || e.type == st2Ph2ToSt2Ph2) {
            if (newCostumersInServer2 == 1) {
                newServer2Status = 0;
                newCostumersInServer2--;
            } else {
                if (e.type == st2Ph1ToSt2Ph1 || e.type == st2Ph2ToSt2Ph1) {
                    newServer2Status = 1;
                } else {
                    newServer2Status = 2;
                }
                if (newServer1Status == 3 && absLevel == 1) {
                    newServer1Status = 0;
                } else if (newServer1Status == 3) {
                    rLevel = -1;
                    newServer1Status = 1;
                } else {
                    newCostumersInServer2--;
                }
            }
        }

        QueuePhPh2State newSubState = new QueuePhPh2State(
                newCostumersInServer2,//
                newServer1Status, newServer2Status);
        GeomRelState<QueuePhPh2State> s;
        if ((((e.type == st1Ph2ToSt2Ph1 || e.type == st1Ph2ToSt2Ph2)//
                && i.getCostumersInServer2() < cptyServ2) && absLevel == 1)
                || //
                /**
                 * Verification for boundary destination when finish phase 2 in
                 * station 1
                 */
                (e.type == st2Ph1ToSt2Ph1 || e.type == st2Ph1ToSt2Ph2
                        || e.type == st2Ph2ToSt2Ph1 //
                || e.type == st2Ph2ToSt2Ph2) && i.getCostumersInServer2() != 1
                && i.getServer1Status() == 3 && absLevel == 1
        /**
         * Verification for boundary destination when finish a service phase in
         * station 1
         */
        ) {
            s = new GeomRelState<QueuePhPh2State>(newSubState);
        } else {
            s = new GeomRelState<QueuePhPh2State>(newSubState, rLevel);
        }
        StatesSet<GeomRelState<QueuePhPh2State>> statesSet //
        = new StatesSet<GeomRelState<QueuePhPh2State>>(s);
        return statesSet.toStateArray();
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
    public double rate(QueuePhPh2State i, int iLevel, QueuePhPh2State j,
            int jLevel, QueuePhPh2Event e) {

        double tasa = 0;

        switch (e.type) {
        case arrival:
            tasa = lambda1 + lambda2;
            break;
        case st1Ph1ToSt1Ph2:
            tasa = s1Mu1;
            break;
        case st1Ph2ToSt2Ph1:
            tasa = prob1 * s1Mu2;
            break;
        case st1Ph2ToSt2Ph2:
            tasa = (1.0 - prob1) * s1Mu2;
            break;
        case st2Ph1ToSt2Ph1:
            tasa = (prob1) * s2Mu1;
            break;
        case st2Ph1ToSt2Ph2:
            tasa = (1.0 - prob1) * s2Mu1;
            break;
        case st2Ph2ToSt2Ph1:
            tasa = (prob1) * s2Mu2;
            break;
        case st2Ph2ToSt2Ph2:
            tasa = (1.0 - prob1) * s2Mu2;
            break;
        }
        return tasa;
    }// end of rate

    /*
     * (non-Javadoc)
     * 
     * @see jmarkov.SimpleMarkovProcess#description()
     */
    @Override
    public String description() {
        return "Ph/Ph/2";
    }// end of description.

    /**
     * Main method
     * @param a Not used
     */
    public static void main(String[] a) {

        Date inicio = new Date();
        QueuePhPh2 theQueue = new QueuePhPh2(2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 2);
        theQueue.showGUI();
        theQueue.printAll();
        Date termina = new Date();

        long diferencia = termina.getTime() - inicio.getTime();

        System.out.println("tiempo de ejecución = " + diferencia);

    }// end of main

}// end of class

/**
 * @author Julio Goez
 * 
 */
class QueuePhPh2State extends PropertiesState {

    /**
     * We identify the states with a vector with the number of costumers in
     * server 2 (0,1,2,3,...), the phase of server in station 1 (1,2), 0 if idle
     * or 3 if blocked, and the phase of server in station 2, 1 if items of type
     * 1, 2 if items of type 2 or 0 if idle.
     * 
     * @param CostumersInServer2
     *            Costumers in station 2.
     * @param Service1Phase
     *            Service phase in server of station 1.
     * @param Service2Phase
     *            Service phase in server of station 1.
     */
    public QueuePhPh2State(int CostumersInServer2, int Service1Phase,
            int Service2Phase) {

        super(new int[] { CostumersInServer2, Service1Phase, Service2Phase });

    }

    @Override
    public void computeMOPs(MarkovProcess mp) {
        setMOP(mp,"Server 1 Utilization", (getServer1Status() != 0) ? 1 : 0);
        setMOP(mp,"Server 2 Utilization", (getServer1Status() != 0) ? 1 : 0);
    }

    /**
     * Returns the status of the first Server 
     * @return Status of the first Server
     */
    public int getServer1Status() {
        return this.prop[1];
    }

    /**
     * Returns the status of the second Server 
     * @return Status of the second Server
     */
    public int getServer2Status() {
        return this.prop[2];
    }

    /**
     * Returns the status of the second Server buffer 
     * @return Status of the second Server buffer
     */
    public int getBuffer2Size() {
        return (this.prop[0] - 1);
    }
    
    /**
     * Returns the number of customers in the second Server buffer 
     * @return Number of customers in the second Server buffer
     */
    public int getCostumersInServer2() {
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

    @Override
    public String label() {
        String stg = "St2:" + this.getCostumersInServer2() + ", " + "Ph1:"
                + this.getServer1Status() + ", " + "Ph2:"
                + this.getServer2Status();
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

/**
 * This class define the events.
 * 
 * @author Julio Goez
 * 
 * To change the template for this generated type comment go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
class QueuePhPh2Event extends Event {

    /**
     * Arrivals to the station 1.
     */
    final static int arrival = 0; // Arrivals to the system.

    /**
     * Finished phase 1 in server 1 and go to phase 2 in server 1.
     */
    final static int st1Ph1ToSt1Ph2 = 1; // Finished phase 1 of server 1.

    /**
     * Finished phase 2 in server 1 and go to phase 1 in server 2.
     */
    final static int st1Ph2ToSt2Ph1 = 2;

    /**
     * Finished phase 2 in server 1 and go to phase 2 in server 2.
     */
    final static int st1Ph2ToSt2Ph2 = 3;

    /**
     * Finished phase 1 in server 2 and go to phase 1 in server 2.
     */
    final static int st2Ph1ToSt2Ph1 = 4; // Finished phase 1 of server 2.

    /**
     * Finished phase 1 in server 2 and go to phase 2 in server 2.
     */
    final static int st2Ph1ToSt2Ph2 = 5; // Finished phase 1 of server 2.

    /**
     * Finished phase 2 in server 2 and go to phase 1 in server 2.
     */
    final static int st2Ph2ToSt2Ph1 = 6; // Finished phase 2 of server 2.

    /**
     * Finished phase 2 in server 2 and go to phase 2 in server 2.
     */
    final static int st2Ph2ToSt2Ph2 = 7; // Finished phase 2 of server 2.

    int type;

    QueuePhPh2Event(int type) {
        this.type = type;
    }

    static EventsSet<QueuePhPh2Event> getAllEvents() {
        EventsSet<QueuePhPh2Event> E = new EventsSet<QueuePhPh2Event>();
        E.add(new QueuePhPh2Event(arrival));
        E.add(new QueuePhPh2Event(st1Ph1ToSt1Ph2));
        E.add(new QueuePhPh2Event(st1Ph2ToSt2Ph1));
        E.add(new QueuePhPh2Event(st1Ph2ToSt2Ph2));
        E.add(new QueuePhPh2Event(st2Ph1ToSt2Ph1));
        E.add(new QueuePhPh2Event(st2Ph1ToSt2Ph2));
        E.add(new QueuePhPh2Event(st2Ph2ToSt2Ph1));
        E.add(new QueuePhPh2Event(st2Ph2ToSt2Ph2));
        return E;

    }

    @Override
    public String label() {
        String stg = "";
        switch (type) {
        case (arrival):
            stg = "Arrival";
            break;
        case (st1Ph1ToSt1Ph2):
            stg = "S1 first stage";
            break;
        case (st1Ph2ToSt2Ph1):
            stg = "S1 second stage to S2 first stage";
            break;
        case (st1Ph2ToSt2Ph2):
            stg = "S1 second stage to S2 second stage";
            break;
        case (st2Ph1ToSt2Ph1):
            stg = "S2 first stage to S2 first stage";
            break;
        case (st2Ph1ToSt2Ph2):
            stg = "S2 first stage to S2 second stage";
            break;
        case (st2Ph2ToSt2Ph1):
            stg = "S2 second stage to S2 fisrt stage";
            break;
        case (st2Ph2ToSt2Ph2):
            stg = "S2 second stage to S2 second stage";
            break;
        }
        return stg;
    }
}
