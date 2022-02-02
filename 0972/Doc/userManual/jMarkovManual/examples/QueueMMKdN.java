package examples.jmarkov;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import jmarkov.MarkovProcess;
import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.EventsSet;
import jmarkov.basic.PropertiesState;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;

/**
 * This class represents is a system with K different
 * exponential servers with rates mu1, mu2, etc, 
 * respectively, and arrival rate lambda. A customer
 * that finds more then one server idle chooses according 
 * to relative intensities
 * <tex txt="$\alpha_1, \alpha_2, \ldots, \alpha_K$">
 * alpha1, alpha2, etc</tex>. The probability of choosing 
 * idle server k will be given by
 * <tex txt="\[\beta_k = \frac{\alpha_k}{\sum_{\ell\in \cal I} \alpha_{\ell}},\] 
 * where $\cal I$ is the set of idle servers.">  
 * alpha(k) / sum( alpha(j)), where the sum is over the set of idle servers.
 * </tex> 
 * @author Germán Riaño. Universidad de los Andes.
 */
public class QueueMMKdN extends SimpleMarkovProcess<QueueMMKdNState,QueueMMKdNEvent> {
	// Events

	private double lambda;
	private double[] mu, alpha;
	private int K; // number of servers
	private int N;
	private static final int NDARRIVAL = QueueMMKdNEvent.NDARRIVAL;
	private static final int DIRARRIVAL = QueueMMKdNEvent.DIRARRIVAL;
	private static final int DEPARTURE = QueueMMKdNEvent.DEPARTURE;

	/**
	 * Constructs a M/M/Kd queue with arrival rate lambda and service 
	 * rates mu, relative probabilities of choosing each server alpha 
	 * @param lambda  Arrival rate
	 * @param mu      Server   rates
	 * @param alpha   Relative probability of an arriving customer choosing each server. 
	 * @param N       Max number in the system
	 */
	public QueueMMKdN(double lambda, double[] mu, double[] alpha, int N) {
		super(
			new QueueMMKdNState(mu.length, alpha),
			QueueMMKdNEvent.getAllEvents(mu.length));
		this.K = mu.length;
		this.lambda = lambda;
		this.mu = mu;
		this.alpha = alpha;
		this.N = N;
	}
	
    /**
     * Returns an QueueMMKdN object with arrival rate 1.0, 
     * service rates of 2.0, 3.0 and 4.0;  
     * and capacity of 8 customers in the system. 
     * Used by GUI
     */
	public QueueMMKdN(){
		this(1.0, new double[]{2,3,4}, new double[]{2,3,4}, 8);
	}

	/**
	 * Determines the active events. 
	 */
	@Override
	public boolean active(QueueMMKdNState i, QueueMMKdNEvent e) {
		boolean result = false;
		switch (e.type) {
			case (NDARRIVAL) : // NDARIIVAL occurs only if servers are busy and there is roon in the Q
				result = (i.allBusy() && (i.getQSize() < N - K));
				break;
			case (DIRARRIVAL) :
				{
					result = (i.getStatus(e.server) == 0);
					//DirARRIVAL occurs if server is EMPTY.
					break;
				}
			case (DEPARTURE) :
				{ // ev.type == DEPARTURE     
					result = (i.getStatus(e.server) == 1);
					//DEPARTURE occurs if server is busy.
				}
		}
		return result;
	}

	/* 
	 * Determines the possible destination event (actually one in this case).
	 */

	@Override
	public States<QueueMMKdNState> dests(QueueMMKdNState i, QueueMMKdNEvent e) {
		int[] status = new int[K];
		for (int k = 0; k < K; k++)
			status[k] = i.getStatus(k); //copy current values
		int Q = i.getQSize();
		switch (e.type) {
			case (NDARRIVAL) :
				Q++; // non-directed ARRIVAL
				break;
			case (DIRARRIVAL) :
				status[e.server] = 1; //directed ARRIVAL, picks a server.
				break;
			case (DEPARTURE) :
				if (Q > 0) { //there is Queue
					status[e.server] = 1; //set (keeps) server busy
					Q--; // reduce queue
				} else
					status[e.server] = 0; //set server idle
		}
		return new StatesSet<QueueMMKdNState>(new QueueMMKdNState(status, Q, alpha));
	}

	/* 
	 * The rate is lambda, or mu for non-directed arrival and for departure. 
	 * For directed arrival rate id lambda 8 prob(server is choosen)
	 * @see jmarkov.SimpleMarkovProcess#rate(jmarkov.State, jmarkov.State, jmarkov.Event)
	 */
	@Override
	public double rate(QueueMMKdNState i, QueueMMKdNState j, QueueMMKdNEvent e) {
		double result = 0;

		switch (e.type) {
			case (DEPARTURE) :
				result = mu[e.server];
				break;
			case (NDARRIVAL) :
				result = lambda;
				break; //non-directed arrival
			case (DIRARRIVAL) :
				result = i.prob(e.server) * lambda;
		}
		return result;
	}

	/**
	 * Main Method. This asks the user for parameters 
     * and tests the program.
     * @param a Not used 
	 */
	public static void main(String[] a) {
		BufferedReader rdr =
			new BufferedReader(new InputStreamReader(System.in));
		try {
			System.out.println("Input Rate: ");
			double lda = Double.parseDouble(rdr.readLine());
			System.out.println("Num Servers: ");
			int K = Integer.parseInt(rdr.readLine());
			double mu[] = new double[K];
			double alpha[] = new double[K];
			for (int k = 0; k < K; k++) {
				System.out.println("Service rate, server " + (k + 1) + " : ");
				mu[k] = Double.parseDouble(rdr.readLine());
			}
			for (int k = 0; k < K; k++) {
				System.out.println(
					"Choosing intensity, server  " + (k + 1) + " : ");
				alpha[k] = Double.parseDouble(rdr.readLine());
			}
			System.out.println("Max in system : ");
			int N = Integer.parseInt(rdr.readLine());
			QueueMMKdN theModel = new QueueMMKdN(lda, mu, alpha, N);
			theModel.showGUI();
			//theModel.setDebugLevel(2);
			theModel.printAll();
		} catch (IOException e) {
		};
	}

	/**
	 * @see jmarkov.SimpleMarkovProcess#description()
	 */
	@Override
	public String description() {
		String stg = "M/M/k/N SYSTEM\n\n";
		stg += "Multiple server queue with " + this.K + " different servers\n";
		stg += "Arrival Rate = " + lambda + ", Max number in system " + N;
		return stg;
	}

} //class end
/**
 * This is a particular case of propertiesState, whith K + 1
 * properties, namely the server 1, 2, ..., K status, plus the queue level.
 * 
 * @author Germán Riaño. Universidad de los Andes.
 */
class QueueMMKdNState extends PropertiesState {

	private int K; // number of servers
	private double sumProb = -1; // sum of relative probabilities
	private double[] alpha; //relative frequency of servers 
	private double[] beta; //probabilities for this state 
	/** 
	 * Constructs a state for an empty system with K servers, and
	 * choosing intensities alpha.
	 * @param K Number of servers.
	 */
	QueueMMKdNState(int K, double[] alpha) {
		this(new int[K], 0, alpha);
	}

	/**
	 * We identify each State with a vector that counts the 
	 * ststus fo the k servers and 
	 * the number in queue (0,1, ..,N-K). 
	 */
	QueueMMKdNState(int[] status, int Qsize, double[] alpha) {
		super(alpha.length + 1);
		this.K = alpha.length;
		this.alpha = alpha;
		this.beta = new double[K];
		int sum = 0; // adds the number of busy server = people in service
		for (int i = 0; i < K; i++) {
			prop[i] = status[i];
			sum += status[i];
		}
		prop[K] = Qsize;
	}

	/**
	 * Computes the MOPs
	 * @see jmarkov.basic.State#computeMOPs(MarkovProcess)
	 */
	@Override
	public void computeMOPs(MarkovProcess mp) {
		double sum = 0.0;
		for (int i = 0; i < K; i++) {
			sum += getStatus(i);
			setMOP(mp,"Server Status " + (i + 1), getStatus(i));
		}
		setMOP(mp,"Queue Length", getQSize());
		setMOP(mp,"Number in System", sum + getQSize());
	}

    /**
     * Returns the status of the kth Server 
     * @param k server index
     * @return Status of the kth Server
     */
	public int getStatus(int k) {
		return prop[k];
	}

    /**
     * Returns the size of the queue 
     * @return Status of the size of the queue
     */
	public int getQSize() {
		return prop[K];
	}
	/**
	 * Determines if all servers are busy
     * @return True, if all servers are busy. False, otherwise
	 */
	public boolean allBusy() {
		boolean result = true;
		for (int k = 0; result && (k < K); k++)
			result = result && (getStatus(k) == 1);
		return result;
	}
    /**
     * Determines if all servers are idle
     * @return True, if all servers are idle. False, otherwise
     */
	public boolean allIdle() {
		boolean result = true;
		for (int k = 0; result && (k < K); k++)
			result = result && (getStatus(k) == 0);
		return result;
	}
    /**
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        // TODO Complete
        return true;
    }
	/*
	 * determines the sum of all intensities for idle servers. The result
	 * is kept in sumProb for future use.
	 */
	private double sum() {
		if (sumProb != -1)
			return sumProb;
		double res = 0;
		for (int k = 0; k < K; k++) {
			res += (1 - getStatus(k)) * alpha[k];
		}
		return (sumProb = res);
	}
	/**
	 * Detemines the probability of an idle server being choosen 
	 * among idle servers. A customer that finds more then one server
	 *  idle chooses according to relative intensities
	 * <tex txt="$\alpha_1, \alpha_2, \ldots, \alpha_K$">
	 * alpha1, alpha2, etc</tex>. The probability of choosing idle
	 *  server k will be given by
	 * <tex txt="\[\beta_k = \frac{\alpha_k}{\sum_{\ell\in \cal I} \alpha_{\ell}},\] 
	 * where $\cal I$ is the set of idle servers.">  
	 * alpha(k) / sum(j, alpha(j)), where the sum is over the set 
	 * of idle servers. </tex> 
     * @param server server index
     * @return probability of an idle server being choosen 
     * among idle servers
	 */
	public double prob(int server) {
		if (beta[server] != 0)
			return beta[server];
		return (
			beta[server] = (((1 - getStatus(server)) * alpha[server]) / sum()));
	}

	/**
	 * Returns a label with the format SxxQz, whre xx is the list of busy servers.
	 * @see jmarkov.basic.State#label()
	 */
	@Override
	public String label() {
		String stg = "S";
		for (int k = 0; k < K; k++) {
			stg += (getStatus(k) == 1) ? "" + (k + 1) : "";
		}
		return stg + "Q" + getQSize();
	}

	/*
	 * This method gives a verbal description of the State.
	 */
	@Override
	public String description() {
		String stg = "";
		if (!allIdle())
			stg += "Busy Servers:";
		else
			stg += "No one in service";
		for (int k = 0; k < K; k++) {
			stg += (getStatus(k) == 1) ? "" + (k + 1) + "," : "";
		}
		stg += " There are " + getQSize() + " customers waiting in queue.";
		return stg;
	}

}
/**
	 * 
	 * This class define the events. 
	 * An event has two components: type which can have three values 
	 * depending whether it represents a directed arrival, a 
	 * non-directed arrival or a departure, and server, which 
	 * represents the choosen server (if arrival) or the finishing 
	 * server. For non-directed arrivals we set server -1 by convention.
	 * 
	 * @author Germán Riaño
	 *
	 */
class QueueMMKdNEvent extends Event {
	final static int NDARRIVAL = 0;
	//Non directed arrival (when all servers are busy)
	final static int DIRARRIVAL = 1; //Directed arrival chooses among server(s)
	final static int DEPARTURE = 2;
	int type; // ARRIVAL or DEPARTURE
	/*	server = chosen server if ARRIVAL finds many available,
	 * server = -1 if no server available 
	 * server = finishing server if DEPARTURE event
	 */
	int server;
	QueueMMKdNEvent(int type, int server) {
		this.type = type;
		this.server = server;
	}

	static EventsSet<QueueMMKdNEvent> getAllEvents(int K) {
		EventsSet<QueueMMKdNEvent> eSet = new EventsSet<QueueMMKdNEvent>();
		eSet.add(new QueueMMKdNEvent(NDARRIVAL, -1));
		for (int i = 0; i < K; i++) {
			eSet.add(new QueueMMKdNEvent(DIRARRIVAL, i));
		}
		for (int i = 0; i < K; i++) {
			eSet.add(new QueueMMKdNEvent(DEPARTURE, i));
		}
		return eSet;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String label() {
		String stg = "";
		switch (type) {
			case (NDARRIVAL) :
				stg += "Non-directed arrival";
				break;
			case (DIRARRIVAL) :
				stg += "Directed arrival to server " + (server + 1);
				break;
			case (DEPARTURE) :
				stg += "Departure from server " + (server + 1);
				break;
		}
		return stg;
	}

} //end class
// Lets start defining the State

// Now we define the main  class

