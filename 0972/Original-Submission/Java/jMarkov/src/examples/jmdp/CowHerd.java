package examples.jmdp;

import jmarkov.MarkovProcess;
import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.FiniteDP;
import jmarkov.jmdp.solvers.FiniteSolver;

/**
 * The cow herd problem. Given the prices on a horizon decide whether keep or
 * sell the cow herd. There is a max number of cows per period. If you keep the
 * cows there will be about 1.5 times the current number of cows in the next
 * period.
 * 
 * @author Germán Riaño
 */
public class CowHerd extends FiniteDP<CowQuantity, CowActions> {

	private double price[] = null;
	private int mxCows[] = null;

	private ActionsSet<CowActions> actions = new ActionsSet<CowActions>();

	/**
	 * @param initial
	 * @param horizon
	 * @param price
	 * @param mxCows
	 */
	public CowHerd(States<CowQuantity> initial, int horizon, double[] price,
			int[] mxCows) {
		super(initial, horizon);
		this.price = price;
		this.mxCows = mxCows;
		// initStates();
		actions.add(CowActions.SELL);
		actions.add(CowActions.KEEP);
	}

	/**
	 * Builds a cow Herd problem with the given initial cows.
	 * @param initialQuantity
	 * @param horizon
	 * @param price
	 * @param mxCows
	 */
	public CowHerd(int initialQuantity, int horizon, double[] price,
			int[] mxCows) {
		this(initStates(initialQuantity), horizon, price, mxCows);
	}

	private static States<CowQuantity> initStates(int quantity) {
		CowQuantity state = new CowQuantity(quantity);
		return new StatesSet<CowQuantity>(state);
	}
	
	/**
	 * @param initQ
	 * @return The optimal value fir this initial quantity
	 * @throws SolverException
	 */
	public double getValue(int initQ) throws SolverException{
		ValueFunction<CowQuantity> vf = getOptimalValueFunction();
		CowQuantity state = new CowQuantity(initQ);
		return vf.get(state);
	}

	// public void initStates() {
	// states = new StatesCollection<CowQuantity>();
	// int maxCows = mxCows[5];
	// for (int i = 0; i <= maxCows; i++) {
	// states.add(new CowQuantity(i));
	// }
	// }

	@Override
	public double immediateCost(CowQuantity i, CowActions a, int t) {
		if (a == CowActions.SELL)
			return (-i.cows * price[t]);
		return 0.0;
	}

	@Override
	public double finalCost(CowQuantity i) {
		return (-i.cows * price[5]);
	}

	/**
	 * @see jmarkov.jmdp.FiniteDP#destination(State , Action , int)
	 */
	@Override
	public CowQuantity destination(CowQuantity i, CowActions a, int t) {
		if (a == CowActions.KEEP) {
			int newCows = (int) Math.floor(1.5 * i.cows);
			if (newCows <= mxCows[t + 1]) {
				return new CowQuantity(newCows);
			}
		} else if (a == CowActions.SELL) {
			return new CowQuantity(0);
		}
		return null;
	}

	@Override
	public Actions<CowActions> feasibleActions(CowQuantity st, int t) {
		if (st.cows == 0) {
			ActionsSet<CowActions> act = new ActionsSet<CowActions>();
			act.add(CowActions.KEEP);
			return act;
		} else if ((int) Math.floor(1.5 * st.cows) > mxCows[t + 1]) {
			ActionsSet<CowActions> act = new ActionsSet<CowActions>();
			act.add(CowActions.SELL);
			return act;
		} else
			return actions;
	}

	/**
     * Small test program
	 * @param a Not used
	 * @throws SolverException
	 */
	public static void main(String a[]) throws SolverException {
		double price[] = { 0, 10, 15, 25, 7, 10 };
		int mxCows[] = { 50, 75, 112, 168, 252, 378 };
		CowQuantity initial = new CowQuantity(50);
		CowHerd prob = new CowHerd(50, 5, price, mxCows);
		FiniteSolver<CowQuantity, CowActions> theSolver = new FiniteSolver<CowQuantity, CowActions>(
				prob);
		theSolver.solve();
		System.out.println(theSolver.bestPolicy(initial));
	}

}

class CowQuantity extends State {
	int cows; // number of cows

	CowQuantity(int cws) {
		cows = cws;
	}

	@Override
    public int compareTo(State s) {
		CowQuantity i = (CowQuantity) s;
		return cows - i.cows;
	}

	@Override
	public String label() {
		return cows + " cows";
	}

   

    /**
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        // TODO Auto-generated method stub
        return true;
    }

    /**
     * @see jmarkov.basic.State#computeMOPs(jmarkov.MarkovProcess)
     */
    @Override
    public void computeMOPs(MarkovProcess<?, ?> model) {
        // TODO Auto-generated method stub
        
    }

}

class CowActions extends Action {
	/**   Sell the cows     */
	public static final CowActions SELL = new CowActions("Sell");
	/**   Keep all the cows     */
	public static final CowActions KEEP = new CowActions("Keep");
	String name;
    
	private CowActions(String str) {
        name = str;
	}
    /**
     * @see jmarkov.basic.Action#label()
     */
    @Override
    public String label() {
        return name;
    }

    /**
     * @see java.lang.Comparable#compareTo(Object)
     */
    public int compareTo(Action o) {
        return name.compareTo(o.label());
    }
    
}
