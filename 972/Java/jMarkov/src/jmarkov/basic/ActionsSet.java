package jmarkov.basic;

import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

/**
 * This class represents a set of objects Action. It uses the TreeSet structure
 * to avoid repeated actions. This class extends the Actions class.
 * 
 * @author Andres Sarmiento, Germán Riaño - Universidad de Los Andes
 * @param <A> The action class.
 * 
 * @see java.util.Collection
 * @see jmarkov.basic.Actions
 * @see jmarkov.basic.Action
 */

public class ActionsSet<A extends Action> implements Actions<A> {

	private Collection<A> actions = null;
	
	/**
	 * Creates a set of Actions from a given set of Actions.
	 * 
	 * @param ac
	 *            set of Actions of type Actions.
	 */
	public ActionsSet(Actions<A> ac) {
		this.actions = new TreeSet<A>();
		for (A a : ac)
			actions.add(a);
	}

    /**
     * Creates a set of actions from any iterable object over actions.
     * @param actIter
     */
	public ActionsSet(Iterable<A> actIter) {
		this.actions = new TreeSet<A>();
		for (A a : actIter)
			actions.add(a);
	}

	/**
	 * Creates a set of Actions from a given array of Actions. This constructor
	 * organizes the actions in a TreeSet.
	 * 
	 * @param acArray  set of Actions of type Actions.
	 */
	public ActionsSet(A[] acArray) {
		this.actions = new TreeSet<A>();
		for (A a : acArray) {
			actions.add(a);
		}
	}
	
	/**
	 * Creates a set of Actions from a given Action. This constructor
	 * organizes the actions in a TreeSet.
	 * 
	 * @param ac  an Action.
	 */
	public ActionsSet(A ac) {
		this.actions = new TreeSet<A>();
		actions.add(ac);		
	}

	/**
	 * Creates an empty set of Actions.
	 */
	public ActionsSet() {
		this.actions = new TreeSet<A>();
	}

	/**
	 * Creates a set of Actions from a given collection of Actions. This
	 * constructor organizes the actions in a TreeSet.
	 * 
	 * @param col set of Actions of type Actions.
	 *//*
	public ActionsSet(Collection<? extends Action> col) {
		this.actions = new TreeSet<Action>();
		for (Action a : col) {
			actions.add(a);
		}
	}
*/
	/**
	 * This method adds a new action to the set.
	 * @param a The action to be added.
	 */
	public void add(A a){
		actions.add(a);
	}
	
	/**
	 * This method returns a safe way to walk along the actions in a particular
	 * set. Collections and their implementations (Set, List, and Map) have
	 * iterators defined by default.
	 * 
	 * @return iterator over the states.
	 */
	public final Iterator<A> iterator() {
		return actions.iterator();
	}

	/**
	 * @return the amount of actions in the set.
	 */
	public int size() {
		return actions.size();
	}
	
    @Override
    public String toString(){
        String stg = "(";
        for (A a : this) stg += a + " ";
        return stg + ")";
    }


}