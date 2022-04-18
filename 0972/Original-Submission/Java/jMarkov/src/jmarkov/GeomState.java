/**
 * GeomState.java
 * Created: Jun 26, 2005
 */
package jmarkov;

import jmarkov.basic.State;

/**
 * The actual Geometric model is build using this class. The user normally does
 * not have to manipulate this class.
 * 
 * @author Julio Goez, Germán Riano. Universidad de los Andes. (C) 2005
 * @param <Sub>
 *            The sub-States class used.
 */
public final class GeomState<Sub extends State> extends State {

    /**
     * This represents the relative level.
     */
    protected int level;

    /**
     * subState represnts the background states in every level.
     */
    protected Sub subState;

    /**
     * Creates a GeomState with the given level, ans subState.
     * 
     * @param level
     * @param subState
     */
    public GeomState(Sub subState, int level) {
        super();
        this.level = level;
        this.subState = subState;
    }

    /**
     * @return Returns the level.
     */
    public int getLevel() {
        return level;
    }

    /**
     * @return true if this state is level 0.
     */
    public boolean isBoundary() {
        return (getLevel() == 0);
    }

    /**
     * @return Returns the subState.
     */
    public Sub getSubState() {
        return subState;
    }

    /**
     * Compares GeomStates according to level first and then according to the
     * subStates comparator.
     * 
     * @param s
     *            state to compare to.
     */
    @SuppressWarnings("unchecked")
    @Override
    public int compareTo(State s) {

        GeomState gs = (GeomState) s;
        // First Compare level
        if (level < gs.getLevel()) {
            return -1;
        } else if (level > gs.getLevel()) {
            return +1;
        }
        // compare the rest
        else {
            return getSubState().compareTo(gs.getSubState());
        }
    }

    // *** MOPs stuff

//    /**
//     * @see jmarkov.basic.State#addMOPName(java.lang.String)
//     */
//    public @Override boolean addMOPName(String mopName) {
//        return subState.addMOPName(mopName);
//    }
//
//    /**
//     * @see jmarkov.basic.State#clearMOPs()
//     */
//    @Override
//    public void clearMOPs() {
//        subState.clearMOPs();
//    }

    /**
     * @see jmarkov.basic.State#computeMOPs(MarkovProcess)
     */
    @Override
    public void computeMOPs(MarkovProcess mp) {
        subState.computeMOPs(mp);
    }

    /**
     * @see jmarkov.basic.State#getMOP(int)
     */
    @Override
    public double getMOP(int index) {
        return subState.getMOP(index);
    }

    

//    /**
//     * @see jmarkov.basic.State#getMOPIndex(java.lang.String)
//     */
//    public @Override int getMOPIndex(String name) {
//        return subState.getMOPIndex(name);
//    }
//
//    /**
//     * @see jmarkov.basic.State#getMOPNames()
//     */
//    public @Override String[] getMOPNames() {
//        return subState.getMOPNames();
//    }

//    /**
//     * @see jmarkov.basic.State#getMOPNames(int)
//     */
//    public @Override String getMOPNames(int i) {
//        return subState.getMOPNames(i);
//    }

    /**
     * @see jmarkov.basic.State#label()
     */
    @Override
    public String label() {
        return "L:" + level + "(" + subState.label() + ")";
    }

//    /**
//     * @see jmarkov.basic.State#mopsNames()
//     */
//    public @Override ArrayList<String> mopsNames() {
//        return subState.mopsNames();
//    }
//
//    /**
//     * @see jmarkov.basic.State#numMOPNames()
//     */
//    public @Override int numMOPNames() {
//        return subState.numMOPNames();
//    }

//    /**
//     * @see jmarkov.basic.State#setMOP(int, double)
//     */
//    @Override
//    public int setMOP(int index, double value) {
//        return subState.setMOP(index, value);
//    }

    /**
     * @see jmarkov.basic.State#setMOP(MarkovProcess,java.lang.String,  double)
     */
    @Override
    public int setMOP(MarkovProcess mp,String mopName, double value) {
        return subState.setMOP( mp,mopName,  value);
    }

//    /**
//     * @see jmarkov.basic.State#setMOPNames(java.lang.String[])
//     */
//    public @Override void setMOPNames(String[] mopNames) {
//        subState.setMOPNames(mopNames);
//    }

    /*
     * (non-Javadoc)
     * 
     * @see jmarkov.State#description()
     */
    @Override
    public String description() {
        return "Level: " + level + ", sub-State: " + subState.description();
    }

    /**
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        return level >= 0 && subState.isConsistent();
    }

}
