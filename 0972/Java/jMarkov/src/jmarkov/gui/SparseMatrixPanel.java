/*
 * Created on 31/08/2004
 */
package jmarkov.gui;

import java.awt.BorderLayout;
import java.awt.Component;

import javax.swing.ImageIcon;
import javax.swing.JEditorPane;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTree;
import javax.swing.ToolTipManager;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

import jmarkov.MarkovProcess;
import jmarkov.SimpleMarkovProcess;
import jmarkov.basic.Event;
import jmarkov.basic.State;
import jmarkov.basic.Transition;
import jmarkov.basic.Transitions;
import jmarkov.basic.exceptions.NotUnichainException;

/**
 * This class show with a tree all states / actions / Events
 * 
 * @author Germán Riaño. Universidad de los Andes.
 */
public class SparseMatrixPanel extends InfoPanel implements
        TreeSelectionListener {

    private static final long serialVersionUID = 1969;

    private JTree jTree = null;
    private DefaultTreeModel treeModel = null;
    private JScrollPane infoScrollPane = null;
    private boolean showEvents = false; // whether events are show in the tree
    private DefaultMutableTreeNode root = null;
    private JSplitPane jSplitPane = null;
    private JScrollPane treeScrollPane = null; // @jve:decl-index=0:visual-constraint="50,-164"
    private JEditorPane treeInfoPanel = null;
    private ImageIcon stateIcon = null;
    private ImageIcon eventIcon = null;
    private ImageIcon modelIcon = null;

    /**
     * This method initializes the view.
     */
    public SparseMatrixPanel() {
        super();
        initialize();
    }

    /**
     * This method initializes the view.
     * 
     * @param showEvents
     *            Whether events are to be shown.
     */
    public SparseMatrixPanel(boolean showEvents) {
        super();
        initialize();
        this.showEvents = showEvents;
    }

    /**
     * This method initializes this control
     */
    private void initialize() {
        this.setLayout(new BorderLayout());
        this.add(getJSplitPane(), java.awt.BorderLayout.CENTER);
        stateIcon = GuiUtils.createIcon("/jmarkov/gui/images/GreenBall.gif");
        eventIcon = GuiUtils.createIcon("/jmarkov/gui/images/RedSquare.gif");
        modelIcon = GuiUtils.createIcon("/jmarkov/gui/images/middle.gif");
    }

    @SuppressWarnings("unchecked")
    @Override
    public void updateMP() {
        MarkovProcess mp = super.getMP();
        root.removeAllChildren();
        if (mp == null) {
            root = new DefaultMutableTreeNode("No model loaded");
            treeModel.reload();
            jTree = null;
            treeScrollPane.setViewportView(getJTree());
            return;
        }
        State[] states = mp.getStates(false).toStateArray();// not generate
        Event[] events = mp.getEvents();
        // String label;
        for (State i : states) {
            DefaultMutableTreeNode iStateNode = new DefaultMutableTreeNode(i);
            root.add(iStateNode);
            if (showEvents) {
                for (Event e : events) {
                    Transitions<?> trans = mp.activeTransitions(i, e);
                    if (trans.size() != 0) {
                        DefaultMutableTreeNode eEventNode = new DefaultMutableTreeNode(
                                e);
                        iStateNode.add(eEventNode);
                        for (Transition tr : trans) {
                            State j = tr.getState();
                            DefaultMutableTreeNode jStateNode = new DefaultMutableTreeNode(
                                    j);
                            eEventNode.add(jStateNode);
                        }
                    }
                }
            } else {
                Transitions<State> reachable = mp.getRates(i);
                for (Transition j : reachable) {
                    DefaultMutableTreeNode jStateNode = new DefaultMutableTreeNode(
                            j.getState());
                    iStateNode.add(jStateNode);
                }
            }
        }
        root.setUserObject(mp.getClass().getName());
        treeModel.reload();
        treeModel.nodeChanged(root);
        jTree.setRootVisible(true);
        TreePath rootPath = new TreePath(root.getPath());
        jTree.expandPath(rootPath);
        jTree.setSelectionPath(rootPath);
        jTree.validate();
    }

    /**
     * This method initializes jTree
     * 
     * @return javax.swing.JTree
     */
    private JTree getJTree() {
        if (jTree == null) {
            root = new DefaultMutableTreeNode("No data yet ...");
            treeModel = new DefaultTreeModel(root);
            jTree = new JTree(treeModel);
            ToolTipManager.sharedInstance().registerComponent(jTree);
            jTree.setShowsRootHandles(true);
            jTree.setName("States Tree");
            jTree.setToolTipText("Rates Tree view");
            jTree.getSelectionModel().setSelectionMode(
                    TreeSelectionModel.SINGLE_TREE_SELECTION);
            jTree.addTreeSelectionListener(this);
            jTree.setCellRenderer(new MarkovTreeCellRenderer());
        }
        return jTree;
    }

    /**
     * This method initializes tableScrollPane
     * 
     * @return javax.swing.JScrollPane
     */
    private JScrollPane getInfoScrollPane() {
        if (infoScrollPane == null) {
            infoScrollPane = new JScrollPane();
            infoScrollPane.setViewportView(getTreeInfoPanel());
        }
        return infoScrollPane;
    }

    private JEditorPane getTreeInfoPanel() {
        if (treeInfoPanel == null) {
            treeInfoPanel = new JEditorPane();
            treeInfoPanel.setContentType("text/html");
        }
        return treeInfoPanel;
    }

    /**
     * This method initializes jSplitPane
     * 
     * @return javax.swing.JSplitPane
     */
    private JSplitPane getJSplitPane() {
        if (jSplitPane == null) {
            jSplitPane = new JSplitPane();
            jSplitPane.setRightComponent(getInfoScrollPane());
            jSplitPane.setLeftComponent(getTreeScrollPane());
            jSplitPane.setDividerLocation(400);
            jSplitPane.setCursor(new java.awt.Cursor(
                    java.awt.Cursor.DEFAULT_CURSOR));
        }
        return jSplitPane;
    }

    /**
     * This method initializes treeScrollPane
     * 
     * @return javax.swing.JScrollPane
     */
    private JScrollPane getTreeScrollPane() {
        if (treeScrollPane == null) {
            treeScrollPane = new JScrollPane();
            treeScrollPane.setViewportView(getJTree());
        }
        return treeScrollPane;
    }

    class MarkovTreeCellRenderer extends DefaultTreeCellRenderer {

        private static final long serialVersionUID = -9073174257320648115L;

        @Override
        public Component getTreeCellRendererComponent(JTree tree, Object value,
                boolean sel, boolean expanded, boolean leaf, int row,
                boolean hasFocus) {
            super.getTreeCellRendererComponent(tree, value, sel, expanded,
                    leaf, row, hasFocus);
            if (value instanceof DefaultMutableTreeNode) {
                Object ob = ((DefaultMutableTreeNode) value).getUserObject();
                if (ob instanceof State) {
                    State s = (State) ob;
                    setToolTipText(s.description());
                    if (stateIcon != null)
                        setIcon(stateIcon);
                } else if (ob instanceof Event) {
                    Event e = (Event) ob;
                    setToolTipText(e.description());
                    if (eventIcon != null)
                        setIcon(eventIcon);
                } else {
                    setToolTipText(null);
                    if (modelIcon != null)
                        setIcon(modelIcon);
                }
            }
            return this;
        }
    }

    /**
     * Called on node selection
     * 
     * @see javax.swing.event.TreeSelectionListener#valueChanged(javax.swing.event.TreeSelectionEvent)
     */
    public void valueChanged(TreeSelectionEvent te) {
        DefaultMutableTreeNode node = (DefaultMutableTreeNode) jTree
                .getLastSelectedPathComponent();
        if (node == null)
            return;
        Object ob = node.getUserObject();
        if (node.isLeaf() && ob instanceof State) {
            State j = (State) ob;
            if (showEvents) {
                DefaultMutableTreeNode eNode = (DefaultMutableTreeNode) node
                        .getParent();
                Event e = (Event) eNode.getUserObject();
                DefaultMutableTreeNode iNode = (DefaultMutableTreeNode) eNode
                        .getParent();
                State i = (State) iNode.getUserObject();
                showInformation(i, j, e);
            } else {// no events
                DefaultMutableTreeNode iNode = (DefaultMutableTreeNode) node
                        .getParent();
                State i = (State) iNode.getUserObject();
                showInformation(i, j);
            }

        } else if (ob instanceof State) {
            State i = (State) ob;
            showInformation(i);
        } else if (ob instanceof Event) { // event node chosen
            Event e = (Event) ob;
            DefaultMutableTreeNode iNode = (DefaultMutableTreeNode) node
                    .getParent();
            State i = (State) iNode.getUserObject();
            showInformation(i, e);
        } else if (mp != null) {
            showInformation(mp);
        }

    }

    /**
     * Show the information of transition from i to j under e.
     * 
     * @param i
     *            origin
     * @param j
     *            destination
     * @param e
     *            event
     */

    @SuppressWarnings("unchecked")
    private void showInformation(State i, State j, Event e) {
        String stg = "<html>";
        stg += "<font color=+2 color=red>Origin: </font><br>";
        stg += stateToHtml(i);
        stg += "<font color=+2 color=red>Destination: </font><br>";
        stg += stateToHtml(j);
        stg += "<font color=+2 color=red>Event: </font><br>";
        stg += eventToHtml(e);
        // TODO Clean this for general MP
        if (mp instanceof SimpleMarkovProcess){
            SimpleMarkovProcess smp = (SimpleMarkovProcess)mp;
            double rate = smp.rate(i,j,e);
             stg += "<font color=+2 color=red>Rate = " + rate + "</font><br>";
        }
        stg += "</html>";
        treeInfoPanel.setText(stg);
    }

    /**
     * Show the information of transition from i to j under e.
     * 
     * @param i
     *            origin
     * @param j
     *            destination
     */
    @SuppressWarnings("unchecked")
    private void showInformation(State i, State j) {
        String stg = "<html>";
        stg += "<font color=+2 color=red>Origin:</font><br>";
        stg += stateToHtml(i);
        stg += "<font color=+2 color=red>Destination:</font><br>";
        stg += stateToHtml(j);
        stg += "<font color=+2 color=red>Rate = " + mp.getRate(i, j)
                + "</font><br>";
        stg += "</html>";
        treeInfoPanel.setText(stg);
    }

    private void showInformation(State i, Event e) {
        String stg = "<html>";
        stg += stateToHtml(i) + eventToHtml(e) + "</html>";
        treeInfoPanel.setText(stg);
    }

    @SuppressWarnings("unchecked")
    private void showInformation(State i) {
        String stg = stateToHtml(i);
        if (mp.isGenerated()) {
            try {
                stg += "Steady State Probability = "
                        + mp.getSteadyState()[i.getIndex()];
            } catch (NotUnichainException e) {
                stg += "ERROR!!";
            }
        }
        treeInfoPanel.setText(stg);
    }

    @SuppressWarnings("unchecked")
    private void showInformation(MarkovProcess mp) {
        treeInfoPanel.setText("<pre>" + mp.description() + "</pre>");
    }

    private String stateToHtml(State s) {
        String stg = "";
        stg += "<font color=blue> State " + s.label() + "</font><br>";
        stg += s.description() + "<br>";
        return stg;
    }

    private String eventToHtml(Event e) {
        String stg = "";
        stg += "<font color=blue> Event " + e.toString() + "</font><br>";
        return stg;
    }

} // @jve:decl-index=0:visual-constraint="51,-79"
