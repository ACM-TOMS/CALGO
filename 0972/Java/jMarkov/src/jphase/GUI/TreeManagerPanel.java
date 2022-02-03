package jphase.GUI;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;

import javax.swing.JTree;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.MutableTreeNode;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;
import javax.swing.border.LineBorder;
import javax.swing.event.TreeModelEvent;
import javax.swing.event.TreeModelListener;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.ImageIcon;
import javax.swing.JTextArea;
import javax.swing.JTabbedPane;
import javax.swing.JComponent;

//import java.net.URL;
//import java.io.IOException;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.Color;
import java.util.ArrayList;

import jphase.PhaseVar;
import jphase.DenseContPhaseVar;
import jphase.AbstractContPhaseVar;
//import jphase.AbstractDiscPhaseVar;


//import org.jfree.data.xy.*;
import org.jfree.data.xy.DefaultTableXYDataset;
import org.jfree.data.xy.XYSeries;
//import org.jfree.chart.plot.*;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;


/**
 * This class manages the PH in the jPhase GUI, which are represented as a tree
 * @author Juan F. Perez
 *
 */
public class TreeManagerPanel extends JPanel 
                          implements TreeSelectionListener {
    
	
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	
	/**
	 * Panel for split the main Panel
	 */
	JSplitPane splitPane;
	
	
	/**
	 * panel with multiple tabs
	 */
	private JTabbedPane varPane;
	
	/**
	 * Scroll for variable Panel
	 */
	private JScrollPane varScroll;
	
	/**
	 * Parameter Panel
	 */
	private JComponent paramPane;
	
	/**
	 * Scroll for parameter Panel
	 */
	//private JScrollPane paramScroll;
	
	/**
	 * Probability Density Function Panel
	 */
	private ChartPanel pdfPane;
	
	/**
	 * Cumulative Function panel
	 */
	private ChartPanel cdfPane;
	
	
	/**
	 * Statistics panel
	 */
	private JPanel statPane;
	
    /**
     * Root node for tree
     */
	protected DefaultMutableTreeNode rootNode;
    
	/**
	 * Tree Model for Tree
	 */
    protected DefaultTreeModel treeModel;
    
    /**
     * Tree of Phase-Type variables
     */
    protected JTree tree;
    
    /**
     * Toolkit
     */
    private Toolkit toolkit = Toolkit.getDefaultToolkit();
    
    /**
     * List of the Phase variables
     */
    private ArrayList<PhaseVarInfo> variables;
		
	
  
    
	/**
	 * Debug
	 */
	private static boolean DEBUG = false;

	/**
	 * Frame width
	 */
	private int width = 500;
	
	/**
	 * Frame height
	 */
	private int height = 300;
	
	
	/**
	 * Tree width 
	 */ 
	private int widthDiv = 200;
	
	

	/**
	 * Tree manager default constructor
	 */
    public TreeManagerPanel() {
        super(new GridLayout(1,0));

        //Create the nodes.
  
        //Create a tree that allows one selection at a time.
        rootNode = new DefaultMutableTreeNode("Phase Var Sets");
        treeModel = new DefaultTreeModel(rootNode);
        treeModel.addTreeModelListener(new MyTreeModelListener());

        tree = new JTree(treeModel);
        tree.setEditable(true);
        tree.getSelectionModel().setSelectionMode
                (TreeSelectionModel.SINGLE_TREE_SELECTION);
        tree.setShowsRootHandles(true);
        createNodes();
        
        //Set the icon for leaf nodes.
        ImageIcon leafIcon = createImageIcon("/jphase/GUI/images/var.jpg");
        ImageIcon closedIcon = createImageIcon("/jphase/GUI/images/closedSet.gif");
        ImageIcon openIcon = createImageIcon("/jphase/GUI/images/openSet.gif");
        if (leafIcon != null) {
            DefaultTreeCellRenderer renderer = new DefaultTreeCellRenderer();
            renderer.setLeafIcon(leafIcon);
            renderer.setOpenIcon(openIcon);
            renderer.setClosedIcon(closedIcon);
            tree.setCellRenderer(renderer);
        } else {
            System.err.println("Leaf icon missing; using default.");
        }

        //Listen for when the selection changes.
        tree.addTreeSelectionListener(this);
        
        //Create the scroll pane and add the tree to it. 
        JScrollPane treeView = new JScrollPane(tree);

        /*Create the HTML viewing pane.*/
        varPane = new JTabbedPane();
        
        //initHelp();
        varScroll = new JScrollPane(varPane);
    
        //Add the scroll panes to a split pane.
        splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
        splitPane.setTopComponent(treeView);
        splitPane.setBottomComponent(varScroll);

        Dimension minimumSize = new Dimension(100, 50);
        varScroll.setMinimumSize(minimumSize);
        treeView.setMinimumSize(minimumSize);
        splitPane.setDividerLocation(widthDiv); //XXX: ignored in some releases
        								//of Swing. bug 4101306
        //workaround for bug 4101306:
        //treeView.setPreferredSize(new Dimension(100, 100)); 

        splitPane.setPreferredSize(new Dimension(width, height));

        //Add the split pane to this panel.
        add(splitPane);
    }

    /** Required by TreeSelectionListener interface. */
public void valueChanged(TreeSelectionEvent e) {
        DefaultMutableTreeNode node = (DefaultMutableTreeNode)
                           tree.getLastSelectedPathComponent();
        if (node == null) {
        	System.out.println("lastSelectedPathComponent es null");
        	return;
        }

        Object nodeInfo = node.getUserObject();
        if (node.isLeaf()) {
            PhaseVarInfo info = (PhaseVarInfo)nodeInfo; 
            displayVar(info);
            if (DEBUG) {
                //add debug print 
            }
        } else {
 
        }
        if (DEBUG) {
            System.out.println(nodeInfo.toString());
        }
    }

   /**
    * Display variable in the tab pane
    * @param info Variable to display
    */
    public void displayVar(PhaseVarInfo info) {
            if (info != null) {
                if(paramPane!=null)paramPane.removeAll();
                if(pdfPane!=null)pdfPane.removeAll();
                if(cdfPane!=null)cdfPane.removeAll();
                if(statPane!=null)statPane.removeAll();
                if(varPane!=null)varPane.removeAll();
            	
            	//PARAMETERS Tab           	
            	JTextArea param = new JTextArea(
                		info.toString()+"\n"+info.var.description());
                param.setEditable(false);
                param.setPreferredSize(new java.awt.Dimension(Math.max(100 + info.var.getVectorArray().length*95,height), height));
                
                ImageIcon icon = createImageIcon("/jphase/GUI/images/var.jpg");
                param.setBorder(new LineBorder(Color.BLACK));
                paramPane = new JPanel();
                
                paramPane.setPreferredSize(
                		new java.awt.Dimension(width, height));
                paramPane.add(param, BorderLayout.CENTER);
                
                varPane.addTab("Parameters", icon,
                		paramPane, "Distribution Parameters");
                
                //PDF Tab
                XYSeries dataPDF = new XYSeries("PDF", true, false);
                double mean = info.var.expectedValue();
                double sigma = info.var.stdDeviation();
        		double xMax = mean + Math.max(3.0 * sigma, 0.5*mean);
        		double cv = info.var.CV();
        		double numPoints = 20 + ( (!Double.isNaN(cv)) ? (int)(5 * cv ) :0 );
        		numPoints = (int)Math.max(400,numPoints);
        		double dx= xMax / numPoints;
        		double x = 0;
        		if(info.var.getNumPhases()==1)x = dx;
                for(int i = 0; i<numPoints; i++){
                	dataPDF.add(x, ((AbstractContPhaseVar)info.var).pdf(x) );
                	x+=dx;
                }
                DefaultTableXYDataset datasetPDF = new DefaultTableXYDataset(); 
                datasetPDF.addSeries(dataPDF);
                JFreeChart chartPDF = ChartFactory.createXYLineChart("PDF",
                		"x", "f(x)", datasetPDF, PlotOrientation.VERTICAL, 
                		true, true, true);
                pdfPane = new ChartPanel(chartPDF);
                pdfPane.setPreferredSize(new java.awt.Dimension(width, height));
                varPane.addTab("PDF", icon,
                		pdfPane, "Probability Density Function");
                
                //CDF Tab TODO
                XYSeries dataCDF = new XYSeries("CDF", true, false);
                x=0;
                for(int i = 0; i<numPoints; i++){
                	dataCDF.add(x, info.var.cdf(x));
                	x+=dx;
                }
                DefaultTableXYDataset datasetCDF = new DefaultTableXYDataset(); 
                datasetCDF.addSeries(dataCDF);

                JFreeChart chartCDF = ChartFactory.createXYLineChart("CDF",
                		"x", "F(x)", datasetCDF, PlotOrientation.VERTICAL, 
                		true, true, true);
                cdfPane = new ChartPanel(chartCDF);
                cdfPane.setPreferredSize(new java.awt.Dimension(width, height));
                varPane.addTab("CDF", icon,
                		cdfPane, "Cumulative Probability Function");
                
                //STAT Tab
                String stats = "Expected Value: " + String.format("%6.4f", info.var.expectedValue());
                stats += "\nVariance: " + String.format("%6.4f", info.var.variance());
                stats += "\nNumber of Phases: " + info.var.getNumPhases();
                stats += "\n--------------------------------------------";
                stats += "\nMoment 1: " + String.format("%6.4f", info.var.moment(1));
                stats += "\nMoment 2: " + String.format("%6.4f", info.var.moment(2));
                stats += "\nMoment 3: " + String.format("%6.4f", info.var.moment(3));
                stats += "\nMoment 4: " + String.format("%6.4f", info.var.moment(4));
                param = new JTextArea(stats);
                statPane = new JPanel();
                statPane.setPreferredSize(
                		new java.awt.Dimension(width, height));
                statPane.add(param, BorderLayout.CENTER);
                varPane.addTab("Statistics", icon,
                		statPane, "Distribution Statistics");
                
                
            } else { 
            	varPane.setBackground(Color.WHITE);
                if (DEBUG) {
                    System.out.println("Attempted to display a null URL.");
                }
            }
    }

    /**
     * Create the initial nodes in the tree
     *
     */
    private void createNodes() {
        /**
         * Initial tree Node
         */
    	DefaultMutableTreeNode category1;
        PhaseVarInfo var = null;
        category1 = addObject(null, "Set 1");
        
        variables = new ArrayList<PhaseVarInfo>();
        var = new PhaseVarInfo("Expo 5", 
        		DenseContPhaseVar.expo(5));
        addObject(category1, var);
        variables.add(var);
        
        var = new PhaseVarInfo("Erlang (5, 2)", 
        		DenseContPhaseVar.Erlang(5,2));
        addObject(category1, var);
        variables.add(var);
        
        var = new PhaseVarInfo("Expo 15", 
        		DenseContPhaseVar.expo(15));
        addObject(category1, var);
        variables.add(var);
        
        var = new PhaseVarInfo("Erlang (15, 3)", 
        		DenseContPhaseVar.Erlang(15,3));
        addObject(category1, var);
        
        variables.add(var);
    }
    
    /**
     * Add new set of nodes to the tree
     * @param name name of the new set
     */
    public void addSet(String name){
        DefaultMutableTreeNode cat;
        PhaseVarInfo info = null;
        cat = addObject(null, name);
        info = new PhaseVarInfo("New Var", DenseContPhaseVar.expo(2));
        addObject(cat, info);
    }
    
    /**
     * Returns True if a variable can be added
     * @return True if a variable can be added
     */
    public boolean canAddVar(){
    	TreePath parentPath = tree.getSelectionPath();
    	if(parentPath == null || 
        		parentPath.getLastPathComponent() == rootNode){
        	toolkit.beep();
        	return false;
        }else return true;
    }
    
    /**
     * Add a new variable to the tree in the selected set
     * @param name variable name
     * @param var PH Variable
     */
    public void addVar(String name, PhaseVar var){
        PhaseVarInfo child = new PhaseVarInfo(name, var);
        DefaultMutableTreeNode parentNode = null;
        TreePath parentPath = tree.getSelectionPath();
        parentNode = (DefaultMutableTreeNode)(parentPath.getLastPathComponent());
        if(parentNode.getParent()!=rootNode)
            	parentNode = (DefaultMutableTreeNode)parentNode.getParent();
        
        addObject(parentNode, child, true);
        variables.add(child);
    }
    
    
    /** Remove all nodes except the root node. */
    public void clear() {
        rootNode.removeAllChildren();
        treeModel.reload();
    }
    
    /** Remove the currently selected node. */
    public void removeCurrentNode() {
        TreePath currentSelection = tree.getSelectionPath();
        if (currentSelection != null) {
            DefaultMutableTreeNode currentNode = (DefaultMutableTreeNode)
                         (currentSelection.getLastPathComponent());
            MutableTreeNode parent = (MutableTreeNode)(currentNode.getParent());
            if (parent != null) {
                treeModel.removeNodeFromParent(currentNode);
                return;
            }
        } 

        // Either there was no selection, or the root was selected.
        toolkit.beep();
    }

    /** 
     * Add child to the currently selected node. 
     * @param child
     *//*
    public DefaultMutableTreeNode addObject(Object child) {
        DefaultMutableTreeNode parentNode = null;
        TreePath parentPath = tree.getSelectionPath();

        if (parentPath == null) {
            parentNode = rootNode;
        } else {
            parentNode = (DefaultMutableTreeNode)
                         (parentPath.getLastPathComponent());
        }

        return addObject(parentNode, child, true);
    }*/
    
    /**
     * Add an object to the tree
     * @param parent Parent that receive the new child
     * @param child Object to be added
     * @return Modified Tree
     */
    public DefaultMutableTreeNode addObject(
    		DefaultMutableTreeNode parent, Object child) {
        return addObject(parent, child, false);
    }
    
    /**
     * Add an object to the tree in the selected parent
     * @param parent Parent that receive the new child
     * @param child Object to be added
     * @param shouldBeVisible True if the new child should be visible
     * @return Modified Tree
     */
    public DefaultMutableTreeNode addObject(DefaultMutableTreeNode parent,
                                            Object child, boolean shouldBeVisible) {
        DefaultMutableTreeNode childNode = 
                new DefaultMutableTreeNode(child);

        if (parent == null) {
            parent = rootNode;
        }
        treeModel.insertNodeInto(childNode, parent, parent.getChildCount());

        if (shouldBeVisible) {
            tree.scrollPathToVisible(new TreePath(childNode.getPath()));
        }
        return childNode;
    }
    
    /** 
     * Returns an ImageIcon, or null if the path was invalid. 
     * @param path Path for the image to be selected
     * @return ImageIcon
     */
    protected static ImageIcon createImageIcon(String path) {
        java.net.URL imgURL = TreeManagerPanel.class.getResource(path);
        if (imgURL != null) {
            return new ImageIcon(imgURL);
        } else {
            System.err.println("Couldn't find file: " + path);
            return null;
        }
    }
    
    /**
     * TreeModel Listener for PH Variable tree
     * @author Juan F. Perez
     *
     */
    class MyTreeModelListener implements TreeModelListener {

		public void treeNodesChanged(TreeModelEvent e) {
            DefaultMutableTreeNode node;
            node = (DefaultMutableTreeNode)
                     (e.getTreePath().getLastPathComponent());
            try {
                int index = e.getChildIndices()[0];
                node = (DefaultMutableTreeNode)
                       (node.getChildAt(index));
            } catch (NullPointerException exc) {}

            System.out.println("The user has finished editing the node.");
            System.out.println("New value: " + node.getUserObject());
        }
        
		public void treeNodesInserted(TreeModelEvent e) {
        }
        
		public void treeNodesRemoved(TreeModelEvent e) {
        }
        
		public void treeStructureChanged(TreeModelEvent e) {
        }
    }
	/**
	 * @return the variables
	 */
	public ArrayList<PhaseVarInfo> getVariables() {
		return variables;
	}   
}
