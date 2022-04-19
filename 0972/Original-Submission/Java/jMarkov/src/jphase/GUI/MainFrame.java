package jphase.GUI;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.io.File;
import java.util.ArrayList;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JDialog;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JToolBar;
import javax.swing.UIManager;
import javax.swing.JOptionPane;

/**
 * <p>Title: MainFrame</p>
 * <p>Description: Graphic User Interface for JPhase</p>
 * <p>Copyright: Copyright (c) 2006 - 2014</p>
 * @author Juan F. Pérez
 * @version 1.0
 */
public class MainFrame extends JFrame {
	/**
	 * Class Version
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * Instance
	 */
	private static MainFrame instance;
	
	/**
	 * Text area to provide information to the user
	 */
	private JTextArea log;
	
	/**
	 * Menu bar
	 */
	private JMenuBar menuBar;
	
	/**
	 * File Menu
	 */
	private JMenu fileMenu;
	
	/**
	 * Item New in File menu 
	 */
	private JMenuItem newItem;

	/**
	 * Item Open in File menu 
	 */
	private JMenuItem openItem;

	/**
	 * Item Save in File menu 
	 */
	private JMenuItem saveItem;

	/**
	 * Item Close in File menu 
	 */
	private JMenuItem closeItem;

	/**
	 * Item Exit in menu File
	 */
	private JMenuItem exitItem;
	
	/**
	 * PhaseVar Menu	
	 */
	private JMenu phVarMenu;
	
	/**
	 * Item newSet in PhaseVar menu
	 */
	private JMenuItem newSetItem;
	
	/**
	 * Item newVar in PhaseVar menu
	 */
	private JMenuItem newVarItem;
	
	/**
	 * Item newVar in PhaseVar menu
	 */
	private JMenuItem newVarGItem;
	
	/**
	 * Item newVar in PhaseVar menu
	 */
	private JMenuItem newQueue;

	/**
	 * Item PDF in PhaseVar menu
	 */
	//private JMenuItem pdfItem;
	
	/**
	 * Item CDF in PhaseVar File
	 */
	//private JMenuItem cdfItem;
	
	
	/**
	 * Fit Menu
	 */
	//private JMenu fitMenu;
	
	/**
	 * Item LoadData in PhaseVar menu
	 */
	//private JMenuItem loadItem;

	/**
	 * Item DoFit in PhaseVar menu
	 */
	//private JMenuItem doFitItem;
	
	/**
	 * Help Menu
	 */
	private JMenu helpMenu;
	
	/**
	 * Item About in Help menu
	 */
	private JMenuItem aboutItem;
	
	/**
	 * Tool Bar
	 */
	private JToolBar toolBar;
	
	
	/**
	 * Main Pane 
	 */
	private JPanel mainPane;
	
	/**
	 * 
	 */
	private TreeManagerPanel treeContentPane;
	
	/**
	 * 
	 */
	private int width = 700;
	
	/**
	 * 
	 */
	private int height = 400;
	
	
	/**
	 * JFileChooser to load data in files
	 */
	//private JFileChooser fileExplorer;
	
	/**
	 * @param args 
	 */
	 public static void main(String[] args) {
	    try {
	      UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
	    }
	    catch(Exception e) {
	      e.printStackTrace();
	    }
	    new MainFrame();
	}
	
	
	
	/**
	 * Main frame instance
	 * @return instance 
	 */
	public static MainFrame getInstance() {
	    if (instance == null) {
	      instance = new MainFrame();
	    }
	    return instance;
	  }
	
	/**
	 * Main frame constructor 
	 */
	public MainFrame(){
		this.setTitle("JPhase - Java Framework for Phase Type Distributions");
		this.setIconImage(new ImageIcon("images/var.jpg").getImage());
		this.setDefaultCloseOperation(EXIT_ON_CLOSE);
		this.setFocusable(true);
		this.setResizable(true);
		
		//Menu bar 
		menuBar = new JMenuBar();
		menuBar.setVisible(true);
		menuBar.setPreferredSize(new java.awt.Dimension(width, 25));
		setJMenuBar(menuBar);
		
		fileMenu = new JMenu("File");
		fileMenu.setMnemonic(KeyEvent.VK_F);
		menuBar.add(fileMenu);
		
		newItem = new JMenuItem("New");
		newItem.addActionListener(
				new MainFrame_newItem_actionAdapter(this));
		fileMenu.add(newItem);
		
		openItem = new JMenuItem("Open");
		fileMenu.add(openItem);
		
		saveItem = new JMenuItem("Save");
		fileMenu.add(saveItem);
		
		closeItem = new JMenuItem("Close");
		closeItem.addActionListener(
				new MainFrame_closeItem_actionAdapter(this));
		fileMenu.add(closeItem);
		
		exitItem = new JMenuItem("Exit");
		exitItem.addActionListener(
				new MainFrame_exitItem_actionAdapter(this));
		fileMenu.add(exitItem);
	
		
		phVarMenu = new JMenu("Phase Var");
		phVarMenu.setMnemonic(KeyEvent.VK_P);
		menuBar.add(phVarMenu);
				
		newSetItem = new JMenuItem("New Set");
		newSetItem.addActionListener(
				new MainFrame_newSetItem_actionAdapter(this));
		phVarMenu.add(newSetItem);
		
		
		newVarItem = new JMenuItem("New Phase-type Variable");
		newVarItem.addActionListener(
				new MainFrame_newVarItem_actionAdapter(this));
		phVarMenu.add(newVarItem);		
		
		newVarGItem = new JMenuItem("New General Variable");
		newVarGItem.addActionListener(
				new MainFrame_newVarGItem_actionAdapter(this));
		phVarMenu.add(newVarGItem);		
		
		newQueue = new JMenuItem("Generate a PH/PH/1 Queue");
		newQueue.addActionListener(
				new MainFrame_newQueue_actionAdapter(this));
		phVarMenu.add(newQueue);
		
		newQueue = new JMenuItem("GOF for a Phase-Type Variable");
		newQueue.addActionListener(
				new MainFrame_GROPhase_actionAdapter(this));
		phVarMenu.add(newQueue);
		
		/*pdfItem = new JMenuItem("PDF");
		phVarMenu.add(pdfItem);*/
		
		/*cdfItem = new JMenuItem("CDF");
		phVarMenu.add(cdfItem);*/
		
		/*fitMenu = new JMenu("Fit");
		fitMenu.setMnemonic(KeyEvent.VK_I);
		menuBar.add(fitMenu);*/
		
		/*doFitItem = new JMenuItem("Do Fit");
		fitMenu.add(doFitItem);*/
		
		helpMenu = new JMenu("Help");
		helpMenu.setMnemonic(KeyEvent.VK_H);
		menuBar.add(helpMenu);
		
		aboutItem = new JMenuItem("About JPhase");
		aboutItem.addActionListener(
				new MainFrame_aboutItem_actionAdapter(this));
		helpMenu.add(aboutItem);
		
		
		toolBar = new JToolBar("ToolBar");
		JButton newButton = new JButton(new ImageIcon("images/new.jpg"));
		newButton.setToolTipText("Cargar archivo caudal entrante");
		newButton.addActionListener(
				new MainFrame_newItem_actionAdapter(this));
		toolBar.add(newButton);
		
		JButton openButton = new JButton(new ImageIcon("images/open.jpg"));
		openButton.setToolTipText("Create a new distribution");
		toolBar.add(openButton);
		
		JButton saveButton = new JButton(new ImageIcon("images/save.jpg"));
		saveButton.setToolTipText("Save distribution");
		toolBar.add(saveButton);
		
		add(toolBar, BorderLayout.PAGE_START);
		
		
		mainPane = new JPanel();
		mainPane.setLayout(null);
		mainPane.setPreferredSize(new java.awt.Dimension(width, height));
		
		/*
		JLabel fechaInicio= new JLabel("Fecha de inicio: ");
		fechaInicio.setForeground(Color.DARK_GRAY);
		fechaInicio.setBounds(new Rectangle(10, 10, 110, 20));
		fechaInicio.setHorizontalAlignment(SwingConstants.RIGHT);
		principalPanel.add(fechaInicio,null);
		*/
		
		/*
		JLabel fechaFin= new JLabel("Fecha de Finalización: ");
		fechaFin.setForeground(Color.DARK_GRAY);
		fechaFin.setBounds(new Rectangle(10, 45, 110, 20));
		fechaFin.setHorizontalAlignment(SwingConstants.RIGHT);
		principalPanel.add(fechaFin,null);
		*/
		
		mainPane.setBackground(Color.LIGHT_GRAY);
		this.getContentPane().add(mainPane, BorderLayout.CENTER);
		
			
		//Text area for logging 
		log = new JTextArea(2, 8);
		log.setFont(new java.awt.Font("@Arial Unicode MS", 0, 11));
		log.setEditable(false);
		JScrollPane logScrollPane = new JScrollPane(log);
		this.getContentPane().add(logScrollPane, BorderLayout.SOUTH);
		this.pack();
		
		//Screen size
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		Dimension frameSize = this.getSize();
		if (frameSize.height > screenSize.height) frameSize.height = screenSize.height;
		if (frameSize.width > screenSize.width) frameSize.width = screenSize.width;
		this.setLocation((screenSize.width - frameSize.width) / 2, (screenSize.height - frameSize.height) / 2);
		this.setVisible(true);
	}
	
	/**
	 * Add a message to the log panel
	 * @param str String to be added to the log panel
	 */
	public void addLog(String str){
		if(str!=null)
			this.log.append(str);
	}
	
	
	  /**
	   * Creates a new workspace
	   * @param e Event.
	   */
	  void newItem_actionPerformed(ActionEvent e) {
		if(treeContentPane!=null){
			  this.getContentPane().remove(treeContentPane);
			  treeContentPane.tree.removeAll();
			  treeContentPane.treeModel.setRoot(null);
			  treeContentPane.removeAll();
			  treeContentPane = null;
		}
	    treeContentPane = new TreeManagerPanel();
	    treeContentPane.setOpaque(true); //content panes must be opaque
	    treeContentPane.setPreferredSize(new java.awt.Dimension(width, height));
	    if(mainPane!=null){this.remove(mainPane);
	    mainPane=null;}
        this.getContentPane().add(treeContentPane, BorderLayout.CENTER);
	    this.pack();
	  }
	
	
	  /**
	   * Creates a new workspace
	   * @param e Event.
	   */
	  void closeItem_actionPerformed(ActionEvent e) {
		  int res = JOptionPane.showConfirmDialog(null, 
				  "Are you sure you want to close the present JPhase calculator?",
				  "JPhase: Close?",
				  JOptionPane.YES_NO_OPTION);
		  if(res == JOptionPane.YES_OPTION){
			  if(treeContentPane!=null){
				  this.getContentPane().remove(treeContentPane);
				  treeContentPane.tree.removeAll();
				  treeContentPane.treeModel.setRoot(null);
				  treeContentPane.removeAll();
				  treeContentPane = null;
			  }
			  mainPane = new JPanel();
			  mainPane.setLayout(null);
			  mainPane.setPreferredSize(new java.awt.Dimension(width, height));
			  mainPane.setBackground(Color.LIGHT_GRAY);
			  this.getContentPane().add(mainPane, BorderLayout.CENTER);
			  this.getContentPane().add(mainPane);
			  this.getContentPane().setBackground(Color.LIGHT_GRAY);
			  this.pack();
		  }	  
	  }
	  
	  /**
	   * Creates a new workspace
	   * @param e Event.
	   */
	  void exitItem_actionPerformed(ActionEvent e) {
		  int res = JOptionPane.showConfirmDialog(null, 
				  "Are you sure you want to exit JPhase?",
				  "JPhase: Exit?",
				  JOptionPane.YES_NO_OPTION);
		  if(res == JOptionPane.YES_OPTION)System.exit(0);
	  }
	  
	  /**
	   * Creates a new workspace
	   * @param e Event.
	   */
	  void newSetItem_actionPerformed(ActionEvent e) {
		  if(treeContentPane==null){
			  JOptionPane.showMessageDialog(null,  
					  "There must be an active file to create a Set", 
					  "JPhase Alert", 
					  JOptionPane.INFORMATION_MESSAGE);
				
		  }else{
			  InputFrame newSetFrame = new InputFrame(
					  "JPhase: New Set Creation",
					  "Enter the name of the new Set");
			  newSetFrame.setVisible(true);
			  newSetFrame.setFocusable(true);
			  if(newSetFrame.getRes()){
				  treeContentPane.addSet(newSetFrame.getValue());
				  this.pack();  
			  }
		  }
	  }
	  
	  
	  /**
	   * Actions performed when newVarItem is selected
	   * @param e
	   */
	  void newVarItem_actionPerformed(ActionEvent e) {
		  if(treeContentPane==null){
			  JOptionPane.showMessageDialog(null,  
					  "There must be an active file to create a Variable", 
					  "JPhase Alert", 
					  JOptionPane.INFORMATION_MESSAGE);
		  }else{ 
			   boolean res = treeContentPane.canAddVar();
			   if(!res){log.append("New variable cannot be created " +
			   		"without specifying a Set.\n");
			   }else{
				   NewVarFrame myFrame = new NewVarFrame();
				   myFrame.setVisible(true);
				   treeContentPane.addVar(myFrame.varName, myFrame.var);
				   this.pack();
			   }
		  }
	  }
	  
	  /**
	   * Actions performed when newVarGItem is selected
	   * @param e
	   */
	  void newVarGItem_actionPerformed(ActionEvent e) {
		  JFrame interfazVariables = new GenerateVarFrame();
		  interfazVariables.setVisible(true);
	  }
	  
	  /**
	   * Actions performed when newVarGItem is selected
	   * @param e
	   */
	  void newQueue_actionPerformed(ActionEvent e) {
		  JDialog dialog = new NewQueueDialog(this);
		  dialog.setVisible(true);
	  }
	  
	  /**
	   * Actions performed when newVarGItem is selected
	   * @param e
	   */
	  void GOFPhase_actionPerformed(ActionEvent e) {
		  JDialog dialog = new PhaseGOF(this);
		  dialog.setVisible(true);
	  }
	  
	  /**
	   * 
	   * @param e
	   */
	  void aboutItem_actionPerformed(ActionEvent e) {
			String s = "JPhase was developed by \n " +
			"  Juan F. Pérez \n" +
			"  Germán Riaño Mendoza \n\n" +
			"  Andrés Sarmiento \n\n" +
			"Universidad de los Andes\n" +
			"Bogotá, Colombia\n" +
			"Version 1.0 (2006-2014)";
			
			JOptionPane.showMessageDialog(null,  s, 
					"JPhase About", JOptionPane.INFORMATION_MESSAGE);
	  }
	  
	  
	
	/**
	 * Load the input data file
	 */ /*
	private void addItemLoadData() {
		loadItem = new JMenuItem("Load Data");
		fileExplorer = new JFileChooser();
		fileExplorer.setFileFilter(new ExcelFilter());
		fileExplorer.setFileFilter(new TextFilter());
		class CargarListener implements ActionListener {
			public void actionPerformed(ActionEvent evt) {
				log.selectAll();
				log.replaceSelection("");
				int returnVal = fileExplorer.showDialog(principalPane, "Cargar");
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					File file = fileExplorer.getSelectedFile();
					
					//Verificacion de extension .txt
					String extension=getSuffix(file);
										
					if(extension.equalsIgnoreCase("txt") || extension.equalsIgnoreCase("xls")){
						//Load data file 
						Fachada.getInstance().loadData(file,extension);
						log.append("Data file "+file.getName()+ " loaded\n");
					}
					else log.append("Error on data loading. The file must be an Excel or Text file\n");
				} 
				else log.append("Loading operation canceled by user.\n");
				log.setCaretPosition(log.getDocument().getLength());
			}
		}
		ActionListener listener = new CargarListener();
		loadItem.addActionListener(listener);
		fitMenu.add(loadItem);
	}*/
	
	
	/**
	  * Returns file extension  
	  * @return String File Extension
	  * @param f File  
	  */ 
	private String getSuffix(File f) {
		String s = f.getPath(), suffix = null;
		int i = s.lastIndexOf('.');
		if (i>0 && i<s.length()-1) suffix = s.substring(i+1).toLowerCase();
		return suffix;
	}
	 
	
	
	/**
	 * Filter for excel files 
	 */ /*
	private class ExcelFilter extends javax.swing.filechooser.FileFilter {
		@Override
		public boolean accept(File f) {
			boolean accept = f.isDirectory();
			if (!accept) {
				String suffix = getSuffix(f);
				if (suffix != null)accept = suffix.equals("xls");
				}
			return accept;
			}
		
		@Override
		public String getDescription() {
			return "Archivos de excel(*.xls)";
			}
		}
	 */
	/**
	 * Filter for text files 
	 */
	class TextFilter extends javax.swing.filechooser.FileFilter {
		
		@Override
		public boolean accept(File f) {
			boolean accept = f.isDirectory();
			if (!accept) {
				String suffix = getSuffix(f);
				if (suffix != null)accept = suffix.equals("txt");
				}
			return accept;	
			}
		
		@Override
		public String getDescription() {
			return "Text files(*.txt)";
			}
		}
	
	/**
	 * Main frame action listener
	 */
	class MainFrame_newItem_actionAdapter
	    implements java.awt.event.ActionListener {
		
		/**
		 * 
		 */
		MainFrame adaptee;
	
		/**
		 * 
		 * @param adaptee
		 */
		MainFrame_newItem_actionAdapter(MainFrame adaptee) {
			this.adaptee = adaptee;
		}
	
		public void actionPerformed(ActionEvent e) {
			adaptee.newItem_actionPerformed(e);
		}
	}
	
	
	
	/**
	 * Close item action listener
	 * @author Juan F. Perez
	 *
	 */
	class MainFrame_closeItem_actionAdapter
    		implements java.awt.event.ActionListener {
		/**
		 * 
		 */
		MainFrame adaptee;
	
		/**
		 * 
		 * @param adaptee
		 */
		MainFrame_closeItem_actionAdapter(MainFrame adaptee) {
			this.adaptee = adaptee;
		}
	
		public void actionPerformed(ActionEvent e) {
			adaptee.closeItem_actionPerformed(e);
		}
	}
	
	
	/**
	 * Exit item action listener
	 * @author Juan F. Perez
	 *
	 */
	class MainFrame_exitItem_actionAdapter
			implements java.awt.event.ActionListener {
		
		
		/**
		 * 
		 */
		MainFrame adaptee;
		
		/**
		 * 
		 * @param adaptee
		 */
		MainFrame_exitItem_actionAdapter(MainFrame adaptee) {
			this.adaptee = adaptee;
		}

		public void actionPerformed(ActionEvent e) {
			adaptee.exitItem_actionPerformed(e);
		}
	}
	
	
	
	/**
	 * New item action listener
	 * @author Juan F. Perez
	 *
	 */
	class MainFrame_newSetItem_actionAdapter
	    implements java.awt.event.ActionListener {
		
		/**
		 * 
		 */
		MainFrame adaptee;
		
		/**
		 * 
		 * @param adaptee
		 */
		MainFrame_newSetItem_actionAdapter(MainFrame adaptee) {
			this.adaptee = adaptee;
		}
	
		public void actionPerformed(ActionEvent e) {
			adaptee.newSetItem_actionPerformed(e);
		}
	}
	
	/**
	 * New var action listener
	 * @author Juan F. Perez
	 *
	 */
	class MainFrame_newVarGItem_actionAdapter
	    implements java.awt.event.ActionListener {
		
		/**
		 * 
		 */
		MainFrame adaptee;
	
		/**
		 * Action Adapter
		 * @param adaptee
		 */
		MainFrame_newVarGItem_actionAdapter(MainFrame adaptee) {
			this.adaptee = adaptee;
		}
	
		public void actionPerformed(ActionEvent e) {
			adaptee.newVarGItem_actionPerformed(e);
		}
	}	
	
	/**
	 * New queue action listener
	 * @author Juan F. Perez
	 *
	 */
	class MainFrame_newQueue_actionAdapter
	    implements java.awt.event.ActionListener {
		
		/**
		 * 
		 */
		MainFrame adaptee;
	
		/**
		 * Action Adapter
		 * @param adaptee
		 */
		MainFrame_newQueue_actionAdapter(MainFrame adaptee) {
			this.adaptee = adaptee;
		}
	
		public void actionPerformed(ActionEvent e) {
			adaptee.newQueue_actionPerformed(e);
		}
	}
	
	/**
	 * GRO phase item action listener
	 * @author Juan F. Perez
	 *
	 */
	class MainFrame_GROPhase_actionAdapter
	    implements java.awt.event.ActionListener {
		
		/**
		 * 
		 */
		MainFrame adaptee;
	
		/**
		 * Action Adapter
		 * @param adaptee
		 */
		MainFrame_GROPhase_actionAdapter(MainFrame adaptee) {
			this.adaptee = adaptee;
		}
	
		public void actionPerformed(ActionEvent e) {
			adaptee.GOFPhase_actionPerformed(e);
		}
	}	
	
	/**
	 * New var item action listener
	 * @author Juan F. Perez
	 *
	 */
	class MainFrame_newVarItem_actionAdapter
	    implements java.awt.event.ActionListener {
		
		/**
		 * 
		 */
		MainFrame adaptee;
	
		/**
		 * Action Adapter
		 * @param adaptee
		 */
		MainFrame_newVarItem_actionAdapter(MainFrame adaptee) {
			this.adaptee = adaptee;
		}
	
		public void actionPerformed(ActionEvent e) {
			adaptee.newVarItem_actionPerformed(e);
		}
	}
	
	/**
	 * About item actions listener
	 * @author Juan F. Perez
	 *
	 */
	class MainFrame_aboutItem_actionAdapter
	    implements java.awt.event.ActionListener {
		
		/**
		 * 
		 */
		MainFrame adaptee;
		
		/**
		 * 
		 * @param adaptee
		 */
		MainFrame_aboutItem_actionAdapter(MainFrame adaptee) {
			this.adaptee = adaptee;
		}
	
		public void actionPerformed(ActionEvent e) {
			adaptee.aboutItem_actionPerformed(e);
		}
	}

	public ArrayList<PhaseVarInfo> getVariables() {
		  if(treeContentPane==null){
			  JOptionPane.showMessageDialog(null,  
					  "There must be an active file to create a Variable", 
					  "JPhase Alert", 
					  JOptionPane.INFORMATION_MESSAGE);
		  }else{ 
			   boolean res = treeContentPane.canAddVar();
			   if(!res){log.append("You must specify a Set before.\n");
			   }else{
				   return treeContentPane.getVariables();
			   }
		  }
		  return null;
	}
}
