package jphase.GUI;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.StringTokenizer;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.filechooser.FileFilter; 
import jphase.fit.MomentsACPH2Fit;

/**
 * <p>Title: FitFrame</p>
 * <p>Description: Graphic User Interface for JPhase</p>
 * <p>Copyright: Copyright (c) 2006-2014</p>
 * @author Juan F. Pérez
 * @version 1.0
 */
public class FitFrame extends JDialog{
	/**
	 * Class Version
	 */
	private static final long serialVersionUID = 1L;

	
	/**
	 * Main Panel
	 */
	private JPanel principalPanel;
	
	/**
	 * Alert label 
	 */
	JLabel alertLabel;
	
	
	/**
	 * File Chooser label 
	 */
	JLabel fileLabel;
	
	
	/**
	 * Button for File Chooser
	 */
	JButton fileButton;
	
	/**
	 * File Chooser to load data
	 */
	JFileChooser fileChooser;
	
	
	/**
	 * Parameter names
	 */
	JLabel[] paramNames;
	
	/**
	 * Parameter Values
	 */
	JTextField[] paramValues;
	
	/**
	 * Accept Button
	 */
	JButton yesButton;
	
	/**
	 * Cancel Button
	 */
	JButton noButton;
	
	/**
	 * True if all the parameters are completely specified
	 */
	private boolean res = false;
		
	/**
	 * Data loaded from file
	 */
	private double[] data;
	
	/**
	 * The boolean is activated if the frame does not requieres aditional parameters
	 */
	private boolean needPar;
	
	/**
	 * Resulting Phase variable
	 */
	//private PhaseVar var = null;
	
	/**
	 * Frame width
	 */
	private int width = 360;
	
	/**
	 * Frame height
	 */
	//private int height = 220;
	private int height = 150;
	
	/**
	 * Minimum allowed frame width 
	 */
	//private int minWidth = 360;
	
	/**
	 * Minimum allowed frame height
	 */
	//private int minHeight = 220;
	
	
	
	/**
	 * Construcción del frame principal
	 * @param fitType Variable definition in String version
	 */
	public FitFrame(String fitType){
			
			
		this.setTitle("JPhase - New Fit Parameters");
		this.setResizable(false);
		this.setModal(true);
		
		principalPanel = new JPanel();
		principalPanel.setLayout(null);
		principalPanel.setPreferredSize(new java.awt.Dimension(width, height));
		principalPanel.setOpaque(true);
		
		
		alertLabel= new JLabel("Parameters of the "+fitType+" Fit");
		//alertLabel.setForeground(Color..DARK_GRAY);
		alertLabel.setBounds(new Rectangle(width/2-110, 10, 260, 20));
		alertLabel.setFont(new java.awt.Font("@Arial Unicode MS", 1, 12));
		alertLabel.setOpaque(true);
		//alertLabel.setBackground(new java.awt.Color(211 , 210, 192));
		//alertLabel.setBorder(new LineBorder(Color.BLACK));
		principalPanel.add(alertLabel, null);
		

		fileLabel = new JLabel("File Load");
		fileLabel.setBounds(new Rectangle(width/2-110, 40, 100, 20));
		fileLabel.setFont(new java.awt.Font("@Arial Unicode MS", 1, 12));
		principalPanel.add(fileLabel, null);
		
		
		fileButton = new JButton("Browse");
		fileButton.setBounds(new Rectangle(width/2, 40, 100, 25));
		fileButton.setFont(new java.awt.Font("@Arial Unicode MS", 1, 12));
		fileButton.addActionListener(
				new ParamFrame_fileButton_actionAdapter(this));
		principalPanel.add(fileButton, null);
		
		/*String[] fits = {"Moments ACPH2", "Moments EC Complete", 
		"Moments EC Positive",	"Moments ACPH", "EMHyperExpoFit", 
		"EMHyperErlangFit", "EMPhaseFit"}; 
		 */
		
		if(fitType.equals("Moments ACPH2")){
			paramNames = new JLabel[1] ;
			paramNames[0] = new JLabel("Precision");
			paramNames[0].setFont(new java.awt.Font("@Arial Unicode MS", 1, 12));
			
			paramValues = new JTextField[1] ;
			paramValues[0] = new JTextField("" + MomentsACPH2Fit.precision);
			
			
			for(int i = 0; i < paramNames.length; i++){
				paramNames[i].setBounds(new Rectangle(width/2 - 100, 40*(i+2), 100, 20));
				principalPanel.add(paramNames[i],null );
				
				paramValues[i].setBackground(Color.WHITE);
				paramValues[i].setBounds(new Rectangle(width/2, 40*(i+2), 100, 20));
				principalPanel.add(paramValues[i],null );
			}
			
					
		}else if(fitType.equals("Moments EC Complete")){
			needPar = true;
			paramNames = new JLabel[1] ;
			paramNames[0] = new JLabel("Precision");
			paramNames[0].setFont(new java.awt.Font("@Arial Unicode MS", 1, 12));
			
			paramValues = new JTextField[1] ;
			paramValues[0] = new JTextField("" + MomentsACPH2Fit.precision);
			
			
			for(int i = 0; i < paramNames.length; i++){
				paramNames[i].setBounds(new Rectangle(width/2 - 100, 40*(i+2), 100, 20));
				principalPanel.add(paramNames[i],null );
				
				paramValues[i].setBackground(Color.WHITE);
				paramValues[i].setBounds(new Rectangle(width/2, 40*(i+2), 100, 20));
				principalPanel.add(paramValues[i],null );
			}
			
			
		}else if(fitType.equals("Moments EC Positive")){
			paramNames = new JLabel[1] ;
			paramNames[0] = new JLabel("Precision");
			paramNames[0].setFont(new java.awt.Font("@Arial Unicode MS", 1, 12));
			
			paramValues = new JTextField[1] ;
			paramValues[0] = new JTextField("" + MomentsACPH2Fit.precision);
			
			
			for(int i = 0; i < paramNames.length; i++){
				paramNames[i].setBounds(new Rectangle(width/2 - 100, 40*(i+2), 100, 20));
				principalPanel.add(paramNames[i],null );
				
				paramValues[i].setBackground(Color.WHITE);
				paramValues[i].setBounds(new Rectangle(width/2, 40*(i+2), 100, 20));
				principalPanel.add(paramValues[i],null );
			}
		}else if(fitType.equals("Moments ACPH")){
			paramNames = new JLabel[1] ;
			paramNames[0] = new JLabel("Precision");
			paramNames[0].setFont(new java.awt.Font("@Arial Unicode MS", 1, 12));
			
			paramValues = new JTextField[1] ;
			paramValues[0] = new JTextField("" + MomentsACPH2Fit.precision);
			
			
			for(int i = 0; i < paramNames.length; i++){
				paramNames[i].setBounds(new Rectangle(width/2 - 100, 40*(i+2), 100, 20));
				principalPanel.add(paramNames[i],null );
				
				paramValues[i].setBackground(Color.WHITE);
				paramValues[i].setBounds(new Rectangle(width/2, 40*(i+2), 100, 20));
				principalPanel.add(paramValues[i],null );
			}
		}else if(fitType.equals("EMHyperExpoFit")){
			needPar = false;			
		}else if(fitType.equals("EMHyperErlangFit")){
			needPar = false;
		}else if(fitType.equals("EMPhaseFit")){
			needPar = false;
		}else{
			System.out.println("Non-known fitting algorithm");		
		}
		
		
		
		
		/*
		for(int i = 0; i < paramNames.length; i++){
			paramNames[i].setBounds(new Rectangle(width/2 - 100, 40*(i+1), 80, 20));
			principalPanel.add(paramNames[i],null );
			
			paramValues[i].setBackground(Color.WHITE);
			paramValues[i].setBounds(new Rectangle(width/2 + 20, 40*(i+1), 80, 20));
			principalPanel.add(paramValues[i],null );
		}*/
		
		
		/*
		String[] vars = {"Expo","Erlang", "HyperEponential",
				"Coxian", "HyperErlang", "General Phase"}; 
		*/
		
	
		
		
		
		yesButton = new JButton("Enter");
		yesButton.setBounds(new Rectangle(width/2 - 100, height - 40, 80, 25));
		yesButton.addActionListener(
				new ParamFrame_yesButton_actionAdapter(this));
		yesButton.setOpaque(true);
		principalPanel.add(yesButton, null);
		
		noButton = new JButton("Cancel");
		noButton.setBounds(new Rectangle(width/2 + 20, height - 40, 80, 25));
		noButton.addActionListener(
				new ParamFrame_noButton_actionAdapter(this));
		principalPanel.add(noButton, null);
		
		//principalPanel.setBackground(Color.LIGHT_GRAY);
		this.getContentPane().add(principalPanel, BorderLayout.CENTER);
		this.centrarFrame();
		pack();
		
		//Dimension en pantalla
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		Dimension frameSize = this.getSize();
		if (frameSize.height > screenSize.height) frameSize.height = screenSize.height;
		if (frameSize.width > screenSize.width) frameSize.width = screenSize.width;
		this.setLocation((screenSize.width - frameSize.width) / 2, (screenSize.height - frameSize.height) / 2);
		this.setVisible(true);

		
	}
	
	/**
	 * 
	 * @return True if all parameter are fully specified
	 */
	public boolean getRes(){
		return this.res;
	}
	
	/**
	 * 
	 * @return data loaded from file
	 */
	public double[] getData(){
		return this.data;
	}
	
	/*
	public String getValue(){
		return this.textInput.getText();
	}*/
	
	/**
	 * Action performed when fileButton is clicked
	 * @param e Event
	 */
	void fileButton_actionPerformed(ActionEvent e) {
		fileChooser = new JFileChooser();
		fileChooser.setBounds(new Rectangle(width/2-20, 40, 100, 20));
		fileChooser.setDialogType(JFileChooser.FILES_AND_DIRECTORIES);
		fileChooser.setDialogTitle("Load data File for fitting procedure");
		fileChooser.setFileFilter(new TextFilter());
		
		int returnVal = fileChooser.showOpenDialog(this);
	    if(returnVal == JFileChooser.APPROVE_OPTION) {
	       System.out.println("You chose to open this file: " +
	            fileChooser.getSelectedFile().getName());
	    }
	    
	    this.data = readTextFile(fileChooser.getSelectedFile().getAbsolutePath());
	    

				
	  }
	
	/**
	 * Reads a text file and stores data in an array 
	 * @param nombreArchivo File name to read 
	 * @return data array
	 */
	private double[] readTextFile(String nombreArchivo){
		ArrayList<Double> data = new ArrayList<Double>();
		
		try{
			FileReader archivo = new FileReader(nombreArchivo);
			BufferedReader entrada = new BufferedReader(archivo);
			String s;
			StringTokenizer str;
											
			while (entrada.ready()){ 
				//Leer el número de instalaciones 
				s=entrada.readLine();
				str=new StringTokenizer (s);
				if (str.countTokens() != 1)throw new Exception ();
				data.add( new Double(Double.parseDouble(str.nextToken())) );
			}
			entrada.close();
			archivo.close();
			
		}catch(Exception e){
			System.out.println("Data file could not be read.");
			return null;
		}
        double[] datos = new double[data.size()];
		for(int i = 0; i < data.size();i++)datos[i] = data.get(i).doubleValue();
		return datos;
	}
	
	
	/**
	 * @param e Event
	 */
	void yesButton_actionPerformed(ActionEvent e) {
		  if(this.paramValues!=null && needPar){
			  for(int i = 0; i < this.paramValues.length; i++)
				  if(this.paramValues[i].getText().equals("")){
				  JOptionPane.showMessageDialog(null,  
						  "You must enter ALL the Parameters", 
						  "JPhase Alert", 
						  JOptionPane.INFORMATION_MESSAGE);
			  }else{
				  this.res = true;
				  this.setVisible(false);
			  }
		  }
		  else if(!needPar){
			  if(data == null){
				  JOptionPane.showMessageDialog(null,  
						  "You must select the data file", 
						  "JPhase Alert", 
						  JOptionPane.INFORMATION_MESSAGE);				  
			  }
			  else{
				  this.res = true;
				  this.paramValues = new JTextField[0];
				  this.setVisible(false);				  
			  }
		  }
	  }

	
	/**
	 * @param e Event
	 */
	  void noButton_actionPerformed(ActionEvent e) {
		  this.setVisible(false);
		  
	  }
	  
  
	  /**
	   * Centers the main frame 
	   */ 
	 private void centrarFrame() {
		  Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		  Dimension frameSize = this.getSize();
		  if (frameSize.height > screenSize.height) 
			  frameSize.height = screenSize.height;
		  if (frameSize.width > screenSize.width) 
			  frameSize.width = screenSize.width;
		  this.setLocation((screenSize.width - width) / 2 , 
				  (screenSize.height - height)/2);
	  }

	 	/**
	 	 * Action Adapter for fileButton
	 	 * @author Juan F. Perez
	 	 *
	 	 */
		class ParamFrame_fileButton_actionAdapter
	    			implements java.awt.event.ActionListener {
			
			/**
			 * 
			 */
			FitFrame adaptee;
			
			/**
			 * 
			 * @param adaptee
			 */
			ParamFrame_fileButton_actionAdapter(FitFrame adaptee) {
				this.adaptee = adaptee;
			}
		
			public void actionPerformed(ActionEvent e) {
				adaptee.fileButton_actionPerformed(e);
			}
		}
	 
	 
	 	/**
	 	 * Yes button action listener
	 	 * @author Juan F. Perez
	 	 *
	 	 */
		class ParamFrame_yesButton_actionAdapter
	    			implements java.awt.event.ActionListener {
			
			/**
			 * 
			 */
			FitFrame adaptee;
			
			/**
			 * 
			 * @param adaptee
			 */
			ParamFrame_yesButton_actionAdapter(FitFrame adaptee) {
				this.adaptee = adaptee;
			}
		
			public void actionPerformed(ActionEvent e) {
				adaptee.yesButton_actionPerformed(e);
			}
		}
		
		/**
		 * No button action listener
		 * @author Juan F. Perez
		 *
		 */
		class ParamFrame_noButton_actionAdapter
				implements java.awt.event.ActionListener {
			
			/**
			 * 
			 */
			FitFrame adaptee;
			
			/**
			 * 
			 * @param adaptee
			 */
			ParamFrame_noButton_actionAdapter(FitFrame adaptee) {
				this.adaptee = adaptee;
			}
			
			public void actionPerformed(ActionEvent e) {
				adaptee.noButton_actionPerformed(e);
			}
		}
	
		
		/**
		 * Text file filter 
		 */
		class TextFilter extends FileFilter {
			
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
				return "Archivos de texto(*.txt)";
			}
		}
		
		
		/**
		  * Returns file extension  
		  * @return String File extension 
		  * @param f File  
		  */ 
		private String getSuffix(File f) {
			String s = f.getPath(), suffix = null;
			int i = s.lastIndexOf('.');
			if (i>0 && i<s.length()-1) suffix = s.substring(i+1).toLowerCase();
			return suffix;
		}
}