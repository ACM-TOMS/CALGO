package jphase.GUI;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.JDialog;
import javax.swing.JRadioButton;
import javax.swing.ButtonGroup;

//import java.awt.color.*;
import jphase.PhaseVar;
import jphase.DenseContPhaseVar;

import jphase.fit.PhaseFitter;
import jphase.fit.MomentsACPH2Fit;
import jphase.fit.MomentsECCompleteFit;
import jphase.fit.MomentsECPositiveFit;
import jphase.fit.MomentsACPHFit;

import jphase.fit.EMHyperExpoFit;
import jphase.fit.EMHyperErlangFit;
import jphase.fit.EMPhaseFit;


/**
 * Frame to create a new PH variable
 * @author Juan F. Pérez
 * @version 1.0
 */
public class NewVarFrame extends JDialog implements ActionListener{
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
	 * Ask the user for a text input
	 */
	JTextField textInput;

	/**
	 * A theoretical variable or a fit procedure
	 */
	ButtonGroup  sourceType;

	/**
	 * 
	 */
	JRadioButton varSource;

	/**
	 * Selected if the source for the new variable is a fitting procedure
	 */
	JRadioButton fitSource;


	/**
	 * Type of theoretical variable
	 */
	JComboBox varType;

	/**
	 * Button to specify the parameters of the variable
	 */
	JButton paramVarButton;

	/**
	 * Type of fitting procedure
	 */
	JComboBox fitType;

	/**
	 * Button to specify the parameters of the fitting procedure
	 */
	JButton paramFitButton;

	/**
	 * Enter button
	 */
	JButton yesButton;

	/**
	 * Cancel Button
	 */
	JButton noButton;

	/**
	 * True if all the parameters are completely specified
	 */
	public boolean res = false;

	/**
	 * Fitting object
	 */
	public PhaseFitter fit;

	/**
	 * Fitted variable
	 */
	public PhaseVar var = null;

	/**
	 * Fitted Variable name
	 */
	public String varName;

	/**
	 * Frame width
	 */
	private int width = 400;

	/**
	 * Frame height
	 */
	private int height = 250;
	/**
	 * Construcción del frame principal
	 */
	public NewVarFrame(){

		this.setTitle("JPhase - New Variable Creation");
		this.setResizable(false);
		this.setModal(true);

		principalPanel = new JPanel();
		principalPanel.setLayout(null);
		principalPanel.setPreferredSize(new java.awt.Dimension(width, height));
		principalPanel.setOpaque(true);


		alertLabel= new JLabel("Enter the name of the new Variable");
		alertLabel.setBounds(new Rectangle(100, 10, 200, 20));
		alertLabel.setFont(new java.awt.Font("@Arial Unicode MS", 1, 12));
		alertLabel.setOpaque(true);
		principalPanel.add(alertLabel, null);

		textInput = new JTextField();
		textInput.setBounds(new Rectangle(110, 40, 150, 20));
		textInput.setBackground(Color.WHITE);
		textInput.setText("");
		principalPanel.add(textInput,BorderLayout.CENTER);


		sourceType = new ButtonGroup();

		varSource = new JRadioButton("Theoretical Variable");
		varSource.setOpaque(true);
		varSource.setBounds(new Rectangle(50, 80, 120, 20));
		varSource.setSelected(true);
		varSource.setActionCommand("SourceVar");
		varSource.addActionListener(this);

		fitSource = new JRadioButton("Fitted Variable");
		fitSource.setOpaque(true);
		fitSource.setBounds(new Rectangle(220, 80, 120, 20));
		fitSource.setActionCommand("SourceFit");
		fitSource.addActionListener(this);
		sourceType.add(varSource);
		sourceType.add(fitSource);
		principalPanel.add(varSource, null);
		principalPanel.add(fitSource, null);

		String[] vars = {"Expo","Erlang", "HyperExponential",
				"Coxian", "HyperErlang", "General Phase"}; 
		varType = new JComboBox(vars);
		varType.setBounds(new Rectangle(70, 110, 120, 20));
		principalPanel.add(varType, null);


		String[] fits = {"Moments ACPH2", "Moments EC Complete", 
				"Moments EC Positive",	"Moments ACPH", "EMHyperExpoFit", 
				"EMHyperErlangFit", "EMPhaseFit"}; 
		fitType = new JComboBox(fits);
		fitType.setBounds(new Rectangle(240, 110, 120, 20));
		fitType.setEnabled(false);
		principalPanel.add(fitType, null);


		paramVarButton = new JButton("Parameters");
		paramVarButton.setBounds(new Rectangle(80, 160, 90, 25));
		paramVarButton.addActionListener(this);
		paramVarButton.setActionCommand("Var");
		paramVarButton.setOpaque(true);
		principalPanel.add(paramVarButton, null);

		paramFitButton = new JButton("Parameters");
		paramFitButton.setBounds(new Rectangle(250, 160, 90, 25));
		paramFitButton.addActionListener(this);
		paramFitButton.setActionCommand("Fit");
		paramFitButton.setOpaque(true);
		paramFitButton.setEnabled(false);
		principalPanel.add(paramFitButton, null);




		yesButton = new JButton("Enter");
		yesButton.setBounds(new Rectangle(100, 220, 90, 25));
		yesButton.setActionCommand("Enter");
		yesButton.addActionListener(this);
		yesButton.setOpaque(true);
		principalPanel.add(yesButton, null);

		noButton = new JButton("Cancel");
		noButton.setBounds(new Rectangle(210, 220, 90, 25));
		noButton.setActionCommand("Cancel");
		noButton.addActionListener(this);
		principalPanel.add(noButton, null);

		//principalPanel.setBackground(Color.LIGHT_GRAY);
		this.getContentPane().add(principalPanel, BorderLayout.CENTER);
		this.centrarFrame();
		pack();

		//Screen size 
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		Dimension frameSize = this.getSize();
		if (frameSize.height > screenSize.height) frameSize.height = screenSize.height;
		if (frameSize.width > screenSize.width) frameSize.width = screenSize.width;
		this.setLocation((screenSize.width - frameSize.width) / 2, (screenSize.height - frameSize.height) / 2);
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
	 * @return Value of Text Field corresponding to Variable name
	 */
	public String getValue(){
		return this.textInput.getText();
	}

	/**
	 * 
	 * @param e
	 */
	void yesButton_actionPerformed( ) {
		if(this.textInput.getText().equals("") || !this.res){
			JOptionPane.showMessageDialog(null,  
					"You need to specify a name and all the parameters", 
					"JPhase Alert", 
					JOptionPane.INFORMATION_MESSAGE);
		}else{
			this.varName = this.textInput.getText();
			this.res = false;
			this.setVisible(false);
		}
	}

	/**
	 * 
	 * @param e
	 */
	void noButton_actionPerformed( ) {
		this.res = false;
		this.setVisible(false);
	}

	void paramVarButton_actionPerformed( ) {
		String reqVar = (String)this.varType.getSelectedItem();
		ParamFrame param = new ParamFrame(reqVar);
		param.setVisible(true);
		if(param.res){
			double[] p = new double[param.paramValues.length];
			for(int i =0; i < p.length; i++)p[i] = Double.valueOf(param.paramValues[i].getText()).doubleValue(); 
			if(reqVar.equals("Expo")){
				this.var = DenseContPhaseVar.expo(p[0]);
			}else if(reqVar.equals("Erlang")){
				this.var = DenseContPhaseVar.Erlang(p[0], (int)p [1]);
			}else if(reqVar.equals("Coxian")){
				int n = (int)p[0];
				double[] rates = new double[n];
				double[] probs = new double[n-1];
				System.arraycopy(p,1,rates,0,n );
				System.arraycopy(p,n+1,probs,0,n-1);
				this.var = DenseContPhaseVar.Coxian(n, rates, probs);
			}else if(reqVar.equals("HyperExponential")){
				int n = (int)p[0];
				double[] rates = new double[n];
				double[] probs = new double[n];
				System.arraycopy(p,1,rates,0,n );
				System.arraycopy(p,n+1,probs,0,n);
				this.var = DenseContPhaseVar.HyperExpo(rates, probs);
			}else if(reqVar.equals("HyperErlang")){
				int n = (int)p[0];
				double[] rates = new double[n];
				double[] probs = new double[n];
				double[] phases = new double[n];
				int[] phasesR = new int[n];
				System.arraycopy(p,1,rates,0,n );
				System.arraycopy(p,n+1,probs,0,n);
				System.arraycopy(p,2*n+1,phases,0,n);
				for (int i = 0; i != phases.length; i++){
					phasesR[i] = (int) phases[i];
				}
				this.var = DenseContPhaseVar.HyperErlang(n, rates, phasesR, probs);
			}
			else{
				System.out.println("NOT YET IMPLEMENTED");
			}
			if(this.var!= null)
				this.res = true;
			System.out.println("trulian " + res);
		} 
	}

	/**
	 * 
	 * @param e
	 */
	void paramFitButton_actionPerformed( ) {
		this.res = false;

		String reqFit = (String)this.fitType.getSelectedItem();

		FitFrame param = new FitFrame(reqFit);
		if(param.getRes()){

			double[] p = new double[param.paramValues.length];
			for(int i =0; i < p.length; i++)
				p[i] = Double.valueOf(param.paramValues[i].getText()).doubleValue();
			double[] data = param.getData();

			if(reqFit.equals("Moments ACPH2")){
				this.fit = new MomentsACPH2Fit(data);
				this.var = this.fit.fit();
			}else if(reqFit.equals("Moments EC Complete")){
				this.fit = new MomentsECCompleteFit(data);
				this.var = this.fit.fit();
			}else if(reqFit.equals("Moments EC Positive")){
				this.fit = new MomentsECPositiveFit(data);
				this.var = this.fit.fit();
			}else if(reqFit.equals("Moments ACPH")){
				this.fit = new MomentsACPHFit(data);
				this.var = this.fit.fit();
			}else if(reqFit.equals("EMHyperExpoFit")){
				this.fit = new EMHyperExpoFit(data);
				this.var = this.fit.fit();
			}else if(reqFit.equals("EMHyperErlangFit")){
				this.fit = new EMHyperErlangFit(data);
				this.var = this.fit.fit();
			}else if(reqFit.equals("EMPhaseFit")){
				this.fit = new EMPhaseFit(data);
				this.var = this.fit.fit();
			}else{
				System.out.println("NOT YET IMPLEMENTED");
			}
			this.varName = this.getValue();
			if(this.var!= null){
				this.res = true;
			}else{
				JOptionPane.showMessageDialog(null,  
						"Variable could not be fitted with the requested method\n More information in log", 
						"JPhase Alert", 
						JOptionPane.INFORMATION_MESSAGE);

				MainFrame.getInstance().addLog("Variable could not be fitted with the requested method");
				this.setVisible(false);

			}
		} 		  
	}


	/**
	 * 
	 * @param e
	 */
	void varSource_actionPerformed( ) {
		this.varType.setEnabled(true);
		this.varType.setForeground(Color.DARK_GRAY);
		this.paramVarButton.setEnabled(true);
		this.fitType.setEnabled(false);
		this.fitType.setForeground(Color.LIGHT_GRAY);
		this.paramFitButton.setEnabled(false);
	}


	/**
	 * 
	 * @param e
	 */
	void fitSource_actionPerformed( ) {
		this.varType.setEnabled(false);
		this.varType.setForeground(Color.LIGHT_GRAY);
		this.paramVarButton.setEnabled(false);
		this.fitType.setEnabled(true);
		this.fitType.setForeground(Color.DARK_GRAY);
		this.paramFitButton.setEnabled(true);
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
		//this.setLocation(100, 100);
	}	

	public void actionPerformed(ActionEvent arg0) {
		String evento = arg0.getActionCommand();

		if(evento.equals("Fit")){
			paramFitButton_actionPerformed();			
		}
		else if(evento.equals("Var")){
			paramVarButton_actionPerformed();
		}
		else if(evento.equals("SourceVar")){
			varSource_actionPerformed();
		}
		else if(evento.equals("SourceFit")){
			fitSource_actionPerformed();
		}
		else if(evento.equals("Enter")){
			yesButton_actionPerformed( );
		}
		else if(evento.equals("Cancel")){
			noButton_actionPerformed( );
		}
	}
}

