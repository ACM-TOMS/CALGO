package jphase.GUI;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;

import jmarkov.jqbd.solvers.CLPAlgorithm;
import jmarkov.jqbd.solvers.ModBoundSolverFromMatrix;
import jmarkov.jqbd.solvers.MtjLogRedSolverFromMatrix;
import jmarkov.jqbd.solvers.QBDPhaseSolver;

public class NewQueueDialog extends JDialog implements ItemListener, ActionListener{
	public final String ACEPT = "Generate";	
	public final String CANCEL = "Cancel";

	public JButton aceptB;	
	public JButton cancelB;

	public ArrayList<PhaseVarInfo> variables;

	public JComboBox arrivals;
	public JComboBox services;

	public ButtonGroup algorith;

	public JRadioButton Ualg;
	public JRadioButton LR;
	public JRadioButton ModBon;
	public NewQueueDialog(MainFrame window){

		setTitle("JPhase - Construction of a PH/PH/1 model");

		algorith=new ButtonGroup();
		Ualg=new JRadioButton("Linear Progresion");
		algorith.add(Ualg);
		LR=new JRadioButton("Logarithmic Reduction");
		algorith.add(LR);
		ModBon=new JRadioButton("Modified Boundary");
		algorith.add(ModBon);
		LR.setSelected(true);

		aceptB = new JButton(ACEPT);
		aceptB.setActionCommand(ACEPT);
		aceptB.addActionListener( this );

		cancelB = new JButton(CANCEL);
		cancelB.setActionCommand(CANCEL);
		cancelB.addActionListener( this );
		setLayout(new BorderLayout());

		JPanel center = new JPanel();
		JPanel ncenter = new JPanel();
		JPanel algo = new JPanel();
		JPanel buttons = new JPanel();
		JPanel left = new JPanel();
		JPanel right = new JPanel();

		variables = window.getVariables();

		arrivals = new JComboBox( getNames(variables).toArray() );
		arrivals.setSelectedIndex( 0 );
		arrivals.addItemListener( this );
		add( arrivals );

		services = new JComboBox( getNames(variables).toArray() );
		services.setSelectedIndex( 0 );
		services.addItemListener( this );

		JLabel label1 = new JLabel("Please, select the Phase-Type distributions");
		JLabel label2 = new JLabel("that correspond to the distributions of");
		JLabel label3 = new JLabel("arrival and service times for a PH/PH/1/GD/∞/∞ model");
		label1.setHorizontalAlignment(javax.swing.SwingConstants.CENTER); 
		label2.setHorizontalAlignment(javax.swing.SwingConstants.CENTER); 
		label3.setHorizontalAlignment(javax.swing.SwingConstants.CENTER); 

		JLabel arrival = new JLabel("Arrivals");
		JLabel service = new JLabel("Services");      

		//left.setLayout(new GridLayout(2,1));
		left.add( arrival );
		left.add( arrivals );

		//right.setLayout(new GridLayout(2,1));
		right.add( service );
		right.add( services );   

		center.setLayout(new BorderLayout());
		center.add( left , BorderLayout.NORTH );
		center.add( right , BorderLayout.SOUTH );

		algo.setLayout(new GridLayout(3,1));
		algo.add( LR );
		algo.add( ModBon );
		algo.add( Ualg );

		ncenter.setLayout(new BorderLayout());
		ncenter.add( center , BorderLayout.WEST );
		ncenter.add( algo , BorderLayout.EAST );

		JPanel labels = new JPanel();
		labels.setLayout(new GridLayout(3,1));
		labels.add(label1);
		labels.add(label2);
		labels.add(label3);

		add( labels , BorderLayout.NORTH );
		add( ncenter , BorderLayout.CENTER );

		buttons.add( cancelB ); 
		buttons.add( aceptB );

		add( buttons , BorderLayout.SOUTH );

		setSize(400, 200);
		setResizable(true);
		setLocationRelativeTo (null);
	}

	public void itemStateChanged(ItemEvent arg0) {		
	}

	public void actionPerformed(ActionEvent arg0) {
		String message = arg0.getActionCommand();
		if(message.equals(ACEPT)){
			long ctime = System.currentTimeMillis();
			int arri = arrivals.getSelectedIndex();
			int serv = services.getSelectedIndex();
			QBDPhaseSolver test = null;
			if(Ualg.isSelected()){
				test = new CLPAlgorithm(variables.get(arri).var, variables.get(serv).var);
			}
			else if(ModBon.isSelected()){
				test = new ModBoundSolverFromMatrix(variables.get(arri).var, variables.get(serv).var);
			}
			else{
				test = new MtjLogRedSolverFromMatrix(variables.get(arri).var, variables.get(serv).var);
			}
			if(test.unstableSystem()){
				JOptionPane.showMessageDialog( this, "Unstable system", "Performance Measures", JOptionPane.ERROR_MESSAGE );
			}
			else{
				double[][] R = test.getRmatrix();
				test.printMatrices();
				String resp = test.performanceMeasures(R);
				ctime = System.currentTimeMillis() - ctime;
				resp += "\nTime duration: " + ctime + " milisecond"; 
				resp += (ctime!=1)?"s":"";
				JOptionPane.showMessageDialog( this, resp, "Performance Measures", JOptionPane.INFORMATION_MESSAGE );


				InputFrame probFrame = new InputFrame(
						"Steady State",
				"Number of states");
				probFrame.setVisible(true);
				probFrame.setFocusable(true);

				if(probFrame.res){
					try{
						int n = Integer.parseInt(probFrame.getValue());						
						ArrayList<Double> pis = test.getLevelSteadyStateProbs();				
						ArrayList<DenseMatrix> pisP = test.getSteadyStateProbsPerLevel();

						String ans = "";

						for(int i = 0 ;i < n && i < pis.size(); i++ ){
							ans += "Pi(" + i + ")=" + String.format("%6.4f", pis.get(i)) + "\t{";
							DenseMatrix doub = pisP.get(i);
							double[][] matrix = Matrices.getArray(doub);
							for(int j = 0 ;j < doub.numColumns(); j++ ){
								ans += String.format("%6.4f", matrix[0][j]) + "|";
							}
							ans += "}\n";
						}
						JOptionPane.showMessageDialog( this, ans, "Steady State Probabilities", JOptionPane.INFORMATION_MESSAGE );
					}
					catch(NumberFormatException e){

					}
				}
			}
		}
		else if(message.equals(CANCEL)){
			dispose();
		}
	}

	public ArrayList<String> getNames(ArrayList<PhaseVarInfo> vars){
		ArrayList<String> names = new ArrayList<String>();
		for(PhaseVarInfo var : vars)
			names.add(var.varName);
		return names;
	}

	public PhaseVarInfo getVar(String name){
		for(PhaseVarInfo var : variables){
			if(var.varName.equals(name))
				return var;			
		}
		return null;
	}

}
