package jphase.GUI;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.JDialog;

/**
 * <p>Title: InputFrame</p>
 * <p>Description: Graphic User Interface for JPhase</p>
 * <p>Copyright: Copyright (c) 2006</p>
 * @author Juan F. Pérez
 * @version 1.0
 */
public class InputFrame extends JDialog implements ActionListener{
				
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
	 * Accept Button
	 */
	JButton yesButton;
	
	/**
	 * Decline Button
	 */
	JButton noButton;
	
	/**
	 * Frame ready
	 */
	public boolean res = false;
	
	/**
	 * Frame width
	 */
	private int width = 250;
	
	/**
	 * Frame Height
	 */
	private int height = 140;
	
	/**
	 * Constructor of the main frame 
	 * @param title Frame Title
	 * @param alert Frame Message 
	 */
	public InputFrame(String title, String alert){
			
			
		this.setTitle("JPhase - " +title);
		this.setResizable(false);
		this.setModal(true);
		
		principalPanel = new JPanel();
		principalPanel.setLayout(null);
		principalPanel.setPreferredSize(new java.awt.Dimension(width, height));
		principalPanel.setOpaque(true);
		
		
		alertLabel= new JLabel(alert);
		alertLabel.setBounds(new Rectangle(40, 20, 200, 20));
		alertLabel.setFont(new java.awt.Font("@Arial Unicode MS", 1, 12));
		alertLabel.setOpaque(true);
		principalPanel.add(alertLabel, null);
		
		textInput = new JTextField();
		textInput.setBounds(new Rectangle(50, 60, 150, 20));
		textInput.setBackground(Color.WHITE);
		textInput.setText("");
		principalPanel.add(textInput,BorderLayout.CENTER);
		
		yesButton = new JButton("Enter");
		yesButton.setBounds(new Rectangle(20, 100, 90, 25));
		yesButton.setActionCommand("Enter");
		yesButton.addActionListener(this);
		yesButton.setOpaque(true);
		principalPanel.add(yesButton, null);
		
		noButton = new JButton("Cancel");
		noButton.setBounds(new Rectangle(130, 100, 90, 25));
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
	 * Returns True if all the info in the Frame is provided, False elsewhere
	 * @return True if all the info in the Frame is provided, False elsewhere
	 */
	public boolean getRes(){
		return this.res;
	}
	
	/**
	 * Returns text input given by the user
	 * @return Text input given by the user
	 */
	public String getValue(){
		return this.textInput.getText();
	}
	
	/**
	 * 
	 * @param e
	 */
	void yesButton_actionPerformed( ) {
		  if(this.textInput.getText().equals("")){
			  JOptionPane.showMessageDialog(null,  
					  "You must enter a name for the new Set", 
					  "JPhase Alert", 
					  JOptionPane.INFORMATION_MESSAGE);
		  }else{
			  this.res = true;
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


	public void actionPerformed(ActionEvent e) { 
		String message = e.getActionCommand();
		if (message.equals("Cancel")){
			noButton_actionPerformed();
		}else if (message.equals("Enter")){
			yesButton_actionPerformed();
		}
	}
}