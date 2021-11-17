package jphase.GUI;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.Toolkit;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.border.LineBorder;

/**
 * <p>Title: FrameMini</p>
 * <p>Description: Graphic User Interface for JPhase</p>
 * <p>Copyright: Copyright (c) 2006-2014</p>
 * @author Juan F. Pérez
 * @version 1.0
 */
public class FrameMini extends JFrame {
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
	JTextArea alertLabel;
	
	/**
	 * Frame Width
	 */
	private int width = 250;
	
	/**
	 * Frame Height
	 */
	private int height = 150;
	
	/**
	 * Main Frame Construction 
	 * @param title Frame Title
	 * @param alert Frame Message 
	 */
	public FrameMini(String title, String alert){
		this.setTitle("JPhase - " +title);
		this.setIconImage(new ImageIcon("images/logo2.jpg").getImage());
		this.setResizable(false);
		
		principalPanel = new JPanel();
		principalPanel.setLayout(null);
		principalPanel.setPreferredSize(new java.awt.Dimension(width, height));
		
		alertLabel= new JTextArea(alert);
		alertLabel.setBounds(new Rectangle(10, 10, 230, 130));
		alertLabel.setFont(new java.awt.Font("@Arial Unicode MS", 1, 12));
		alertLabel.setBackground(new java.awt.Color(211 , 210, 192));
		alertLabel.setBorder(new LineBorder(Color.BLACK));
		principalPanel.add(alertLabel,null);
		
		this.getContentPane().add(principalPanel, BorderLayout.CENTER);
		this.centrarFrame();
		pack();
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
}

