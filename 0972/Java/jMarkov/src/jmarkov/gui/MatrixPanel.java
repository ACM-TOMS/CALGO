/*
 * Created on 3/09/2004
 */
package jmarkov.gui;

import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.ButtonGroup;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JRadioButton;


/**
 * @author Germán Riaño. Universidad de los Andes.
 *
 */
public class MatrixPanel extends InfoPanel 
implements ItemListener{

    private static final long serialVersionUID = 1969;

	private JRadioButton jRadioButtonDense = null;
	private JRadioButton jRadioButtonSparseRow = null;
	private JRadioButton jRadioButtonSparseRowsEvents = null;
	private JPanel jMatrixCardPanel = null;
	private JPanel ratioSelPanel = null;  //  @jve:decl-index=0:visual-constraint="12,21"
	private DenseRatesPanel denseRatesPanel = null;
	private SparseMatrixPanel sparseRatesPanelRows = null;
	private SparseMatrixPanel sparseRatesPanelRowsEvents = null;
	
	final static String DENSE = "Dense Matrix";
	final static String SPARSEROW = "Sparse by rows";
	final static String SPARSEROWEVENTS = "Sparse by rows with events";

	
	/**
	 * This method initializes the panel. 
	 * 
	 */
	public MatrixPanel() {
		super();
		initialize();
	}
	
	/**
	 * This method initializes this class
	 * 
	 */
	private void initialize() {
        this.setLayout(new BorderLayout());
        this.add(getJMatrixCardPanel(), java.awt.BorderLayout.CENTER);
        this.add(getRatioSelPanel(), java.awt.BorderLayout.NORTH);
	}
	
	/**
	 * This method initializes jRadioButtonDense	
	 * 	
	 * @return javax.swing.JRadioButton	
	 */    
	private JRadioButton getJRadioButtonDense() {
		if (jRadioButtonDense == null) {
			jRadioButtonDense = new JRadioButton();
			jRadioButtonDense.setText(DENSE);
			jRadioButtonDense.setSelected(true);
			jRadioButtonDense.setToolTipText("Shows the matrix in dense format");
			jRadioButtonDense.addItemListener(this);
		}
		return jRadioButtonDense;
	}
	/**
	 * This method initializes jRadioButtonSparseRow	
	 * 	
	 * @return javax.swing.JRadioButton	
	 */    
	private JRadioButton getJRadioButtonSparseRow() {
		if (jRadioButtonSparseRow == null) {
			jRadioButtonSparseRow = new JRadioButton();
			jRadioButtonSparseRow.setText(SPARSEROW);
			jRadioButtonSparseRow.setToolTipText("Shows the matrix in sparse format by rows");
			jRadioButtonSparseRow.addItemListener(this);
		}
		return jRadioButtonSparseRow;
	}
	/**
	 * This method initializes jRadioButtonSparseColumn	
	 * 	
	 * @return javax.swing.JRadioButton	
	 */    
	private JRadioButton getJRadioButtonSparseRowsEvents() {
		if (jRadioButtonSparseRowsEvents == null) {
			jRadioButtonSparseRowsEvents = new JRadioButton();
			jRadioButtonSparseRowsEvents.setText(SPARSEROWEVENTS);
			jRadioButtonSparseRowsEvents.setToolTipText("Shows the matrix in sparse format by rows, with Events");
			jRadioButtonSparseRowsEvents.setActionCommand("Sparse by rows, with events");
			jRadioButtonSparseRowsEvents.addItemListener(this);
		}
		return jRadioButtonSparseRowsEvents;
	}
	/**
	 * This method initializes jUpperPanel	
	 * 	
	 * @return javax.swing.JPanel	
	 */    
	private JPanel getJMatrixCardPanel() {
		if (jMatrixCardPanel == null) {
			jMatrixCardPanel = new JPanel();
			jMatrixCardPanel.setLayout(new CardLayout());
			denseRatesPanel = new DenseRatesPanel(); 
			sparseRatesPanelRows = new SparseMatrixPanel(); 
			sparseRatesPanelRowsEvents = new SparseMatrixPanel(true); 
			jMatrixCardPanel.add(denseRatesPanel,DENSE);
			jMatrixCardPanel.add(sparseRatesPanelRows,SPARSEROW);
			jMatrixCardPanel.add(sparseRatesPanelRowsEvents,SPARSEROWEVENTS);
		}
		return jMatrixCardPanel;
	}
	/**
	 * This method initializes jLowerPanel	
	 * 	
	 * @return javax.swing.JPanel	
	 */    
	private JPanel getRatioSelPanel() {
		if (ratioSelPanel == null) {
			ratioSelPanel = new JPanel();
			ratioSelPanel.add(getJRadioButtonDense(), null);
			ratioSelPanel.add(getJRadioButtonSparseRow(), null);
			ratioSelPanel.add(getJRadioButtonSparseRowsEvents(), null);
			ratioSelPanel.setSize(100,100);
			ButtonGroup group = new ButtonGroup();
		    group.add(jRadioButtonDense);
		    group.add(jRadioButtonSparseRow);
		    group.add(jRadioButtonSparseRowsEvents);
		}
		return ratioSelPanel;
	}
	
	public void itemStateChanged(ItemEvent evt) {
	    CardLayout cl = (CardLayout)(jMatrixCardPanel.getLayout());
	    cl.show(jMatrixCardPanel, ((JRadioButton)evt.getItem()).getText() );
	}

	
	
	/**
	 * Tests this class
	 * @param args
	 */
      	public static void main(String[] args) {
      		JFrame test = new JFrame();
      		test.setContentPane(new MatrixPanel());
      		test.setSize(850,750);
      		test.		addWindowListener(new WindowAdapter() {//listener to close window
    			@Override
					public void windowClosing(WindowEvent e) {
    				System.exit(0);
    				//me.setVisible(false);
    			}
    		});
      		test.setVisible(true);

	}
      	
      	
      	
	/** 
	 * @see jmarkov.gui.InfoPanel#updateMP()
	 */
	@Override
	public void updateMP() {
		denseRatesPanel.setMP(mp);
		sparseRatesPanelRows.setMP(mp);
		sparseRatesPanelRowsEvents.setMP(mp);
	}

}//@jve:decl-index=0:visual-constraint="10,10"
