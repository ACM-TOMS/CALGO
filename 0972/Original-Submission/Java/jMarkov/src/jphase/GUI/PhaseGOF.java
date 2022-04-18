package jphase.GUI;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.WindowConstants;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.xy.DefaultTableXYDataset;
import org.jfree.data.xy.XYSeries;

import jphase.AbstractContPhaseVar;
import jphase.FittingLEcuyer;
import jphase.PhaseVar;
import jphase.values.FitResult;

/**
 * This class is used to as the Interface for the Goodness of
 * Fit Tests for Phase-Type Variables
 * @author Andres Sarmiento Romero
 */
@SuppressWarnings("serial")
public class PhaseGOF extends JDialog implements ItemListener, ActionListener{
	
	/**
	 * String used as constant for the Select button
	 */
	private final String SELECT = "Select";
	
	/**
	 * String used as constant for the Acept button
	 */
	public final String ACEPT = "Acept";	
	
	/**
	 * String used as constant for the Close button
	 */
	public final String CANCEL = "Close";
	
	public JLabel sizeV;
	public JTextField cantF;

	public JLabel percentage;
	public JTextField perc;

	public JLabel save;
	public JButton saveB;

	public JButton ACEPTB;	
	public JButton CANCELB;

	private FittingLEcuyer fitting;

	private String ruta;

	public ArrayList<PhaseVarInfo> variables;

	public JComboBox varCombo;

	public PhaseGOF(MainFrame window){

		variables = window.getVariables();
		
		setTitle("JPhase - GOF for a Phase Type Variable");

		setSize(400,190); 
		setLocationRelativeTo(null);

		ACEPTB = new JButton(ACEPT);
		ACEPTB.setActionCommand(ACEPT);
		ACEPTB.addActionListener( this );

		CANCELB = new JButton(CANCEL);
		CANCELB.setActionCommand(CANCEL);
		CANCELB.addActionListener( this );


		fitting = new FittingLEcuyer();

		JPanel izq = new JPanel();
		JPanel der = new JPanel();

		sizeV = new JLabel("Number of groups for the Chi-Sqare test");
		cantF = new JTextField("10");

		percentage = new JLabel("Percentaje of data shown in the graph");
		perc = new JTextField("90");

		save = new JLabel("Select the file with the data");
		saveB = new JButton("Select");
		saveB.setActionCommand(SELECT);
		saveB.addActionListener( this );

		JLabel change = new JLabel("Select a Phase Variable");

		varCombo = new JComboBox( getNames(variables).toArray() );
		varCombo.setSelectedIndex( 0 );
		varCombo.addItemListener( this );

		setLayout(new BorderLayout());
		izq.setLayout(new GridLayout(4,1));
		der.setLayout(new GridLayout(4,1));

		izq.add(sizeV);
		der.add(cantF);
		izq.add(percentage);
		der.add(perc);
		izq.add(save);
		der.add(saveB);
		izq.add(change);
		der.add(varCombo);

		JPanel norte = new JPanel();
		norte.setLayout(new BorderLayout());

		JPanel sur = new JPanel();
		sur.add(CANCELB,BorderLayout.WEST);
		sur.add(ACEPTB,BorderLayout.EAST);

		norte.add(izq,BorderLayout.WEST);
		norte.add(der,BorderLayout.EAST);

		add(sur,BorderLayout.SOUTH);
		add(norte,BorderLayout.NORTH);
	}

	public void actionPerformed(ActionEvent e) {
		String evento = e.getActionCommand();
		if (evento.equals(ACEPT)){
			int serv = varCombo.getSelectedIndex();

			PhaseVarInfo variable = variables.get(serv);
			int n = 0;
			double p = 95;
			try{
				n = Integer.parseInt(cantF.getText());
				p = Double.parseDouble(perc.getText());
			}
			catch(NumberFormatException ex){
				JOptionPane.showMessageDialog( this, "Wrong number of groups", "Error", JOptionPane.ERROR_MESSAGE );
			}
			if(n<=0 || p <= 0 || p>100)
				JOptionPane.showMessageDialog( this, "Wrong number of groups", "Error", JOptionPane.ERROR_MESSAGE );
			else{
				FitResult results = fitting.fitPhases(variable, ruta, n, p);

				PhaseFitResultDialog res = new PhaseFitResultDialog(results);

				res.setVisible(true);
			}
		}
		else if (evento.equals(CANCEL))
			this.dispose();		
		else if (evento.equals(SELECT))
			selectFile();	
	}

	private void selectFile() {
		JFileChooser selector = new JFileChooser(); 
		//selector.setCurrentDirectory(new java.io.File("./data/"));
		selector.setDialogTitle("Select a File");
		selector.setFileSelectionMode(JFileChooser.FILES_ONLY);
		selector.setAcceptAllFileFilterUsed(false);

		if (selector.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) { 
			ruta = selector.getSelectedFile().toString();
		}
	}

	public void itemStateChanged(ItemEvent arg0) {		
	}

	public ArrayList<String> getNames(ArrayList<PhaseVarInfo> vars){
		ArrayList<String> names = new ArrayList<String>();
		for(PhaseVarInfo var : vars)
			names.add(var.varName);
		return names;
	}
}

/**
 * @author Andrés Sarmiento 
 *	This dialog is composed by the results of a phase type dialog
 */
class PhaseFitResultDialog extends JDialog implements ActionListener{
	
	private static final String ACEPTAR = "Ok"; 

	private static final String SAVE = "Save";

	private static double max; 
	
	private FitResult resultado;
	
	private JPanel estadisticos;
	
	private ChartPanel grafico;
	
	private JFreeChart chart;
	
	private DefaultTableXYDataset dataset2;
	
	private JButton aceptar;
	
	private JButton salvar;
	
	public PhaseFitResultDialog(FitResult resul){
		setDefaultCloseOperation( WindowConstants.DISPOSE_ON_CLOSE );
		resultado = resul;
		
		setTitle("JPhase - GOF for a Phase Type Variable");
		
		aceptar = new JButton(ACEPTAR);
		aceptar.setActionCommand(ACEPTAR);
		aceptar.addActionListener( this );
		
		salvar = new JButton(SAVE + " histogram");
		salvar.setActionCommand(SAVE);
		salvar.addActionListener( this );
		
		inicializarEstadisticos( );
		inicializarHistograma( );
		
		JPanel ajusteN = new JPanel();
		ajusteN.add(estadisticos,BorderLayout.WEST);
		ajusteN.add(grafico,BorderLayout.EAST);
		//ajusteN.add(grafico2,BorderLayout.CENTER);
		
		add(ajusteN,BorderLayout.NORTH);
		
		JPanel aux = new JPanel();
		aux.add(aceptar,BorderLayout.EAST);
		aux.add(salvar,BorderLayout.WEST);
		add(aux,BorderLayout.SOUTH);

		setSize(800, 400); 
		setLocationRelativeTo(null);
	}

	private void inicializarHistograma( ) {
		max = resultado.getP();
		HistogramDataset dataset = new HistogramDataset();
		//vecto almacena los ingresos quincenales de 45 personas
		double vector[] = darMax(resultado.getData());
		dataset.addSeries("Data", vector, resultado.getGroups());
		
		chart = ChartFactory.createHistogram( "Data", "x", "F(x)",
				dataset, PlotOrientation.VERTICAL, true, true, false);
		XYPlot plot = (XYPlot) chart.getPlot();
		XYBarRenderer renderer = (XYBarRenderer) plot.getRenderer();
		
		double min = chart.getXYPlot().getDomainAxis().getLowerBound();
//		chart.getXYPlot().getDomainAxis().setUpperBound(darMax(vector));
		double max = chart.getXYPlot().getDomainAxis().getUpperBound();
		
		PhaseVar var = resultado.getVar().var;
		
        XYSeries dataPDF = new XYSeries(resultado.getVar().varName, true, false);
		int numPoints = 200;
		double dx= (max-min) / numPoints;
		double x = min;
        for(int i = 0; i<numPoints; i++){
        	dataPDF.add(x, ((AbstractContPhaseVar) var).pdf(x));
        	x+=dx;
        }
        dataset2 = new DefaultTableXYDataset();
        dataset2.addSeries(dataPDF);
		
		renderer.setDrawBarOutline(false);
		
		JFreeChart chart2 =  ChartFactory.createXYLineChart("Phase Var " + resultado.getVar().varName + " vs Data",
        		"x", "F(x)", dataset2, PlotOrientation.VERTICAL, 
        		true, true, true);  
		
		XYPlot plot2 = (XYPlot) chart2.getPlot();
		ValueAxis domain2 = new NumberAxis( );
		
		plot2.setDataset(1,dataset);
		plot2.setRenderer(1,renderer);
		plot2.setRangeAxis(1,domain2);
		plot2.mapDatasetToRangeAxis(1, 1);
		
		grafico = new ChartPanel(chart2);
		chart = grafico.getChart();
		//grafico2 = new ChartPanel(charte);
		
		grafico.setPreferredSize(new java.awt.Dimension(500, 300)); 
		//grafico2.setPreferredSize(new java.awt.Dimension(300, 230)); 
	}

	private void inicializarEstadisticos( ) {
		estadisticos = new JPanel();
		estadisticos.setLayout(new GridLayout(6,1));
		estadisticos.add(new JLabel("Square Error: " + Math.rint(resultado.getSqrError()*1000)/1000));
		estadisticos.add(new JLabel("\t(Reject if P-value < 0.05)"));
		estadisticos.add(new JLabel("\tChi2 P-Value: " + Math.rint(resultado.getChi2()*1000)/1000));
	}

	public void actionPerformed(ActionEvent e) {
		String evento = e.getActionCommand();
		if(evento.equals(ACEPTAR))
			dispose();
		else if(evento.equals(SAVE))
			saveGraph();
	}
	
	public void saveGraph(){
		try{
			char sys[] = ("" + System.currentTimeMillis()).toCharArray();
			String name = ".\\data\\Histogram-" + resultado.getVar().varName + "-" + sys[sys.length-4] + sys[sys.length-3] + sys[sys.length-2] + sys[sys.length-1] + ".jpg";
			File file = new File(name);
			ChartUtilities.saveChartAsJPEG(file , chart, 500, 475);
	        JOptionPane.showMessageDialog( this, "The histogram was succesfully created:\n\t"+ name, "Done", JOptionPane.INFORMATION_MESSAGE );
		}
		catch(IOException e){
	        JOptionPane.showMessageDialog( this, "There were problems creating the histogram", "Error", JOptionPane.ERROR_MESSAGE );
		}		
	}

	public static double[] darMax (double[] unift){
		
        for( int i = unift.length; i > 0; i-- )
        {
            for( int j = 0; j < i - 1; j++ )
            {
                double p1 = unift[j];
                double p2 = unift[j + 1 ];

                if( p1 > p2 )
                {
                	unift[j] = p2;
                	unift[j + 1] =  p1;
                }
            }
        }
        int pos = (int)(unift.length * max/100);
        double list[] = new double[pos];
        for(int i = 0; i != pos; i++){
        	list[i] = unift[i];
        }
        return list;
//        return unift;
	}
}