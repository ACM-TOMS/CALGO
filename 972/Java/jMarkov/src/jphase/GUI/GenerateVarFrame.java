/**
 * 
 */
package jphase.GUI;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;

import jphase.FittingLEcuyer;
import jphase.GeneratorLEcuyer;
import jphase.distributions.IDistribution;
import jphase.values.FitResult;
import jphase.values.GenerationResult;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.xy.DefaultTableXYDataset;


import umontreal.iro.lecuyer.randvar.RandomVariateGen;


/**
 * This Frame contains the variate generation and the Goodness Of fit tests tabs
 * @author Andrés Sarmiento Romero. Universidad de los Andes. (C) 2013
 */
public class GenerateVarFrame extends JFrame{

	/**
	 * Serializable id
	 */
	private static final long serialVersionUID = 6181837692334903139L;
	
	/**
	 * Panel of the variate generation
	 */
	public GeneratorPanel panelGenerador;
	
	/**
	 * Acceptance panel
	 */
	public AceptPanel panelAceptar;
	
	/**
	 * Panel of the goodness of fit tests
	 */
	public GOFPanel panelFOG;
	
	/**
	 * Variate Generator class
	 */
	public GeneratorLEcuyer generador;
	 
	/**
	 * Fitting class
	 */
	public FittingLEcuyer fitting;
	
	/**
	 * Tabbed panel
	 */
	public JTabbedPane tabble;

	/**
	 * Constructor of the class
	 */
	public GenerateVarFrame(){
		tabble = new JTabbedPane();
		generador = new GeneratorLEcuyer();
		fitting = new FittingLEcuyer();

		setTitle( "Goodness of Fit test and R.V. Generator" );
		setSize( new Dimension( 500, 255 ) );
		setResizable(true);

		panelGenerador = new GeneratorPanel(this, generador);
		panelFOG = new GOFPanel(this, fitting);
		panelAceptar = new AceptPanel(this);

		JPanel pan = new JPanel();
		pan.setLayout(new BorderLayout());
		pan.add(panelGenerador,BorderLayout.NORTH);
		pan.add(panelAceptar,BorderLayout.SOUTH);

		tabble.addTab("Generator", pan);
		tabble.addTab("Goodness Of Fit", panelFOG);

		add(tabble);

		setLocationRelativeTo (panelAceptar.aceptarB);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		GenerateVarFrame i = new GenerateVarFrame();
		i.setVisible(true);
	}

	/**
	 * Get the number of variates to be generated
	 * @return the number of variates to be generated
	 */
	public int darNumeroIteraciones() {
		int n = 0;
		try{
			n = Integer.parseInt(panelGenerador.darNumero());
		}
		catch (NumberFormatException e) {
		}
		return n;
	}

	/**
	 * Generates the variates
	 * @param n Number of variates to be generated
	 */
	public void generar(int n) {
		String ubicacion = panelGenerador.darUbicacion();
		if(ubicacion != null){
			generador.setN(n);
			generador.setUbicacion(ubicacion);
			generador.setNombre(panelGenerador.darNombre());
			try {
				FitResult result = generador.generar();
				GenerationResult gen = generador.darInfo( );
				if(result!= null){
					JDialog dia = new FitResultDialog(result);
					dia.setVisible(true);
					JOptionPane.showMessageDialog( this, gen.toString(), "Done", JOptionPane.INFORMATION_MESSAGE );
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
		else{
			JOptionPane.showMessageDialog( this, "Please select a location for the data", "Error", JOptionPane.ERROR_MESSAGE );
			panelGenerador.cambiarUbicacion();
		}
	}

	/**
	 * Generates the moments according the parameters
	 */
	public void generateMoments() {
		String ubicacion = panelGenerador.darUbicacion();
		System.out.println(ubicacion);
		if(ubicacion != null){
			generador.setUbicacion(ubicacion);
			try {
				double[] dob = generador.generateMoments();				 
				JOptionPane.showMessageDialog( this, "Moments:" +
						"\n\tMoment 1: " + Math.rint(dob[0]*1000)/1000 +
						"\n\tMoment 2: " + Math.rint(dob[1]*1000)/1000 + 
						"\n\tMoment 3: " + Math.rint(dob[2]*1000)/1000 , "Done", JOptionPane.INFORMATION_MESSAGE );
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
	}
}


/**
 * @author Andrés Sarmiento Romero. Universidad de los Andes. (C) 2013
 * This class content the Goodness Of Fit test Panel
 */
class GOFPanel extends JPanel implements ActionListener{

	/**        
	 * Serial id
	 */
    private static final long serialVersionUID = 1L;
    
    /**
     * Select Constant for the panel buttons
     */
    private final String SELECT = "Select";
    
    /**
     * Change Constant for the panel buttons
     */
	private final String CHANGE = "Change";
    
    /**
     * Accept Constant for the panel buttons
     */
	public final String ACEPT = "Acept";
    
    /**
     * Cancel Constant for the panel buttons
     */	
	public final String CANCEL = "Close";
    
    /**
     * Principal interface panel
     */
	public GenerateVarFrame interfaz;

	/**
	 * Dialog panel
	 */
	public JDialog dialogDist;

	/**
	 * Label for the number of variates
	 */
	public JLabel sizeV;
	
    /**
     *Input Label for the number of variates
     */
	public JTextField cantF;

    /**
     * Label for the percentage of variates
     */
	public JLabel percentage;
    
    /**
     *Input Label for the percentage of variates
     */
	public JTextField perc;

    /**
     * Label for the change of location
     */
	public JLabel change;
    
    /**
     * Button for the change of location
     */
	public JButton changeB;

	/**
	 * Save label
	 */
	public JLabel save;
	
	/**
	 * Save Button
	 */
	public JButton saveB;

	/**
	 * Accept button
	 */
	public JButton ACEPTB;
	
	/**
	 * Cancel button
	 */
	public JButton CANCELB;

	/**
	 * Fitting class
	 */
	private FittingLEcuyer fitting;

	/**
	 * Location of the file
	 */
	private String ruta;
	
	/**
	 * Distribution to be fitted
	 */
	private String distribucion;

	/**
	 * Constructor of the class
	 * @param ventanaX Principal panel
	 * @param fiting fitting class
	 */
	public GOFPanel(GenerateVarFrame ventanaX, FittingLEcuyer fiting){
		interfaz = ventanaX;		
		distribucion = "Normal";

		ACEPTB = new JButton(ACEPT);
		ACEPTB.setActionCommand(ACEPT);
		ACEPTB.addActionListener( this );

		CANCELB = new JButton(CANCEL);
		CANCELB.setActionCommand(CANCEL);
		CANCELB.addActionListener( this );


		fitting = fiting;

		JPanel izq = new JPanel();
		JPanel der = new JPanel();

		sizeV = new JLabel("Number of groups fot the Chi-Sqare test");
		cantF = new JTextField("10");

		percentage = new JLabel("Percentaje of data shown in the graph");
		perc = new JTextField("90");

		save = new JLabel("Select the file with the data");
		saveB = new JButton("Select");
		saveB.setActionCommand(SELECT);
		saveB.addActionListener( this );

		change = new JLabel("Select a Distribucion");
		changeB = new JButton("Change");
		changeB.setActionCommand(CHANGE);
		changeB.addActionListener( this );

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
		der.add(changeB);

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
			else
				fit(n, p);
		}
		else if (evento.equals(CANCEL))
			interfaz.dispose();		
		else if (evento.equals(SELECT))
			selectFile();	
		else if (evento.equals(CHANGE))
		    selectDistribution();	
	}

	/**
	 * Performs the fitting procedure
	 * @param n Number of variates
	 * @param p Percentage to be display
	 */
	private void fit(int n, double p) {
		if (ruta == null || distribucion == null){
			JOptionPane.showMessageDialog( this, "Select a file and a distribution", "Error", JOptionPane.ERROR_MESSAGE );
			if(ruta == null)
				selectFile();
			if (distribucion == null)
			    selectDistribution();
		}
		else{

			FitResult result = fitting.fit(n, ruta, distribucion, p);
			FitResultDialog dial = new FitResultDialog(result);
			dial.setVisible(true);
		}			
	}

	/**
	 * Displays a Distribution panel to select one
	 */
	private void selectDistribution() {
		if(dialogDist == null)
			dialogDist = new FitDistributionDialog(fitting, this);
		dialogDist.setVisible(true);
	}

	/**
	 * Displays a file chooser to select a location
	 */
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

	/**
     * Sets the new distribution
	 * @param distri Name of the distribution
	 */
	public void setDistribucion(String distri) {
		distribucion = distri;
	}

	/**
	 * Returns the name of the actaul distribution
	 * @return the name of the actual distribution
	 */
	public String getDistribucion( ) {
		return distribucion;
	}
}

/**
 * @author Andrés Sarmiento Romero. Universidad de los Andes. (C) 2013
 * This class content the Random Variate Generation Panel
 */
class GeneratorPanel extends JPanel implements ActionListener{

	/**
	 * Serialization id
	 */
	private static final long serialVersionUID = -3025450020220478855L;
	
	/**
	 * Constants for the labels
	 */
	public final static String cantidad = "Number of variables to generate: ";
	public final static String algU = "The uniform random variable generator is ";
	public final static String distC = "The distribution is set on ";
	public final static String genC = "The generator to use is ";
	public final static String CAMBIAR = "Change";
	public final static String UBICACION = "Specify a location to save the file";
	public final static String NAME = "Specify a file name";

	/**
	 * Principal interface
	 */
	public GenerateVarFrame ventana;

	/**
	 * Number of variables
	 */
	public JLabel cantidadV;
	public JTextField cantF;

	/**
	 * Name of the file
	 */
	public JLabel nombre;
	public JTextField nombreF;

	/**
	 * Uniform generator algorithm
	 */
	public JLabel algoUni;
	public JButton algoUniB;

	/**
	 * Distribution to be generated
	 */
	public JLabel dist;
	public JButton distB;

	/**
	 * Generator name
	 */
	public JLabel modGen;
	public JButton modGenB;

	/**
	 * Generator saving file
	 */
	public JLabel guardar;
	public JButton guardarB;

	/**
	 * Generator class
	 */
	public GeneratorLEcuyer generador;

	/**
	 * Uniforms distribution dialog
	 */
	public UniformDialog diaUnif;
	
	/**
	 * Distributions dialog
	 */
	public DistributionDialog diaDiast;
	
	/**
	 * File chooser selector
     */
	private JFileChooser selector;

	/**
	 * Possible location for the file
	 */
	public String ubicacion;

	/**
	 * Generator Panel constructor
	 * @param ventanaX Generator panel
	 * @param generador2 Generator class
	 */
	public GeneratorPanel(GenerateVarFrame ventanaX, GeneratorLEcuyer generador2){
		ubicacion = "./data/";

		ventana = ventanaX;
		generador = generador2;

		JPanel izq = new JPanel();
		JPanel der = new JPanel();

		cantidadV = new JLabel(cantidad);
		cantF = new JTextField("1000");

		nombre = new JLabel(NAME);
		nombreF = new JTextField("" + System.currentTimeMillis());

		algoUni = new JLabel(algU + generador.uniforme);
		algoUniB = new JButton(CAMBIAR);
		algoUniB.setActionCommand(CAMBIAR+"u");
		algoUniB.addActionListener( this );

		dist = new JLabel(distC + generador.distribucion);
		distB = new JButton(CAMBIAR);
		distB.setActionCommand(CAMBIAR+"d");
		distB.addActionListener( this );

		modGen = new JLabel(genC + generador.getGenerador().toString().split("with")[0].trim());
		modGenB = new JButton(CAMBIAR);
		modGenB.setActionCommand(CAMBIAR+"g");
		modGenB.addActionListener( this );

		guardar = new JLabel(UBICACION);
		guardarB = new JButton("Select");
		guardarB.setActionCommand(CAMBIAR+"s");
		guardarB.addActionListener( this );

		setLayout(new BorderLayout());
		izq.setLayout(new GridLayout(6,1));
		der.setLayout(new GridLayout(6,1));

		izq.add(cantidadV);
		der.add(cantF);
		izq.add(algoUni);
		der.add(algoUniB);
		izq.add(dist);
		der.add(distB);
		izq.add(modGen);
		der.add(modGenB);
		izq.add(guardar);
		der.add(guardarB);
		izq.add(nombre);
		der.add(nombreF);

		add(izq,BorderLayout.WEST);
		add(der,BorderLayout.EAST);
	}

	/**
	 * Name of the file
	 * @return name of the file
	 */
	public String darNombre() {
		String nombre = nombreF.getText();
		return (nombre.equals(""))?"" + System.currentTimeMillis():nombre;
	}

	public void actionPerformed(ActionEvent arg0) {
		String evento = arg0.getActionCommand();
		if(evento.equals(CAMBIAR+"u"))
			cambiarUni();
		else if(evento.equals(CAMBIAR+"d"))
			cambiarDistribucion();
		else if(evento.equals(CAMBIAR+"g"))
			cambiarGenerador();
		else if(evento.equals(CAMBIAR+"s"))
			cambiarUbicacion();
	}

	/**
	 * Changes the generator
	 */
	private void cambiarGenerador() {

		GeneratorDialog diaGen = new GeneratorDialog(generador, this);
		diaGen.setVisible(true);
		setEnabled(false);
	}

	/**
	 * Displays the dialog to change the distribution
	 */
	private void cambiarDistribucion() {
		if(diaDiast == null){
			diaDiast = new DistributionDialog(generador, this);
		}
		diaDiast.setVisible(true);
		this.setEnabled(false);

	}

	/**
	 * Displays a dialog to change the uniform generator
	 */
	private void cambiarUni() {
		if(diaUnif == null){
			diaUnif = new UniformDialog(generador, this);
		}
		diaUnif.setVisible(true);
		setEnabled(false);
	}

	/**
	 * Sets the uniform algorithm to be used
	 * @param algoUniT NBame of the algorithm
	 */
	public void setAlgoUni(String algoUniT) {
		generador.setUniforme(algoUniT);
		algoUni.setText(algU + algoUniT);
	}

	/**
	 * Sets the distribution
	 * @param distT Name of the distribution
	 * @param parametros Parameters of the distribution
	 */
	public void setDist(String distT, ArrayList<Double> parametros) {
		dist.setText(distC + distT);
		actualizarGen(distT, parametros);
	}

	/**
	 * Actualizes the generator
	 * @param distT Name of the distribution
	 * @param parametros parameters of the distribution
	 */
	private void actualizarGen(String distT, ArrayList<Double> parametros) {
		IDistribution dist = generador.darDistribucion(distT, parametros);
		RandomVariateGen g = dist.darGenerador(generador.randomStream, null);
		setGen(g.toString());
	}

	/**
	 * Sets the generator according the name
	 * @param modGenT Name of the generator
	 */
	public void setGen(String modGenT) {
		RandomVariateGen di = generador.getDistribucion().darGenerador(generador.randomStream,modGenT);
		modGen.setText(genC + (di.toString()).split("with")[0].trim());
	}

	/**
	 * Returns the generator object
	 * @return Generator
	 */
	public GeneratorLEcuyer getGenerator (){
		return generador;
	}
	
	@Override
    public String toString(){
		return "Generador de V.A.";
	}

	/**
	 * Returns the information of the uniform generator
	 * @param evento name of the generator
	 * @return the description of the generator
	 */
	public String darInfoUnif(String evento) {
		return generador.darInfoUnif(evento);
	}

	/**
	 * Enables the panel
	 */
	public void reEnable() {
		setEnabled(true);
	}

	/**
	 * Enables all the panels
	 * @param bool
	 */
	public void setEnable(Boolean bool){
		algoUniB.setEnabled(bool);
		distB.setEnabled(bool);
		modGenB.setEnabled(bool);
		guardarB.setEnabled(bool);
	}

	/**
	 * Returns the number of variates
	 * @return number of variates
	 */
	public String darNumero() {
		return cantF.getText();
	}

	/**
	 * This method calls a file chooser to select a new location
	 */
	public void cambiarUbicacion( ) {
		selector = new JFileChooser(); 
		selector.setCurrentDirectory(new java.io.File("./data/"));
		selector.setDialogTitle("Select the folder");
		selector.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		selector.setAcceptAllFileFilterUsed(false);

		if (selector.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) { 
			ubicacion = selector.getSelectedFile().toString();
		}
	}

	/**
	 * Returns the actual destination location of the file
	 * @return Actual location
	 */
	public String darUbicacion() {
		return ubicacion;
	}
}

/**
 * @author Andrés Sarmiento Romero. Universidad de los Andes. (C) 2013
 *	This dialog is composed by the available distributions to select
 */
class DistributionDialog extends JDialog implements ActionListener{

	/**
	 * Serialization long
	 */
    private static final long serialVersionUID = 1L;
    
    /**
     * Acceptance constants for the panel buttons
     */
    public final String ACEPTAR = "Acept";	
    
    /**
     * Cancel constant for the panel button
     */
	public final String CANCELAR = "Cancel";

	/**
	 * Generation panel
	 */
	public GeneratorPanel panel;

	/**
	 * Generator class
	 */
	public GeneratorLEcuyer gen;

	/**
	 * Group of buttons
	 */
	public ButtonGroup grupo;

	/**
	 * Acceptance button
	 */
	public JButton aceptarB;

	/**
	 * Canceling button
	 */
	public JButton cancelarB;

	/**
	 * Group of radio buttons
	 */
	public ArrayList<JRadioButton> botones;
	
	/**
	 * List of buttons
	 */
	public ArrayList<JButton> info;

	/**
	 * List of distributions
	 */
	public ArrayList<String> distribuciones;

	/**
	 * Name of actual distribution
	 */
	public String distri;

	/**
	 * Constructor of the class
	 * @param genle generation class
	 * @param pan generation panel
	 */
	public DistributionDialog(GeneratorLEcuyer genle , GeneratorPanel pan){
		panel = pan;		
		gen = genle;

		aceptarB = new JButton(ACEPTAR);
		aceptarB.setActionCommand(ACEPTAR);
		aceptarB.addActionListener( this );

		cancelarB = new JButton(CANCELAR);
		cancelarB.setActionCommand(CANCELAR);
		cancelarB.addActionListener( this );

		distribuciones = gen.darDistribuciones();

		botones = new ArrayList<JRadioButton>();

		grupo = new ButtonGroup();

		JPanel pan2 = new JPanel();

		pan2.setLayout(new GridLayout(distribuciones.size()+1,1));
		setSize(200,(distribuciones.size())*30+50);
		setResizable(false);

		setLayout(new BorderLayout());

		for(int i = 0; i != distribuciones.size(); i++){
			String nom = distribuciones.get(i);
			JRadioButton boton = new JRadioButton(nom);
			grupo.add(boton);
			botones.add(boton);
			pan2.add(boton);
			boton.setSelected(false);
			if(i == 0)
				boton.setSelected(true);
		}
		JPanel pan1 = new JPanel();
		pan1.setLayout(new GridLayout(1,2));
		pan1.setSize(250,85);
		pan1.add(cancelarB);
		pan1.add(aceptarB);
		add(pan2, BorderLayout.NORTH);
		add(pan1, BorderLayout.SOUTH);
		setLocationRelativeTo(null);
	}

	/**
	 * Generation performer
	 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	public void actionPerformed(ActionEvent e) {
		String evento = e.getActionCommand();
		if (evento.equals(ACEPTAR)){
			int j = 0;
			for(int i = 0; i != botones.size(); i++){
				if(botones.get(i).isSelected()){
					j = i;
					break;
				}					
			}
			distri = distribuciones.get(j);
			setEnabled(false);
			IDistribution iDist = gen.darDistribucion(distri, null);
			gen.setDistribucion(iDist);
			ArrayList<String> lista = iDist.darParametros();
			ParametersDialog diap = new ParametersDialog(lista, this);
			diap.setVisible(true);
		}
		else if (evento.equals(CANCELAR))
			dispose();			
	}

	/**
	 * Re-enables the panel
	 */
	public void reactivar() {
		setEnabled(true);
	}

	public void establecerParametros(ArrayList<Double> valores) {
		gen.setDistribucion(distri, valores);
		panel.setDist(gen.getDistribucion().toString(), valores);
		dispose();
	}

	@Override
    public void dispose(){
		panel.reEnable();
		super.dispose();
	}
}

/**
 *  @author Andrés Sarmiento Romero. Universidad de los Andes. (C) 2013
 *	This dialog is composed by the distributions to select for the GOF tests
 */
class FitDistributionDialog extends JDialog implements ActionListener{

	/**        
	 * Serialization long
	 */
    private static final long serialVersionUID = 1L;

    
    /**
     * Acceptance constants for the panel buttons
     */
    public final String ACEPTAR = "Acept";  
    
    /**
     * Cancel constant for the panel button
     */
    public final String CANCELAR = "Cancel";

    /**
     * Goodness of fit panel
     */
	public GOFPanel panel;

	/**
	 * Fitting class
	 */
	public FittingLEcuyer fit;

	/**
	 * Gorup of buttons
	 */
	public ButtonGroup grupo;

	/**
	 * Acceptance button
	 */
	public JButton aceptarB;

	/**
	 * Canceling button
	 */
	public JButton cancelarB;

	/**
	 * jRadio buttons
	 */
	public ArrayList<JRadioButton> botones;
	
	/**
	 * Information buttons
	 */
	public ArrayList<JButton> info;

	/**
	 * Distributions
	 */
	public ArrayList<String> distribuciones;

	/**
	 * Actual distribution
	 */
	public String distri;

	/**
	 * Constructor of the fitting panel
	 * @param fiti Fitting class
	 * @param pan Goodness of fit test panel
	 */
	public FitDistributionDialog(FittingLEcuyer fiti , GOFPanel pan){
		panel = pan;		
		fit = fiti;

		aceptarB = new JButton(ACEPTAR);
		aceptarB.setActionCommand(ACEPTAR);
		aceptarB.addActionListener( this );

		cancelarB = new JButton(CANCELAR);
		cancelarB.setActionCommand(CANCELAR);
		cancelarB.addActionListener( this );

		distribuciones = fit.darDistribuciones();

		botones = new ArrayList<JRadioButton>();

		grupo = new ButtonGroup();

		JPanel pan2 = new JPanel();

		pan2.setLayout(new GridLayout(distribuciones.size()+1,1));
		setSize(200,(distribuciones.size())*30+50);
		setResizable(false);

		setLayout(new BorderLayout());

		for(int i = 0; i != distribuciones.size(); i++){
			String nom = distribuciones.get(i);
			JRadioButton boton = new JRadioButton(nom);
			grupo.add(boton);
			botones.add(boton);
			pan2.add(boton);
			boton.setSelected(true);
		}
		JPanel pan1 = new JPanel();
		pan1.setLayout(new GridLayout(1,2));
		pan1.setSize(250,70);
		pan1.add(cancelarB);
		pan1.add(aceptarB);
		add(pan2, BorderLayout.NORTH);
		add(pan1, BorderLayout.SOUTH);
		setLocationRelativeTo(null);
	}

	public void actionPerformed(ActionEvent e) {
		String evento = e.getActionCommand();
		if (evento.equals(ACEPTAR)){
			int j = 0;
			for(int i = 0; i != botones.size(); i++){
				if(botones.get(i).isSelected()){
					j = i;
					break;
				}					
			}
			distri = distribuciones.get(j);
			panel.setDistribucion(distri);
		}
		dispose();
	}
}

/**
 *  @author Andrés Sarmiento Romero. Universidad de los Andes. (C) 2013
 *	This dialog shows the result of the GOF test
 */
class FitResultDialog extends JDialog implements ActionListener{

	/**        */
    private static final long serialVersionUID = 1L;

    /**
     * Acceptance constant fo the button
     */
    private static final String ACEPTAR = "Ok"; 

    /**
     * Saving constant for the panel buttons 
     */
	private static final String SAVE = "Save";

	/**
	 * Maximum display bars
	 */
	private static double max; 

	/**
	 * Value class that drives the result of the fitting procedure
	 */
	private FitResult resultado;

	/**
	 * Statistics summary panel
	 */
	private JPanel estadisticos;

	/**
	 * Graph of the histogram
	 */
	private ChartPanel grafico;

	/**
	 * Graph of the PDF 
	 */
	private ChartPanel grafico2;

	/**
	 * Chart to be shown in the panel
	 */
	private JFreeChart chart;

	/**
	 * Data set
	 */
	private DefaultTableXYDataset dataset2;

	/**
	 * Acceptance button
	 */
	private JButton aceptar;

	/**
	 * Saving button
	 */
	private JButton salvar;

	/**
	 * Constructor of the class
	 * @param resul Object that mangaes the summary of the fitting procedure
	 */
	public FitResultDialog(FitResult resul){
		setDefaultCloseOperation( JFrame.DISPOSE_ON_CLOSE );
		resultado = resul;

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

	/**
	 * Initializes the histogram (put it in a nice way) 
	 */
	private void inicializarHistograma( ) {
		max = resultado.getP();
		HistogramDataset dataset = new HistogramDataset();
		//		dataset.
		//vecto almacena los ingresos quincenales de 45 personas
		double vector[] = darMax(resultado.getData());
		dataset.addSeries("Data (Histogram)", vector, resultado.getGroups());

		chart = ChartFactory.createHistogram( "Data (Histogram)", "x", "f(x)",
				dataset, PlotOrientation.VERTICAL, true, true, false);
		XYPlot plot = (XYPlot) chart.getPlot();
		XYBarRenderer renderer = (XYBarRenderer) plot.getRenderer();

		double min = chart.getXYPlot().getDomainAxis().getLowerBound();
		double max = chart.getXYPlot().getDomainAxis().getUpperBound();

		dataset2 = resultado.getDistribucion().getPDF(min, max, 200);

		renderer.setDrawBarOutline(false);

		JFreeChart chart2 =  ChartFactory.createXYLineChart("Fit result",
				"x", "f(x)", dataset2, PlotOrientation.VERTICAL, 
				true, true, true);  
		XYPlot plot2 = (XYPlot) chart2.getPlot();
		XYItemRenderer renderer2 = plot2.getRenderer();
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

	/**
	 * Summaries the statistical information to be shown
	 */
	private void inicializarEstadisticos( ) {
		estadisticos = new JPanel();
		estadisticos.setLayout(new GridLayout(6,1));
		estadisticos.add(new JLabel("Square Error: " + Math.rint(resultado.getSqrError()*1000)/1000));
		estadisticos.add(new JLabel("\t(Reject if P-value < 0.05)"));
		estadisticos.add(new JLabel("\tChi2 P-Value: " + Math.rint(resultado.getChi2()*1000)/1000));
		estadisticos.add(new JLabel("\tK-S P-Value: " + Math.rint(resultado.getKs()*1000)/1000));
		estadisticos.add(new JLabel("File name: " + resultado.getNombre()));
		IDistribution distrib = resultado.getDistribucion();
		setTitle(distrib.aString());
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
			File file = new File(".\\data\\histogram" + resultado.getNombre() + ".jpg");
			ChartUtilities.saveChartAsJPEG(file , chart, 500, 475);
			JOptionPane.showMessageDialog( this, "The histogram was succesfully created:\n\t"+ file.getAbsolutePath(), "Done", JOptionPane.INFORMATION_MESSAGE );
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

/**
 *  @author Andrés Sarmiento Romero. Universidad de los Andes. (C) 2013 
 *	This panel contains the buttons used to generate variates, moments or cancel
 */
class AceptPanel extends JPanel implements ActionListener{

	/**        
	 * Serial version
	 */
    private static final long serialVersionUID = 1L;
    
    /**
     * Constant used to generate the variables
     */
    public final String GENVARIATES = "Generate Variables";	
    
    /**
     * Constant used to generate the moments
     */
	public final String GENMOMENTS = "Generate Moments";
	
	/**
	 * Constant to close the panel
	 */
	public final String CANCELAR = "Close";

	/**
	 * Generation interface
	 */
	public GenerateVarFrame interfaz;

	/**
	 * Var Generation button button
	 */
	public JButton aceptarB;

	/**
	 * Moment generation moment
	 */
	public JButton aceptarB2;

	/**
	 * Cancel button
	 */
	public JButton cancelarB;

	/**
	 * Constructor of the class
	 * @param pan Principal panel
	 */
	public AceptPanel(GenerateVarFrame pan){
		interfaz = pan;		

		aceptarB = new JButton(GENVARIATES);
		aceptarB.setActionCommand(GENVARIATES);
		aceptarB.addActionListener( this );

		aceptarB2 = new JButton(GENMOMENTS);
		aceptarB2.setActionCommand(GENMOMENTS);
		aceptarB2.addActionListener( this );

		cancelarB = new JButton(CANCELAR);
		cancelarB.setActionCommand(CANCELAR);
		cancelarB.addActionListener( this );

		add(cancelarB);
		add(aceptarB);
		add(aceptarB2);
	}

	public void actionPerformed(ActionEvent e) {
		String evento = e.getActionCommand();
		if (evento.equals(GENVARIATES)){
			int n = interfaz.darNumeroIteraciones();
			if(n<=0)
				JOptionPane.showMessageDialog( this, "Wrong number of variables", "Number error", JOptionPane.ERROR_MESSAGE );
			else
				interfaz.generar(n);
		}
		if (evento.equals(GENMOMENTS)){
			interfaz.generateMoments( );
		}
		else if (evento.equals(CANCELAR))
			interfaz.dispose();		
	}
}

/**
 *  @author Andrés Sarmiento Romero. Universidad de los Andes. (C) 2013
 *  This dialog is useg to generate the variables
 */
class GeneratorDialog extends JDialog implements ActionListener{

    /**
     * Serialization long
     */
    private static final long serialVersionUID = 1L;
    
    /**
     * Acceptance constants for the panel buttons
     */
    public final String ACEPTAR = "Acept";  
    
    /**
     * Cancel constant for the panel button
     */
    public final String CANCELAR = "Cancel";

    /**
     * Generation panel
     */
    public GeneratorPanel panel;

    /**
     * Generator class
     */
    public GeneratorLEcuyer gen;

    /**
     * Group of buttons
     */
    public ButtonGroup grupo;

    /**
     * Acceptance button
     */
    public JButton aceptarB;

    /**
     * Canceling button
     */
    public JButton cancelarB;

    /**
     * Group of radio buttons
     */
    public ArrayList<JRadioButton> botones;
    
    /**
     * List of generators 
     */
	public ArrayList<String> generadores;

	/**
	 * Constructor of the panel
	 * @param genle Generator class
	 * @param pan Generator panel
	 */
	public GeneratorDialog(GeneratorLEcuyer genle , GeneratorPanel pan){
		panel = pan;		
		gen = genle;


		aceptarB = new JButton(ACEPTAR);
		aceptarB.setActionCommand(ACEPTAR);
		aceptarB.addActionListener( this );

		cancelarB = new JButton(CANCELAR);
		cancelarB.setActionCommand(CANCELAR);
		cancelarB.addActionListener( this );

		IDistribution dist = gen.getDistribucion();

		generadores = dist.darGeneradores();

		botones = new ArrayList<JRadioButton>();

		grupo = new ButtonGroup();

		setLayout(new BorderLayout());
		JPanel p1 = new JPanel();
		p1.setLayout(new GridLayout(generadores.size()+1,1));
		setSize(200,(generadores.size()*30)+50);
		setResizable(false);

		for(int i = 0; i != generadores.size(); i++){
			String nom = generadores.get(i);
			JRadioButton boton = new JRadioButton(nom);
			grupo.add(boton);
			botones.add(boton);
			p1.add(boton);
			if(i == 0)
				boton.setSelected(true);
		}
		JPanel p = new JPanel();
		p.setLayout(new BorderLayout());

		p.add(cancelarB, BorderLayout.WEST);
		p.add(aceptarB, BorderLayout.EAST);

		add(p, BorderLayout.SOUTH);
		add(p1, BorderLayout.NORTH);
		setLocationRelativeTo(null);
		setTitle("Generators");
	}

	public void actionPerformed(ActionEvent e) {
		String evento = e.getActionCommand();
		if (evento.equals(ACEPTAR)){
			int j = 0;
			for(int i = 0; i != botones.size(); i++){
				if(botones.get(i).isSelected()){
					j = i;
					break;
				}					
			}
			panel.setGen(generadores.get(j));
			dispose();
		}
		else if (evento.equals(CANCELAR))
			this.dispose();
		else
			darInfo(evento);			
	}

	/**
	 * Returns information of the uniform distribution
	 * @param evento
	 */
	public void darInfo(String evento) {
		String info = panel.darInfoUnif(evento);
		JOptionPane.showMessageDialog( this, info, evento + info, JOptionPane.INFORMATION_MESSAGE );
	}

	@Override
    public void dispose(){
		panel.reEnable();
		super.dispose();
	}
}


/**
 * @author Andrés Sarmiento Romero. Universidad de los Andes. (C) 2013
 *  This dialog is composed by the available uniform distributions to select
 */
class UniformDialog extends JDialog implements ActionListener{

    /**
     * Serialization long
     */
    private static final long serialVersionUID = 1L;
    
    /**
     * Acceptance constants for the panel buttons
     */
    public final String ACEPTAR = "Acept";  
    
    /**
     * Cancel constant for the panel button
     */
    public final String CANCELAR = "Cancel";

    /**
     * Generation panel
     */
    public GeneratorPanel panel;

    /**
     * Generator class
     */
    public GeneratorLEcuyer gen;

    /**
     * Group of buttons
     */
    public ButtonGroup grupo;

    /**
     * Acceptance button
     */
    public JButton aceptarB;

    /**
     * Canceling button
     */
    public JButton cancelarB;

    /**
     * Group of radio buttons
     */
    public ArrayList<JRadioButton> botones;
    
    /**
     * List of buttons
     */
    public ArrayList<JButton> info;

    /**
     * List of uniform distributions
     */
	public ArrayList<String> uniformes;

	/**
	 * constructor of the class
	 * @param genle Generation class
	 * @param pan Generation panel
	 */
	public UniformDialog(GeneratorLEcuyer genle , GeneratorPanel pan){
		panel = pan;		
		gen = genle;

		aceptarB = new JButton(ACEPTAR);
		aceptarB.setActionCommand(ACEPTAR);
		aceptarB.addActionListener( this );

		cancelarB = new JButton(CANCELAR);
		cancelarB.setActionCommand(CANCELAR);
		cancelarB.addActionListener( this );

		uniformes = gen.darGeneradoresUniformes();

		botones = new ArrayList<JRadioButton>();

		grupo = new ButtonGroup();

		setSize(200,(uniformes.size())*30+20);
		setResizable(false);

		JPanel izq = new JPanel();
		izq.setLayout(new GridLayout(uniformes.size()+1,1));
		JPanel der = new JPanel();
		der.setLayout(new GridLayout(uniformes.size()+1,1));

		for(int i = 0; i != uniformes.size(); i++){
			String nom = uniformes.get(i);
			JRadioButton boton = new JRadioButton(nom);
			grupo.add(boton);
			botones.add(boton);
			izq.add(boton);
			boton.setSelected(false);
			if(i == 0)
				boton.setSelected(true);

			JButton bot = new JButton("i");
			bot.setActionCommand(nom);
			bot.addActionListener( this );
			der.add(bot);
		}


		JPanel sur = new JPanel();
		sur.add(cancelarB);
		sur.add(aceptarB);
		setLayout(new BorderLayout());
		add(sur,BorderLayout.SOUTH);
		add(izq,BorderLayout.WEST);
		add(der,BorderLayout.EAST);

		setLocationRelativeTo(null);		
		setDefaultCloseOperation( JFrame.DISPOSE_ON_CLOSE );
		setTitle("Uniform variable generators");
	}

	public void actionPerformed(ActionEvent e) {
		String evento = e.getActionCommand();
		if (evento.equals(ACEPTAR)){
			int j = 0;
			for(int i = 0; i != botones.size(); i++){
				if(botones.get(i).isSelected()){
					j = i;
					break;
				}					
			}
			panel.setAlgoUni(uniformes.get(j));
			this.dispose();
		}
		else if (evento.equals(CANCELAR))
			this.dispose();
		else
			darInfo(evento);			
	}

	/**
	 * REturns a description of the uniform generator
	 * @param evento
	 */
	public void darInfo(String evento) {
		String info = panel.darInfoUnif(evento);
		JOptionPane.showMessageDialog( this, info, evento + " info", JOptionPane.INFORMATION_MESSAGE );
	}

	@Override
    public void dispose(){
		panel.setEnable(true);
		panel.reEnable();
		super.dispose();
	}
}


/**
 * @author Andrés Sarmiento Romero. Universidad de los Andes. (C) 2013
 *  This dialog displays the availeable parametrs of a distribution
 */
class ParametersDialog extends JDialog implements ActionListener{
    
    /**
     * Serialization long
     */
    private static final long serialVersionUID = 1L;
    
    /**
     * Acceptance constants for the panel buttons
     */
    public final String ACEPTAR = "Acept";  
    
    /**
     * Cancel constant for the panel button
     */
    public final String CANCELAR = "Cancel";

    /**
     * Acceptance button
     */
    public JButton aceptarB;

    /**
     * Canceling button
     */
    public JButton cancelarB;

    /**
     * List of the names of the paramets
     */
	public ArrayList<String> nombre;

	/**
	 * List of labels that models the parameters
	 */
	public ArrayList<JLabel> labels;

	/**
	 * Input labels for the parameters
	 */
	public ArrayList<JTextField> textos;
	
	/**
	 * Distribution dialog
	 */
	private DistributionDialog dialogo;

	/**
	 * Constructor of the class
	 * @param parametros List of the parametrs of the distribution
	 * @param dialogoX Dialog panel
	 */
	public ParametersDialog(ArrayList<String> parametros, DistributionDialog dialogoX){
		dialogo = dialogoX;

		setSize(200,(parametros.size())*30+50);
		setResizable(false);
		setLayout(new BorderLayout());
		nombre = parametros;

		JPanel panUp = new JPanel();
		panUp.setLayout(new GridLayout(parametros.size(),2));

		labels = new ArrayList<JLabel>();
		textos = new ArrayList<JTextField>();

		for(String param: nombre){
			JLabel label = new JLabel(param + ": ");
			labels.add(label);
			panUp.add(label);
			JTextField texto = new JTextField();
			textos.add(texto);
			panUp.add(texto);
		}

		aceptarB = new JButton(ACEPTAR);
		aceptarB.setActionCommand(ACEPTAR);
		aceptarB.addActionListener( this );

		cancelarB = new JButton(CANCELAR);
		cancelarB.setActionCommand(CANCELAR);
		cancelarB.addActionListener( this );

		JPanel panDo = new JPanel();
		panDo.setLayout(new BorderLayout());
		panDo.add(cancelarB, BorderLayout.WEST);
		panDo.add(aceptarB, BorderLayout.EAST);

		add(panUp,BorderLayout.NORTH);
		add(panDo,BorderLayout.SOUTH);

		setLocationRelativeTo(null);
		setTitle("Parameters");
	}

	public void actionPerformed(ActionEvent e) {
		String evento = e.getActionCommand();
		if (evento.equals(ACEPTAR)){
			ArrayList<Double> valores = new ArrayList<Double>();
			boolean number = true;
			for(JTextField texto : textos){
				try{
					double valor = Double.parseDouble(texto.getText());
					valores.add(valor);
				}
				catch(NumberFormatException e1){
					number = false;
					break;
				}
			}
			if(number){
				dialogo.establecerParametros(valores);
				this.dispose();
			}
			else
				JOptionPane.showMessageDialog( this, "There are problems with some parameters", "Error", JOptionPane.ERROR_MESSAGE );
		}
		else if (evento.equals(CANCELAR)){
			dialogo.reactivar();
			this.dispose();			
		}
	}

	@Override
    public void dispose(){
		dialogo.reactivar();
		super.dispose();
	}

}
