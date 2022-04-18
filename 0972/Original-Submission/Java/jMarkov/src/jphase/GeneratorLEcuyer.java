/**
 * 
 */
package jphase;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;


import jphase.distributions.DistributionManager;
import jphase.distributions.IDistribution;
import jphase.distributions.UniformManager;
import jphase.values.FitResult;
import jphase.values.GenerationResult;
import jphase.values.ReaderWriter;

import cern.colt.list.DoubleArrayList;

import umontreal.iro.lecuyer.randvar.RandomVariateGen;
import umontreal.iro.lecuyer.rng.RandomStream;

/**
 * This class is in charge of perform the random variate generator
 * according the parameters provided by the used in the graphical
 * interface
 * @author Andres Sarmiento Romero
 * 
 */
public class GeneratorLEcuyer {

	/**
	 * Current random variated generator object
	 * according the actual distribution
	 */
	public RandomVariateGen generador;

	/**
	 * Current Uniform random variate generator
	 */
	public RandomStream randomStream;

	/**
	 * Current distribution to be generated
	 */
	public IDistribution distribucion;

	/**
	 * Name of the current Uniform random variate generator
	 */
	public String uniforme;
	
	/**
	 * System location to store the generated data
	 */
	public String ruta;

	/**
	 * Manager of the uniform random variate generators
	 */
	public UniformManager manejadorUni;

	/**
	 * Manager of the random variate generators
	 */
	public DistributionManager manejadorDist;

	/**
	 * Number of variates to generate
	 */
	public int n;

	/**
	 * Name of the generation procedure
	 */
	private String nombre;

	/**
	 * Constructor of the class
	 */
	public GeneratorLEcuyer() {
		manejadorUni = new UniformManager();
		manejadorDist = new DistributionManager();

		uniforme = "RandRijndael(Default)";

		randomStream = manejadorUni.darDistribucion(uniforme);		
		distribucion = manejadorDist.darDistribucion("Normal",null);		
		generador = distribucion.darGenerador(randomStream, null);

		n = 1000;
	}

	/**
	 * Generates the random variates according the class atributes 
	 * and returns a FitResults value 
	 */
	public FitResult generar( ) throws FileNotFoundException{		
		generador = distribucion.darGenerador( );

		DoubleArrayList lista = new DoubleArrayList();
		for (int i = 0; i != n; i++){
			double dub; 
			//			if(generador.toString().equals("LogNormalGen"))
			//				dub = Lognormal.nextDouble((LogNormal2Gen) generador);				
			//			else
			dub = generador.nextDouble();
			lista.add(dub);
		}

		//String nombre = Long.toString(System.currentTimeMillis());

		PrintWriter out = obtenerRuta(ruta, nombre +".txt");
		ReaderWriter.escritor(out, lista);

		FittingLEcuyer fit = new FittingLEcuyer();
		return fit.fit(10, ruta + "\\" + nombre + ".txt", distribucion.toString(), 95.0);		
	}

	/**
	 * @return the n
	 */
	public int getN() {
		return n;
	}

	/**
	 * @param n the n to set
	 */
	public void setN(int n) {
		this.n = n;
	}

	/**
	 * @return the generator
	 */
	public RandomVariateGen getGenerador() {
		return generador;
	}

	/**
	 * @param generador the generador to set
	 */
	public void setGenerador(RandomVariateGen generador) {
		this.generador = generador;
	}

	/**
	 * @return the distribuc
	 * iont
	 */
	public IDistribution getDistribucion() {
		return distribucion;
	}

	/**
	 * @param Xdistribucion the distribution to set
	 */
	public void setDistribucion(IDistribution Xdistribucion) {
		this.distribucion = Xdistribucion;
	}

	/**
	 * @return the uniform
	 */
	public String getUniforme() {
		return uniforme;
	}

	/**
	 * @param uniforme the uniform to set
	 */
	public void setUniforme(String uniforme) {
		this.uniforme = uniforme;
		randomStream = manejadorUni.darDistribucion(uniforme);
	}

	/**
	 * @return the ruta
	 */
	public String getRuta() {
		return ruta;
	}

	/**
	 * @param ruta the ruta to set
	 */
	public void setRuta(String ruta) {
		this.ruta = ruta;
	}

	public ArrayList<String> darGeneradoresUniformes() {
		return manejadorUni.darDistribuciones();
	}

	public ArrayList<String> darDistribuciones() {
		return manejadorDist.darDistribuciones();
	}

	public String darInfoUnif(String evento) {
		return manejadorUni.darInfoUnif(evento);
	}

	public IDistribution darDistribucion(String distT, ArrayList<Double> paramatros) {
		return manejadorDist.darDistribucion(distT, paramatros);
	}

	public void setDistribucion(String distri, ArrayList<Double> valores) {
		IDistribution iDis = manejadorDist.darDistribucion(distri, valores);
		setDistribucion(iDis);
	}

	public void setUbicacion(String ubicacion) {
		ruta = ubicacion;
	}

	public void setRandomStream( RandomStream rs){
		randomStream = rs;
	}

	public PrintWriter obtenerRuta(String ruta, String nombre) throws FileNotFoundException{

		File directorioFacturas = new File( ruta );
		if( !directorioFacturas.exists( ) )
			directorioFacturas.mkdirs( );
		File archivoFactura = new File( directorioFacturas, nombre );
		PrintWriter out = new PrintWriter( archivoFactura );

		return out;
	}

	public GenerationResult darInfo() {
		return new GenerationResult(distribucion.darDistribucion(), n, ruta);
	}

	/**
	 * Writes the moments of a distribution in a data file
	 */
	public double[] generateMoments() throws FileNotFoundException {
		double[] dob = distribucion.getMoments();
		//String nombre = Long.toString(System.currentTimeMillis());

		PrintWriter out = obtenerRuta(ruta, nombre +".txt");
		ReaderWriter.escritor(out, new DoubleArrayList(dob));
		
		return dob;
	}
	
	/**
	 * Sets the name of the distribution
	 */
	public void setNombre(String nombre){
		this.nombre = nombre;
	}
}
