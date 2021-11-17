package jphase.values;

import umontreal.iro.lecuyer.probdist.Distribution;

/**
 * This class stores the information generated in a Random Variate
 * Generation Process.
 * @author Andres Sarmiento Romero
 * 
 */
public class GenerationResult {
	
	/**
	 * Distribution used to generate the random variates
	 */
	private Distribution distribucion;
	
	/**
	 * Number of generated variates
	 */
	private int cantidad;
	
	/**
	 * Data file location where the data is going
	 * to be stored
	 */
	private String ruta;

	/**
	 * @param distribucion Distribution used to generate the variates
	 * @param cantidad Number of variates to generate
	 * @param ruta File location
	 */
	public GenerationResult(Distribution distribucion, int cantidad,
			String ruta) {
		this.distribucion = distribucion;
		this.cantidad = cantidad;
		this.ruta = ruta;
	}

	/**
	 * @return the distribution
	 */
	public Distribution getDistribucion() {
		return distribucion;
	}

	/**
	 * @param distribucion the distribution to set
	 */
	public void setDistribucion(Distribution distribucion) {
		this.distribucion = distribucion;
	}

	/**
	 * @return the cantidad
	 */
	public int getCantidad() {
		return cantidad;
	}

	/**
	 * @param cantidad the cantidad to set
	 */
	public void setCantidad(int cantidad) {
		this.cantidad = cantidad;
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

	@Override
	public String toString( ){
		return cantidad + " variables has been created\nwith " +
				distribucion + "\nEn " + ruta;
	}
}
