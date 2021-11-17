package jphase.values;

import jphase.GUI.PhaseVarInfo;
import jphase.distributions.IDistribution;

/**
 * Class that stores the result of a 
 * goodness of fit test procedure
 * @author Andres Sarmiento Romero
 * 
 */
public class FitResult {
	
	/**
	 * Square Error of the GOF test
	 */
	private double sqrError;
	
	/**
	 * Chi-Square value result of the GOF test
	 */
	private double chi2;
	
	/**
	 * Kolmogorov-Smirnov value result of the GOF test
	 */
	private double ks;
	
	/**
	 * Number of groups used in the Chi-Square test
	 */
	private int groups;
	
	/**
	 * Name assigned to the GOF procedure
	 */
	private String nombre;
	
	/**
	 * Distribution used in the GOF test as initial hypothesis
	 */
	private IDistribution distribucion;
	
	/**
	 * If the GOF test was perform using a Phase-Type distribution
	 * in the initial hypothesis this atribute will store that
	 * distribution
	 */
	private PhaseVarInfo var;
	
	/**
	 * This value is used to set the percentage of data displayed
	 */
	private double p;
	
	/**
	 * Data array used in the procedure
	 */
	private double[] data;

	/**
	 * @param sqrError
	 * @param p p-value
	 * @param groups
	 */
	public FitResult(double sqrError, double chi2, double ks, 
			int groups, String nombre, IDistribution distribucion, double[] data, double p) {
		this.sqrError = sqrError;
		this.chi2 = chi2;
		this.ks = ks;
		this.groups = groups;
		this.data = data;
		this.distribucion = distribucion;
		this.nombre = nombre;
		this.p = p;
	}

	/**
	 * @param sqrError
	 * @param groups
	 * @param p 
	 */
	public FitResult(double sqrError, double chi2, int groups, PhaseVarInfo var, double[] data, double p) {
		this.sqrError = sqrError;
		this.chi2 = chi2;
		this.groups = groups;
		this.data = data;
		this.var = var;
		this.p = p;
	}

	/**
	 * @return the sqrError
	 */
	public double getSqrError() {
		return sqrError;
	}

	/**
	 * @param sqrError the sqrError to set
	 */
	public void setSqrError(double sqrError) {
		this.sqrError = sqrError;
	}

	/**
	 * @return the pValue
	 */
	public double getChi2() {
		return chi2;
	}

	/**
	 */
	public void setChi2(double chi2) {
		this.chi2 = chi2;
	}

	/**
	 * @return the ks
	 */
	public double getKs() {
		return ks;
	}

	/**
	 * @param ks the ks to set
	 */
	public void setKs(double ks) {
		this.ks = ks;
	}

	/**
	 * @return the groups
	 */
	public int getGroups() {
		return groups;
	}

	/**
	 * @param groups the groups to set
	 */
	public void setGroups(int groups) {
		this.groups = groups;
	}

	/**
	 * @return the name
	 */
	public String getNombre() {
		return nombre;
	}

	/**
	 * @param nombre the name to set
	 */
	public void setNombre(String nombre) {
		this.nombre = nombre;
	}

	/**
	 * @return the distribution
	 */
	public IDistribution getDistribucion() {
		return distribucion;
	}

	/**
	 * @param distribucion the distribution to set
	 */
	public void setDistribucion(IDistribution distribucion) {
		this.distribucion = distribucion;
	}

	/**
	 * @return the var
	 */
	public PhaseVarInfo getVar() {
		return var;
	}

	/**
	 * @param var the var to set
	 */
	public void setVar(PhaseVarInfo var) {
		this.var = var;
	}

	/**
	 * @return the data
	 */
	public double[] getData() {
		return data;
	}

	/**
	 * @param data the data to set
	 */
	public void setData(double[] data) {
		this.data = data;
	}
	
	/**
	 * @return the p
	 */
	public double getP() {
		return p;
	}

	/**
	 * @param p the p to set
	 */
	public void setP(double p) {
		this.p = p;
	}

	@Override
	public String toString(){
		return "Square Error: " + Math.rint(getSqrError()*1000)/1000+
		"\n\t(Reject if P-value < 0.05)"+
		"\n\tChi2 estimator: " + Math.rint(getChi2()*1000)/1000 + 
		"\n\tK-S estimator: " + Math.rint(getKs()*1000)/1000 + 
		"\nSaved in file: " + getNombre();		
	}
}
