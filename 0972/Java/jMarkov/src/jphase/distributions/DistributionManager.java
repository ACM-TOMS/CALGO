package jphase.distributions;

import java.util.ArrayList;

import umontreal.iro.lecuyer.probdist.ContinuousDistribution;
import umontreal.iro.lecuyer.probdist.ErlangDist;
import umontreal.iro.lecuyer.probdist.ExponentialDist;
import umontreal.iro.lecuyer.probdist.GammaDist;
import umontreal.iro.lecuyer.probdist.LognormalDist;
import umontreal.iro.lecuyer.probdist.NormalDist;
import umontreal.iro.lecuyer.probdist.UniformDist;
import umontreal.iro.lecuyer.probdist.WeibullDist;
import umontreal.iro.lecuyer.randvar.ErlangGen;
import umontreal.iro.lecuyer.randvar.ExponentialGen;
import umontreal.iro.lecuyer.randvar.GammaGen;
import umontreal.iro.lecuyer.randvar.LognormalGen;
import umontreal.iro.lecuyer.randvar.NormalACRGen;
import umontreal.iro.lecuyer.randvar.NormalBoxMullerGen;
import umontreal.iro.lecuyer.randvar.NormalKindermannRamageGen;
import umontreal.iro.lecuyer.randvar.NormalPolarGen;
import umontreal.iro.lecuyer.randvar.RandomVariateGen;
import umontreal.iro.lecuyer.randvar.UniformGen;
import umontreal.iro.lecuyer.randvar.WeibullGen;
import umontreal.iro.lecuyer.rng.RandomStream;
import cern.colt.list.DoubleArrayList;

/**
 * This class manages the distribution used to be fitted or generated
 * @author Andrés Sarmiento Romero. Universidad de los Andes. (C) 2013
 */
public class DistributionManager {

	/**
	 * List of Distributions
	 */
	public ArrayList<String> distribuciones;
	
	/**
	 * Initializes the distributions
	 */
	public DistributionManager(){
		distribuciones = new ArrayList<String>();
		distribuciones.add("Normal");
		distribuciones.add("Lognormal");
		distribuciones.add("Gamma");
		distribuciones.add("Weibull");
		distribuciones.add("Uniform");
		distribuciones.add("Exponential");
		distribuciones.add("Erlang");
	}
	
	/**
	 * Returns the list of distributions
	 * @return List of distributions
	 */
	public ArrayList<String> darDistribuciones() {
		return distribuciones;
	}
	
	/**
	 * Returns the distribution given the name and the distribution
	 * @param distribucion Name of the distribution
	 * @param parametros Parameters of the distribution
	 * @return The distribution
	 */
	public IDistribution darDistribucion(String distribucion, ArrayList<Double> parametros) {
		IDistribution cont;
		String distO = distribucion.split(":")[0].trim();
		if(distO.equals("NormalDist") || distribucion.equals("Normal"))
			cont = new Normal( parametros );
		else if(distO.equals("LognormalDist") || distribucion.equals("Lognormal"))
			cont = new Lognormal( parametros );
		else if(distO.equals("GammaDist") || distribucion.equals("Gamma"))
			cont = new Gamma( parametros );
		else if(distO.equals("WeibullDist") || distribucion.equals("Weibull"))
			cont = new Weibull( parametros );
		else if(distO.equals("ExponentialDist") || distribucion.equals("Exponential"))
			cont = new Exponential( parametros );
		else if(distO.equals("ErlangDist") || distribucion.equals("Erlang"))
			cont = new Erlang( parametros );
		else
			cont = new Uniform( parametros );
		return cont;
	}

	/**
	 * Performs a fitting procedure
	 * @param distribucion Name of the distribution
	 * @param data List of parameters
	 * @return The distribution
	 */
	public ContinuousDistribution ajustarParametros(String distribucion,
			DoubleArrayList data) {
		IDistribution cont;
		distribucion = distribucion.split(":")[0].trim();
		if(distribucion.equals("Normal") || distribucion.equals("NormalDist"))
			cont = new Normal( null );
		else if(distribucion.equals("Lognormal") || distribucion.equals("LognormalDist"))
			cont = new Lognormal( null );
		else if(distribucion.equals("Gamma") || distribucion.equals("GammaDist"))
			cont = new Gamma( null );
		else if(distribucion.equals("Weibull") || distribucion.equals("WeibullDist"))
			cont = new Weibull( null );
		else if(distribucion.equals("Exponential") || distribucion.equals("ExponentialDist"))
			cont = new Exponential( null );
		else if(distribucion.equals("Erlang") || distribucion.equals("ErlangDist"))
			cont = new Erlang( null );
		else
			cont = new Uniform( null );
		cont.ajustarParametros(data);
		return cont.darDistribucion();
	}
}

class Erlang extends EDistribution implements IDistribution{

	public Erlang(ArrayList<Double> paramatros){
		double k = 1;
		double lambda = 1;
		if (paramatros != null){
			k = paramatros.get(0);
			lambda = paramatros.get(1);
		}
		generadores = new ArrayList<String>();
		generadores.add("ErlangGen");
		
		parametros = new ArrayList<String>();
		parametros.add("k");
		parametros.add("lambda");
		
		distribucion = new ErlangDist((int) k, lambda);
	}

	public RandomVariateGen darGenerador(RandomStream s, String nombre,
			ArrayList<Double> paramatros) {
		double k = 1;
		double lambda = 1;
		if (paramatros != null){
			k = paramatros.get(0);
			lambda = paramatros.get(1);
		}
		if (nombre.equals("ErlangGen"))
			generador = new ErlangGen(s, (int) k, lambda);
		return generador;
	}

	public RandomVariateGen darGenerador(RandomStream s, String d) {
		
		generador = new ErlangGen(s, (ErlangDist) darDistribucion());
		return generador;
	}

	public void ajustarParametros(DoubleArrayList data) {
		double[] para = ErlangDist.getMLE(data.elements(), data.size());		
		distribucion = new ErlangDist( (int) para[0], para[1]);
	}
	
	public double[] getMoments(){
		double[] mom = new double[3];
		mom[0]  = distribucion.getMean();
		mom[1] = distribucion.getStandardDeviation();
		mom[2] = 2/Math.sqrt(((ErlangDist)distribucion).getK());
		return mom;
	}

	public String aString() {
		ErlangDist dis = (ErlangDist)distribucion;
		return "Erlang with k: " + nSig(dis.getK(), 3) +  
		", lambda:  " + nSig(dis.getLambda(), 3);
	}
	
	public String getName(){
		return "Erlang";
	}
}

class Exponential extends EDistribution implements IDistribution{

	public Exponential(ArrayList<Double> paramatros){
		double lambda = 1;
		if (paramatros != null){
			lambda = paramatros.get(0);
		}
		generadores = new ArrayList<String>();
		generadores.add("ExponentialGen");
		
		parametros = new ArrayList<String>();
		parametros.add("lambda");
		
		distribucion = new ExponentialDist(lambda);
	}

	public RandomVariateGen darGenerador(RandomStream s, String nombre,
			ArrayList<Double> paramatros) {
		double lambda = 1;
		if (paramatros != null){
			lambda = paramatros.get(0);
		}
		if (nombre.equals("ExponentialGen"))
			generador = new ExponentialGen(s, lambda);
		return generador;
	}

	public RandomVariateGen darGenerador(RandomStream s, String d) {
		
		generador = new ExponentialGen(s, (ExponentialDist) darDistribucion());
		return generador;
	}

	public void ajustarParametros(DoubleArrayList data) {
		double[] para = ExponentialDist.getMLE(data.elements(), data.size());		
		distribucion = new ExponentialDist( para[0]);
	}
	
	public double[] getMoments(){
		double[] mom = new double[3];
		mom[0]  = distribucion.getMean();
		mom[1] = distribucion.getStandardDeviation();
		mom[2] = 2.0;
		return mom;
	}

	public String aString() {
		ExponentialDist dis = (ExponentialDist)distribucion;
		return "Exponencial with lambda:  " + nSig(dis.getLambda(), 3);
	}
	
	public String getName(){
		return "Exponential";
	}
}

class Gamma extends EDistribution implements IDistribution{

	public Gamma(ArrayList<Double> paramatros){
		double alpha = 0.5;
		double lambda = 0.5;
		if (paramatros != null){
			alpha = paramatros.get(0);
			lambda = paramatros.get(1);
		}
		generadores = new ArrayList<String>();
		generadores.add("GammaGen");
		
		parametros = new ArrayList<String>();
		parametros.add("alpha");
		parametros.add("lambda");
		
		distribucion = new GammaDist(alpha, lambda);
	}

	public RandomVariateGen darGenerador(RandomStream s, String nombre,
			ArrayList<Double> paramatros) {
		double alpha = 0.5;
		double lambda = 0.5;
		if (paramatros != null){
			alpha = paramatros.get(0);
			lambda = paramatros.get(1);
		}
		if (nombre.equals("GammaGen"))
			generador = new GammaGen(s, alpha, lambda);
		return generador;
	}

	public RandomVariateGen darGenerador(RandomStream s, String d) {
		
		generador = new GammaGen(s, (GammaDist) darDistribucion());
		return generador;
	}

	public void ajustarParametros(DoubleArrayList data) {
		double[] para = GammaDist.getMLE(data.elements(), data.size());		
		distribucion = new GammaDist( para[0], para[1]);
	}
	
	public double[] getMoments(){
		double[] mom = new double[3];
		mom[0]  = distribucion.getMean();
		mom[1] = distribucion.getStandardDeviation();
		mom[2] = 2/Math.sqrt(((GammaDist)distribucion).getAlpha());
		return mom;
	}

	public String aString() {
		GammaDist dis = (GammaDist)distribucion;
		return "Gamma with alpha: " + nSig(dis.getAlpha(), 3) +  
		", lambda:  " + nSig(dis.getLambda(), 3);
	}
	
	public String getName(){
		return "Gamma";
	}
}

class Lognormal extends EDistribution implements IDistribution{

//	private static NormalGen normalGen;

	public Lognormal(ArrayList<Double> paramatros){
		double media = 0;
		double desv = 1;
		if (paramatros != null){
			media = paramatros.get(0);
			desv = paramatros.get(1);
		}
		generadores = new ArrayList<String>();
		generadores.add("LogNormalGen");

		parametros = new ArrayList<String>();
		parametros.add("mu");
		parametros.add("sigma");

		distribucion = new LognormalDist(media, desv);
	}

	public RandomVariateGen darGenerador(RandomStream s, String nombre,
			ArrayList<Double> paramatros) {
		double media = 30;
		double desv = 1;
		if (paramatros != null){
			media = paramatros.get(0);
			desv = paramatros.get(1);
		}
		if (nombre.equals("LogNormalGen"))
			generador = new LognormalGen(s, media, desv);
		return generador;
	}

	public RandomVariateGen darGenerador(RandomStream s, String d) {

		generador = new LognormalGen(s, (LognormalDist) darDistribucion());
		return generador;
	}

	public void ajustarParametros(DoubleArrayList data) {
		double[] para = LognormalDist.getMLE(data.elements(), data.size());		
		distribucion = new LognormalDist( para[0], para[1]);
	}

//	public static double nextDouble(LogNormal2Gen generador) {
//		if(normalGen == null){
//			normalGen = new NormalKindermannRamageGen(generador.getStream(),generador.getMu(),generador.getSigma());
//		}
//		System.out.println(normalGen + " mu:" + normalGen.getMu() + " sigma:" + normalGen.getSigma());
//		return Math.exp(normalGen.nextDouble());
//	}
	
	public double[] getMoments(){
		double[] mom = new double[3];
		mom[0]  = distribucion.getMean();
		mom[1] = distribucion.getStandardDeviation();
		mom[2] = (Math.exp(mom[1]*mom[1])+2)*(Math.sqrt(Math.exp(mom[1]*mom[1])-1));
		return mom;
	}

	public String aString() {
		LognormalDist dis = (LognormalDist)distribucion;
		return "Lognormal with mu: " + nSig(dis.getMu(), 3) +  
		", sigma:  " + nSig(dis.getSigma(), 3);
	}
	
	public String getName(){
		return "Lognormal";
	}
}

class Normal extends EDistribution implements IDistribution{

	public Normal(ArrayList<Double> paramatros){
		double media = 0;
		double desv = 1;
		if (paramatros != null){
			media = paramatros.get(0);
			desv = paramatros.get(1);
		}
		generadores = new ArrayList<String>();
		generadores.add("KindermannRamage(Default)");
		generadores.add("Polar");
		generadores.add("BoxMuller");
		generadores.add("NormalACR");
		
		parametros = new ArrayList<String>();
		parametros.add("Media");
		parametros.add("Desviacion");
		
		distribucion = new NormalDist(media, desv);
	}

	public RandomVariateGen darGenerador(RandomStream s, String nombre,
			ArrayList<Double> paramatros) {
		double media = 0;
		double desv = 1;
		if (paramatros != null){
			media = paramatros.get(0);
			desv = paramatros.get(1);
		}
		if (nombre.equals("Polar"))
			generador = new NormalPolarGen(s, media, desv);
		else if (nombre.equals("KindermannRamage(Default)"))
			generador = new NormalKindermannRamageGen(s, media, desv);
		else if (nombre.equals("BoxMuller"))
			generador = new NormalBoxMullerGen(s, media, desv);
		else 
			generador = new NormalACRGen(s, media, desv);
		return generador;
	}

	public RandomVariateGen darGenerador(RandomStream s, String d) {
		
		if(d == null || d.equals("KindermannRamage(Default)"))
			generador = new NormalKindermannRamageGen(s, (NormalDist) darDistribucion());
		else if(d.equals("Polar"))
			generador = new NormalPolarGen(s, (NormalDist) darDistribucion());
		else if(d.equals("BoxMuller"))
			generador = new NormalBoxMullerGen(s, (NormalDist) darDistribucion());
		else 
			generador = new NormalACRGen(s, (NormalDist) darDistribucion());
		return generador;
	}

	public void ajustarParametros(DoubleArrayList data) {
		double[] para = NormalDist.getMLE(data.elements(), data.size());		
		distribucion = new NormalDist( para[0], para[1]);
	}
	
	public double[] getMoments(){
		double[] mom = new double[3];
		mom[0]  = distribucion.getMean();
		mom[1] = distribucion.getStandardDeviation();
		mom[2] = 0.0;
		return mom;
	}

	public String aString() {
		NormalDist dis = (NormalDist)distribucion;
		return "Normal with mu: " + nSig(dis.getMean(), 3) +  
		", sigma:  " + nSig(dis.getSigma(), 3);
	}
	
	public String getName(){
		return "Normal";
	}
}

class Uniform extends EDistribution implements IDistribution{

	public Uniform(ArrayList<Double> paramatros){
		double min = 0;
		double max = 1;
		if (paramatros != null){
			min = paramatros.get(0);
			max = paramatros.get(1);
		}
		generadores = new ArrayList<String>();
		generadores.add("UniformGen");
		
		parametros = new ArrayList<String>();
		parametros.add("min");
		parametros.add("max");
		
		distribucion = new UniformDist(min, max);
	}

	public RandomVariateGen darGenerador(RandomStream s, String nombre,
			ArrayList<Double> paramatros) {
		double min = 0;
		double max = 1;
		if (paramatros != null){
			min = paramatros.get(0);
			max = paramatros.get(1);
		}
		if (nombre.equals("UniformGen"))
			generador = new UniformGen(s, min, max);
		return generador;
	}

	public RandomVariateGen darGenerador(RandomStream s, String d) {		
		generador = new UniformGen(s, (UniformDist) darDistribucion());
		return generador;
	}

	
	public void ajustarParametros(DoubleArrayList data) {
		double[] para = UniformDist.getMLE(data.elements(), data.size());		
		distribucion = new UniformDist( para[0], para[1]);
	}
	
	public double[] getMoments(){
		double[] mom = new double[3];
		mom[0]  = distribucion.getMean();
		mom[1] = distribucion.getStandardDeviation();
		mom[2] = 0.0;
		return mom;
	}

	public String aString() {
		UniformDist dis = (UniformDist)distribucion;
		return "Uniforme with min: " + nSig(dis.getA(), 3) +  
		", max:  " + nSig(dis.getB(), 3);
	}
	
	public String getName(){
		return "Uniform";
	}
}

class Weibull extends EDistribution implements IDistribution{

	public Weibull(ArrayList<Double> paramatros){
		double alpha = .5;
		double lambda = 1;
		double delta = 0;
		if (paramatros != null){
			alpha = paramatros.get(0);
			lambda = paramatros.get(1);
			delta = paramatros.get(2);
		}
		generadores = new ArrayList<String>();
		generadores.add("WeibulGen");
		
		parametros = new ArrayList<String>();
		parametros.add("alpha");
		parametros.add("lambda");
		parametros.add("delta");
		
		distribucion = new WeibullDist(alpha, lambda, delta);
	}

	
	public RandomVariateGen darGenerador(RandomStream s, String nombre,
			ArrayList<Double> paramatros) {
		double alpha = 0.5;
		double lambda = 1;
		double delta = 0;
		if (paramatros != null){
			alpha = paramatros.get(0);
			lambda = paramatros.get(1);
			delta = paramatros.get(2);
		}
		if (nombre.equals("WeibullGen"))
			generador = new WeibullGen(s, alpha, lambda, delta);
		return generador;
	}

	public RandomVariateGen darGenerador(RandomStream s, String d) {
		
		generador = new WeibullGen(s, (WeibullDist) darDistribucion());
		return generador;
	}

	public void ajustarParametros(DoubleArrayList data) {
		double[] para = WeibullDist.getMLE(data.elements(), data.size());		
		distribucion = new WeibullDist( para[0], para[1], para[2]);
	}
	
	public double[] getMoments(){
		double[] mom = new double[3];
		mom[0]  = distribucion.getMean();
		mom[1] = distribucion.getStandardDeviation();
		mom[2] = 1.0;
		return mom;
	}

	public String aString() {
		WeibullDist dis = (WeibullDist)distribucion;
		return "Weibull with alpha: " + nSig(dis.getAlpha(), 3) +  
		", lambda:  " + nSig(dis.getLambda(), 3) + ", delta: " + nSig(dis.getDelta(), 3);
	}
	
	public String getName(){
		return "Weibull";
	}
}
