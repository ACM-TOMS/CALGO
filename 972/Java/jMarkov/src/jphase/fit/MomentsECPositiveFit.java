package jphase.fit;

import jphase.DenseContPhaseVar;
import jphase.ContPhaseVar;


/**
 * This class implements the Matching Moments method 
 * proposed by Osogami and Harchol in "Closed form 
 * solutions for mapping general distributions to 
 * quasi-minimal PH distributions", 2005. 
 * This class implements the Positive solution. 
 * This method matches the first three moments. 
 * This method returns an Erlang-Coxian distributions, 
 * a subclass of Phase-Type distributions.
 * @author Juan F. Perez
 * @version 1.0 
 */
public class MomentsECPositiveFit extends MomentsContPhaseFitter{

	/**
	 * @see jphase.fit.MomentsContPhaseFitter#MomentsContPhaseFitter(double[])
	 */
	public MomentsECPositiveFit(double[] data){
		super(data);
	}
	
	/**
     * @see jphase.fit.MomentsContPhaseFitter#MomentsContPhaseFitter(double,double,double)
	 */
	public MomentsECPositiveFit(double m1, double m2, double m3){
		super(m1, m2, m3);
	}
	
	/**
	 * precision for the calculations 
	 */
	private static double precision = 10E-6; 
	
	/**
	 * Fits a Phase Type distribution with no mass at zero 
	 * from a set of moments with the
	 * Positive method described by Osogami et al.  
	 * @return Phase variable found
	 */
	@Override
	public ContPhaseVar fit() {
		System.out.println("m1: "+m1);
		System.out.println("m2: "+m2);
		System.out.println("m3: "+m3);
		try{
			double[] param = getParam();
			DenseContPhaseVar EC = DenseContPhaseVar.ErlangCoxian((int)param[0],param[1],param[2],param[3],param[4],param[5]);
			int n = param.length;
			if(n == 7){
				return EC.sum(DenseContPhaseVar.expo(param[6]), new DenseContPhaseVar((int)param[0]+1));
			}else{
				return EC.mix(param[6], DenseContPhaseVar.expo(param[7]),  new DenseContPhaseVar((int)param[0]+1));
			}
			
		}catch (IllegalArgumentException e){
			return null;
		}
	}
	
	/**
	 * Determines if the set of normalized non-centered 
	 * moments {m1, m2, m3} are in the set of 
	 * distributions U, as defined by Osogami and Harchol.
	 * @param m2 Second normalized non-centered moment
	 * @param m3 Third normalized non-centered moment
	 * @return True, if the set of normalized non-centered 
	 * moments {m1, m2, m3} are in the set of 
	 * distributions U, as defined by Osogami and Harchol. 
	 * False, otherwise. 
	 */
	private boolean inU(double m2, double m3) {
		if(m3 <= 2*m2 + 1 || m2 > 2 || m2 < 1.0  )return false;
		else{
			int i = 1;
			while(i < 1000){
				double limInf = (i+2.0)/(i+1.0);
				double limSup = (i+1.0)/(i);
				if(limInf < m2 && m2 < limSup)return true;
				i++;
			}
		}
		return false;
	}


	/**
	 * Solves the equation system to get the parameters of the distribution, 
	 * if the moments are feasible  
	 * @return Parameters of the Acyclic Continuous Phase Distribution
	 * of order n. param[0]=n, param[1]=p, param[2]=lambdaY,
	 * param[3]=lambdaX1, param[4]=lambdaX2, param[5]=pX. If it is
	 * necessary, a new exponential phase with parameter lambda 
	 * will be aggregated through convolution: param[6]=lambda.
	 * If it is necessary, a new exponential phase with parameter lambda
	 * will be mixed up with probability 1-pMix: param[6]=lambda, 
	 * param[7] = pMix.  
	 * @throws IllegalArgumentException    
	 */
	public double[] getParam() throws IllegalArgumentException{
		double n2 = m2/(m1*m1);
		double n3 = m3/(m1*m2);
		double r = n3/n2;
		MomentsECCompleteFit complete;
		if(n2 >1 && n3 > n2){
			if( inU(n2,n3) || n3 == 2*n2-1 || (r != 1.5 && n3 < 2*n2 - 1) ){
				if( r < 1.5  && r < 1.5){
				
					double k = FitterUtils.floor((2*n2-n3)/(n3-n2),precision);
					if(n3 >= ((k+1)*n2 + (k+4))/(2*(k+2)) * n2 ){
						System.out.println("m3 >= ((k+1)*m2 + (k+4))/(2*(k+2)) * m2");
						double[] param = new double[8];
						double[] paramComplete = new double[6];
						double w = (2-n2)/(4*(1.5-r));
						double a = (2-n2)*(2-n2);
						double pMix = a/(a + 4*(2*n2-1-n3));
						double m1X = m1/(pMix+(1-pMix)*w);
						double m2X = 2*w;
						double m3X = 2*m2X - 1;
						double lambda = 1/(w*m1X);
						complete = new MomentsECCompleteFit(m1X, m2X*m1X*m1X, m3X*m1X*m2X);
						paramComplete = complete.getParam();
						System.arraycopy(paramComplete, 0, param, 0, 6);
						param[6]=lambda;
						param[7]=pMix;
						return param;
						
					}else if(n3 < ((k+1)*n2 + (k+4))/(2*(k+2)) * n2 && Math.abs(n2 - 2) <= precision){
						System.out.println("m3 < ((k+1)*m2 + (k+4))/(2*(k+2))*m2  AND m2 = 2");
						double[] param = new double[7];
						double[] paramComplete = new double[6];
						double z = (n3 - 2*(k+3)/(k+2))/(3-n3);
						double m1X = m1/(1+z);
						double m2X = 2*(1+z);
						double m3X = (k+3)/(k+2)*m2X;
						double lambda = 1/(z*m1X);
						complete = new MomentsECCompleteFit(m1X, m2X*m1X*m1X, m3X*m1X*m2X);
						paramComplete = complete.getParam();
						System.arraycopy(paramComplete, 0, param, 0, 6);
						param[6]=lambda;
						return param;
					}else{
						System.out.println("m3 < ((k+1)*m2 + (k+4))/(2*(k+2))*m2  AND m2 != 2");
						double[] param = new double[7];
						double[] paramComplete = new double[6];
						double a = 2*(k+3)/(k+2)*(n2-2);
						double z = ( n2*(n3-3-a) + n2*FitterUtils.sqrt((n3-3)*(n3-3) + 4*a*(1.5-r) ,precision) )/( a*(n2-2) );
						double m1X = m1/(1+z);
						double m2X = (1+z)*(n2*(1+z)-2*z);
						double m3X = (k+3)/(k+2)*m2X;
						double lambda = 1/(z*m1X);
						complete = new MomentsECCompleteFit(m1X, m2X*m1X*m1X, m3X*m1X*m2X);
						paramComplete = complete.getParam();
						System.arraycopy(paramComplete, 0, param, 0, 6);
						param[6]=lambda;
						return param;
					}
				
				}else{
					double[] param = new double[6];
					complete = new MomentsECCompleteFit(m1,n2,n3);
					param = complete.getParam();
					return param;
				}	
									
			}else{
				throw new IllegalArgumentException("The set of moments is not in the set  " +
						" U || m3 == 2*m2-1 || (r != 1.5 && m3 < 2*m2 - 1)");
			}
		}else{
			throw new IllegalArgumentException("The set of moments is not representable by the PH type distributions\n " +
			"It is not in m2 > 1 AND m3 > m2");
		}	
	}	
}