package jphase.fit;

import no.uib.cipr.matrix.DenseVector;
import jphase.DenseContPhaseVar;
import jphase.ContPhaseVar;


/**
 * This class implements the Matching Moments method 
 * proposed by Osogami and Harchol in "Closed form 
 * solutions for mapping general distributions to 
 * quasi-minimal PH distributions", 2005. 
 * This class implements the Complete solution. 
 * This method matches the first three moments. 
 * This method returns an Erlang-Coxian distributions, 
 * a subclass of Phase-Type distributions. 
 * @author Juan F. Perez
 * @version 1.0 
 */
public class MomentsECCompleteFit extends MomentsContPhaseFitter{

	/**
	 * @see jphase.fit.MomentsContPhaseFitter#MomentsContPhaseFitter(double[])
	 */
	public MomentsECCompleteFit(double[] data){
		super(data);
	}

	/**
	 * @see jphase.fit.MomentsContPhaseFitter#MomentsContPhaseFitter(double,double,double)
	 */
	public MomentsECCompleteFit(double m1, double m2, double m3){
		super(m1, m2, m3);
	}
	
	/**
	 * precision for calculations
	 */
	private static double precision = 10E-6; 
		
	/**
	 * Fits a Phase Type distribution from a set of moments with the
	 * Complete method described by Osogami et al.  
	 * @return Phase variable found
	 */
	@Override
	public ContPhaseVar fit() {
		System.out.println("m1: "+m1);
		System.out.println("m2: "+m2);
		System.out.println("m3: "+m3);
		double[] param = getParam();
		if(param[4]==0 && param[5]==0){
			DenseContPhaseVar var = (DenseContPhaseVar)DenseContPhaseVar.Erlang(param[2],(int)param[0]-2).sum(
					DenseContPhaseVar.expo(param[3]), new DenseContPhaseVar((int)param[0]-1));
			double[] initCond = new double[(int)param[0]-1];
			initCond[0]=param[1];
			var.setVector(new DenseVector(initCond));
			return var;
			
		}else{
			return DenseContPhaseVar.ErlangCoxian((int)param[0],param[1],param[2],param[3],param[4],param[5]);
		}
		
	}
	
	
	/**
	 * Solves the equation system to get the parameters of the distribution, 
	 * if the moments are feasible  
	 * @return Parameters of the Acyclic Continuous Phase Distribution
	 * of order n. param[0]=n, param[1]=p, param[2]=lambdaY,
	 * param[3]=lambdaX1, param[4]=lambdaX2, param[5]=pX    
	 */
	public double[] getParam(){
		double[] param = new double[6];
		double n2 = m2/(m1*m1);
		double n3 = m3/(m1*m2);
		
	if(n2 <= 1 || n3 <=n2){
			System.out.println("Second or Third Moment non Feasible (too little) "+n3+" (0.5<=Cx^2<=1");
			return param;
		}else{
			//Determination of p: param[1]
			if(n3 > 2*n2 - 1 &&  Math.abs((int)(1/(n2-1)) - 1.0/(n2-1)) <= precision){
				System.out.println("m3 > 2*m2 -1 and 1/(m2-1) isInteger");
				param[1] = ( n2*n2 + 2*n2 - 1 )/( 2*n2*n2 );
			}else if(n3 < 2*n2 -1){
				System.out.println("m3 < 2*m2 -1");
				param[1] = 1.0/( 2*n2 - n3 );
			}else{
				System.out.println("m3 = 2*m2 -1 or m3 > 2*m2 -1 and 1/(m2-1) is NOT Integer");
				param[1] = 1.0;
			}
			//W distribution moments
			double m1W = m1/param[1];
			double m2W = param[1]*n2;
			double m3W = param[1]*n3;
			
			//determination of n: param[0]
			if(m3W ==2*m2W-1 && m2W <= 2){
				System.out.println("m3W ==2*m2W-1 and m2W <= 2");
				param[0] = FitterUtils.ceil(m2W/(m2W-1), precision) - 1.0;
			}else{
				System.out.println("m3W != 2*m2W-1 or m2W > 2");
				param[0] = FitterUtils.floor( m2W/(m2W-1)  + 1.0, precision);
			}
			
			//Exponential distribution case
			if(param[0] == 1.0){
				param[2] = 0.0;//irrelevant
				
				param[3] = 1.0 / m1W;
				param[4] = 0.0;//irrelevant
				param[5] = 0.0;
				return param;
			}
			
			//X: distribution moments
			int n = (int)param[0];
			double m2X = ( (n-3)*m2W - (n-2) )/((n-2)*m2W - (n-1)  ) ;
			double m1X = m1W/( (n-2)*m2X -(n-3) );
			
			double alpha = (n-2)*(m2X-1)*( n*(n-1)*m2X*m2X - n*(2*n-5)*m2X + (n-1)*(n-3) );
			double beta = ( (n-1)*m2X-(n-2) )*( (n-2)*m2X - (n-3) )*( (n-2)*m2X - (n-3) );
			double m3X = ( beta*m3W - alpha )/m2X;
			
			//determination of u and v
			double u = 0.0;
			double v = 0.0;
			if(Math.abs(3*m2X - 2*m3X)<0.01){ 
				System.out.println("3*m2X == 2*m3X");
				u = 1.0;
				v = 0.0;
			}else{
				System.out.println("3*m2X != 2*m3X");
				u = (6 - 2*m3X)/(3*m2X - 2*m3X);
				v = (12 - 6*m2X)/(m2X*(3*m2X - 2*m3X));
			}
			
			//determination of parameters 
			param[2]=  1.0/( m1X*(m2X-1) );//lambdaY
			param[3]= (u + FitterUtils.sqrt(u*u - 4*v, precision))/(2*m1X);//lambdaX1
			param[4]= (u - FitterUtils.sqrt(u*u - 4*v, precision))/(2*m1X);//lambdaX2
			param[5]= param[4]*(param[3]*m1X - 1)/param[3] ;//pX
			
			return param;
			
		}
	}

}
