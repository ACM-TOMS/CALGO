package jphase.fit;

import jphase.DenseContPhaseVar;
import jphase.ContPhaseVar;


/**
 * This class implements the moment-matching method 
 * proposed by Telek and Heindl in "Matching Moments for Acyclic
 * discrete and continuous Phase-Type distributions of 
 * Second order", 2002. The method returns a 
 * continuous Phase-type variable.
 * @author Juan F. Perez
 * @version 1.0  
 */
public class MomentsACPH2Fit extends MomentsContPhaseFitter{

	/**
	 *@see jphase.fit.MomentsContPhaseFitter#MomentsContPhaseFitter(double[]) 
	 */
	public MomentsACPH2Fit(double[] data){
		super(data);
	}
	
	/**
	 * This constructor does not allow the determination
     * of any parameter of the fit process
	 * @see jphase.fit.MomentsContPhaseFitter#MomentsContPhaseFitter(double,double,double)
	 */
	public MomentsACPH2Fit(double m1, double m2, double m3){
		super(m1, m2, m3);
	}
	
	
	/**
     * TODO get precision
	 * precision for calculations and convergence criterion 
	 */
	public static double precision = 10E-6; 
	
    
    
	@Override
	public ContPhaseVar fit(){
		double[] param = getParam(m1,m2,m3);
		double[] vector = new double[2];
		vector[0]=param[0];
		vector[1]=1-param[0];
		double[][] matriz = new double[2][2];
		matriz[0][0]=-param[1];
		matriz[0][1]=param[1];
		matriz[1][1]=-param[2];
		return new DenseContPhaseVar(vector, matriz);
	}
		
	
	/**
	 * This method verifies the bounds of the moments to be matched,
	 * and, if they are not feasible, it corrects them if possible.
	 * After the correction, it calls the method to return the 
	 * parameters of the Acyclic Continuous Phase Distribution
	 * of order 2: param[0]=p, param[1]=lambda1, param[2]=lambda[2] 
	 * @param m1 First moment to be matched
	 * @param m2 Second moment to be matched
	 * @param m3 Third moment to be matched
	 * @return Parameters of the Acyclic Continuous Phase Distribution
	 * of order 2: param[0]=p, param[1]=lambda1, param[2]=lambda[2]  
	 */
	private double[] getParam(double m1, double m2, double m3){
		double[] param = new double[3];
		double cx2 = m2/(m1*m1)-1; //squared coefficient of variation
		double m13 = m1*m1*m1;//m1^3
		if(m1>0){
			if(cx2>=0.5){
				if(cx2<=1){
					System.out.println("HypoExponential");
					if(m3 < 3*m13*(3*cx2-1+Math.sqrt(2)*Math.pow(1-cx2,1.5)) ) {
						System.out.println("Third Moment non Feasible (too small) "+m3+" (0.5<=Cx^2<=1");
						m3=3*m13*(3*cx2-1+Math.sqrt(2)*Math.pow(1-cx2,1.5));
						System.out.println("Third Moment fixed on "+m3+" for feasibility");
					}else if(m3 > 6*m13*cx2){
						System.out.println("Third Moment non Feasible (too large) "+m3+" (0.5<=Cx^2<=1");
						m3=6*m13*cx2;
						System.out.println("Third Moment fixed on "+m3+" for feasibility");
					}
					getParamFactible(m1,m2,m3, param);
				}else{
					System.out.println("HyperExponencial");
					if(m3 < 1.5*m13*(1+cx2)*(1+cx2) ){
						System.out.println("Third Moment non Feasible (too small) "+m3+" (Cx^2>=1");
						m3 = 1.5*m13*(1+cx2)*(1+cx2);
						System.out.println("Third Moment fixed on "+m3+" for feasibility");
					}	
					getParamFactible(m1,m2,m3, param);
				}
			}else{
				System.out.println("Squared Coefficient of Variance non Feasible "+cx2+" <0.5");
			}
		}else{
			System.out.println("First moment non feasible "+m1+" < 0");
		}
		
		return param;
	}

	/**
	 * Solve the equation system to get the parameters of the distribution, 
	 * given that the moments are feasible. This are stored in the vector
	 * param: param[0]=p, param[1]=lambda1, param[2]=lambda2   
	 * @param m1 First moment to be matched
	 * @param m2 Second moment to be matched
	 * @param m3 Third moment to be matched
	 * @param param Vector to store the parameters of the Acyclic 
	 * Continuous Phase Distribution of order 2: param[0]=p, 
	 * param[1]=lambda1, param[2]=lambda2  
	 */
	private void getParamFactible(double m1, double m2, double m3, double[] param) {
		double d = 2*m1*m1 - m2;
		double c = 3*m2*m2-2*m1*m3;
		double b = 3*m1*m2-m3;
		double a = b*b-6*c*d;
		if(Math.abs(c)>precision){
			if(c > 0){
				System.out.println("c > 0");
				double sqa=Math.sqrt(a);
				param[0]=(-b+6*m1*d+sqa)/( b+sqa );
				param[1]=(b-sqa)/c;
				param[2]=(b+sqa)/c;
			}else{
				System.out.println("c < 0");
				double sqa=Math.sqrt(a);
				param[0]=(b-6*m1*d+sqa)/( -b+sqa );
				param[1]=(b+sqa)/c;
				param[2]=(b-sqa)/c;
			}
		}else{
			System.out.println("Adjusted to exponential");
			param[0]=0;
			param[2]=1/m1;			
		}
	}

}
