package jphase.fit;

import jphase.DenseDiscPhaseVar;
import jphase.DiscPhaseVar;


/**
 * This class implements the Matching Moments method 
 * proposed by Telek and Heindl in "Matching Moments for Acyclic
 * discrete and continuous Phase-Type distributions of 
 * Second order", 2002. This method returns a 
 * Discrete Phase-type distribution. 
 * @author Juan F. Perez
 * @version 1.0 
 */
public class MomentsADPH2Fit extends MomentsDiscPhaseFitter{

	/**
     * @see jphase.fit.MomentsDiscPhaseFitter#MomentsDiscPhaseFitter(int[])
     */
	public MomentsADPH2Fit(int[] data){
		super(data);
	}

	/**
     * @see jphase.fit.MomentsDiscPhaseFitter#MomentsDiscPhaseFitter(double, double, double)
     */
	public MomentsADPH2Fit(double m1, double m2, double m3){
		super(m1, m2, m3);
	}

	
	/**
	 * Precision for calculations and convergence criterion 
	 */
	private static double precision = 10E-6; 
	
	
	@Override
	public DiscPhaseVar fit() {		
		double[] param = getParam(m1, m2, m3);
		double[] vector = new double[2];
		vector[0]=param[0];
		vector[1]=1-param[0];
		double[][] matriz = new double[2][2];
		matriz[0][0]=1-param[1];
		matriz[0][1]=param[1];
		matriz[1][1]=1-param[2];
		return new DenseDiscPhaseVar(vector, matriz);
	}
	
	
	/**
	 * This method verifies the bounds of the moments to be matched,
	 * and, if they are not feasible, it corrects them if possible.
	 * After the correction, it calls the method to return the 
	 * parameters of the Acyclic Discrete Phase Distribution 
	 * of order 2: param[0]=p, param[1]=beta1, param[2]=beta2  
	 * @param f1 First factorial moment to be matched
	 * @param f2 Second factorial moment to be matched
	 * @param f3 Third factorial moment to be matched
	 * @return Parameters of the Acyclic Discrete Phase Distribution
	 * of order 2. param[0]=p, param[1]=beta1, param[2]=beta2  
	 */
	private double[] getParam(double f1, double f2, double f3){
		double[] param = new double[3];
		double d = 2*f1*f1-2*f1-f2;
		d = Math.sqrt(2*d);
		double g  = 6/(Math.pow(2*f1+d,3))*( f1*(2*f1+d)*(3*f2+2*f1)*(f2-2*f1+2)-2*f2*f2*(f2-d) );
		double f12 = f1*f1;//f1^2
		if(f1>=1){
			if(f1<2){
				if(f2>=2*f1-2){
					if(f2<2*f12-2*f1){
						if(f3 < g ) {
							System.out.println("Third Factorial Moment non Feasible (too small) "+f3+" < g("+g+")");
							f3 = g;
							System.out.println("Third Factorial Moment fixed on "+f3+" for feasibility");
						}else if(f3 > 3*f2*(f2-2*f1+2)/(2*f1-2)){
							System.out.println("Third Factorial Moment non Feasible (too large) "+f3+" > 3*f2*(f2-2*f1+2)/(2*f1-2) ("+(3*f2*(f2-2*f1+2)/(2*f1-2))+")");
							f3 = 3*f2*(f2-2*f1+2)/(2*f1-2);
							System.out.println("Third Factorial Moment fixed on "+f3+" for feasibility");
						}
					}else{
						if(f3 < 3*f2*(f2-2*f1+2)/(2*f1-2)){
							System.out.println("Third Factorial Moment non Feasible (too small) "+f3+" < 3*f2*(f2-2*f1+2)/(2*f1-2) ("+(3*f2*(f2-2*f1+2)/(2*f1-2))+")");
							f3 = 3*f2*(f2-2*f1+2)/(2*f1-2);
							System.out.println("Third Factorial Moment fixed on "+f3+" for feasibility");
						}	
					}
					getParamFactible(f1,f2,f3, param);
				}else{
					System.out.println("Second factorial moment non feasible "+f2+" < 2*f1-2 ("+(2*f1-2)+")");
				}
			}else{
				if(f2>=1.5*f12-2*f1){
					if(f2<2*f1-2){
						if(f3 < g ) {
							System.out.println("Third Factorial Moment non Feasible (too small) "+f3+" < g("+g+")");
							f3 = g;
							System.out.println("Third Factorial Moment fixed on "+f3+" for feasibility");
						}else if(f3 > 6*(f1-1)*(f2+f1-f12) ){
							System.out.println("Third Factorial Moment non Feasible (too large) "+f3+" > 6*(f1-1)*(f2+f1-f12) ("+(6*(f1-1)*(f2+f1-f12))+")");
							f3 = 6*(f1-1)*(f2+f1-f12);
							System.out.println("Third Factorial Moment fixed on "+f3+" for feasibility");
						}
					}else if(f2<2*f12-2*f1){
						if(f3 < g ) {
							System.out.println("Third Factorial Moment non Feasible (too small) "+f3+" < g("+g+")");
							f3 = g;
							System.out.println("Third Factorial Moment fixed on "+f3+" for feasibility");
						}else if(f3 > 3*f2*(f2-2*f1+2)/(2*f1-2) ){
							System.out.println("Third Factorial Moment non Feasible (too large) "+f3+" > 3*f2*(f2-2*f1+2)/(2*f1-2) ("+(3*f2*(f2-2*f1+2)/(2*f1-2))+")");
							f3 = 3*f2*(f2-2*f1+2)/(2*f1-2);
							System.out.println("Third Factorial Moment fixed on "+f3+" for feasibility");
						}
					}else{
						if(f3 < 3*f2*(f2-2*f1+2)/(2*f1-2) ){
							System.out.println("Third Factorial Moment non Feasible (too small) "+f3+" < 3*f2*(f2-2*f1+2)/(2*f1-2) ("+(3*f2*(f2-2*f1+2)/(2*f1-2))+")");
							f3 = 3*f2*(f2-2*f1+2)/(2*f1-2);
							System.out.println("Third Factorial Moment fixed on "+f3+" for feasibility");
						}
					}
					getParamFactible(f1,f2,f3, param);
				}else{
					System.out.println("Second factorial moment non feasible "+f2+" < 1.5*f1^2-2*f1 ("+(1.5*f12-2*f1)+")");
				}
			}
		}else{
			System.out.println("First factorial moment non feasible "+f1+" < 1");
		}
		return param;
	}

	/**
	 * Solve the equation system to get the parameters of the distribution, 
	 * given that the moments are feasible.  This are stored in the vector
	 * param: param[0]=p, param[1]=beta1, param[2]=beta2  
	 * @param f1 First moment to be matched
	 * @param f2 Second moment to be matched
	 * @param f3 Third moment to be matched
	 * @param param Vector to store the parameters of the 
	 * Acyclic Discrete Phase Distribution of order 2: param[0]=p, 
	 * param[1]=beta1, param[2]=beta2  
	 */
	private void getParamFactible(double f1, double f2, double f3, double[] param) {
		double d = 2*f1*f1 - 2*f1 - f2;
		double c = 3*f2*f2-2*f1*f3;System.out.println("c: "+c);
		double b = 3*f1*f2 - 6*(f1+f2-f1*f1)-f3;
		double a = b*b-6*c*d;
		if(Math.abs(c)>precision){
			if(c > 0){
				// Case c > 0
				double sqa=Math.sqrt(a);
				param[0]=(-b+6*f1*d+sqa)/( b+sqa );
				param[1]=(b-sqa)/c;
				param[2]=(b+sqa)/c;
			}else{
				// Case c < 0
				double sqa=Math.sqrt(a);
				param[0]=(b-6*f1*d+sqa)/( -b+sqa );
				param[1]=(b+sqa)/c;
				param[2]=(b-sqa)/c;
			}
		}else{
			System.out.println("Adjusted to geometric");
			param[0]=0;
			param[2]=1/f1;			
		}
	}
}