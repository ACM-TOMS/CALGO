package jphase.fit;

import jphase.DenseContPhaseVar;
import jphase.ContPhaseVar;
import no.uib.cipr.matrix.DenseVector;


/**
 * This class implements the moment-matching method 
 * proposed by Bobbio, Horvath and Telek in "Matching
 * threee moments with minimal acyclic Phase-Type 
 * distributions", 2005. The method matches the first 
 * three moments with a subclass of Phase-type 
 * distributions known as Acyclic Phase-type distributions.  
 * @author Juan F. Perez
 * @version 1.0 
 */
public class MomentsACPHFit extends MomentsContPhaseFitter{

	/**
	 * @see jphase.fit.MomentsContPhaseFitter#MomentsContPhaseFitter(double[])
	 */
	public MomentsACPHFit(double[] data){
		super(data);
	}
	
	/**
	 * This constructor does not allow the determination
	 * of any parameter of the fit process
	 * @see jphase.fit.MomentsContPhaseFitter#MomentsContPhaseFitter(double,double,double)
	 */
	public MomentsACPHFit(double m1, double m2, double m3){
		super(m1, m2, m3);
	}
		
	/**
	 * Solves the equation system to get the 
	 * parameters of the distribution, 
	 * if the moments are feasible  
	 * @return Parameters of the Acyclic Continuous Phase Distribution
	 * of order n. param[0]=n, param[1]=p, param[2]=lambdaY,
	 * param[3]=lambdaX1, param[4]=lambdaX2, param[5]=pX    
	 */
	@Override
	public ContPhaseVar fit(){
		double n2 = m2/(m1*m1);
		double n3 = m3/(m1*m2);
        ContPhaseVar res = new DenseContPhaseVar();

        int n = getSize(n2, n3);
		if(n==0)return null;
		else if(n == 1){
		    double[] alpha = {2/n2};
            double[][] A ={{-1/m1}};
            res = new DenseContPhaseVar(alpha, A);
        }else{
            
			if(n2 <= n/(n-1.0) || n3 <= 2*n2-1.0){
                double c = 12*n2*n2*(n+1.0) + 16*n3*(n+1.0) + n2*(n*(n3-15)*(n3+1) - 8*(n3+3));
				double b = 2*(4 - n*(3*n2-4.0)) /( n2*(4.0+n-n*n3) + Math.sqrt(n*n2*c)  )   ; 
				double a = (b*n2 - 2)*(n-1.0)*b/((b-1.0)*n);
				double p = (b-1.0)/a;
				
				double lambda = (1.0 + a*p)/m1;
				double mu = lambda*(n-1.0)/a;
                
				ContPhaseVar temp = DenseContPhaseVar.Erlang(mu,n-1);
                double[] initCond = new double[n-1];
                initCond[0]=p;
                temp.setVector(new DenseVector(initCond));
                res = temp.sum(DenseContPhaseVar.expo(lambda),new DenseContPhaseVar(n));
			}else if(n2 > n/(n-1.0)){
				double f = calcF(n, n2, n3);
                double a = 2*(f-1.0)*(n-1.0)/((n-1.0)*(n2*f*f-2*f+2)-n);
                double p = (f-1.0)*a;
                
                double lambda = (p+a)/m1;
                double mu = lambda*(n-1)/a;
                
                ContPhaseVar temp = DenseContPhaseVar.expo(lambda);
                double[] initCond = new double[1];
                initCond[0]=p;
                temp.setVector(new DenseVector(initCond));
                res = temp.sum(DenseContPhaseVar.Erlang(mu,n-1),new DenseContPhaseVar(n));
			}else{
                System.out.println("Error when computing resulting distribution: moments cannot be matched.");
			}
		}
		
        return res;
	}
	
	/**
     * @param n number of phases
     * @param n2 second normalized moment
     * @param n3 third normalized moment
     * @return f value
     */
    private double calcF(int n, double n2, double n3) {
        double K1 = n - 1;
        double K2 = n - 2;
        double K3 = 3*n2 - 2*n3;
        double K4 = n3 - 3;
        double K5 = n - n2;
        double K6 = 1 + n2 - n3;
        double K7 = n + n2 - n*n2;
        double K8 = 3 + 3*n2*n2 + n3 - 3*n2*n3;
        double K9 = 108*K1*K1*(4*K2*K2*K3*n*n*n2 + K1*K1*K2*K4*K4*n*n2*n2 + 
                4*K1*K5*(K5*K5 - 3*K2*K6*n*n2)+Math.sqrt(-16*K1*K1*Math.pow(K7,6)
                        +Math.pow( 4*K1*K5*K5*K5 + K1*K1*K2*K4*K4*n*n2*n2 + 
                                4*K2*n*n2*(K4*n*n-3*K6*n2+K8*n) ,2) ));
        double K10 = K4*K4/(4*K3*K3) - K5/(K1*K3*n2);
        double K11 = Math.pow(2, 1.0/3.0)*(3*K5*K5 + K2*(K3+2*K4)*n*n2)/(K3 * Math.pow(K9, 1.0/3.0) * n2);
        
        double K12 = Math.pow(K9, 1.0/3.0)/(3*Math.pow(2,7.0/3.0)*K1*K1*K3*n2);
        double K13 = Math.sqrt(K10 + K11 + K12);
        double K14 = (6*K1*K3*K4*K5 + 4*K2*K3*K3*n - K1*K1*K4*K4*K4*n2)/
                        (4*K1*K1*K3*K3*K3*K13*n2);
        double K15 = -K4/(2*K3);
        double K16 = Math.sqrt(2*K10 - K11 - K12 - K14);
        double K17 = Math.sqrt(2*K10 - K11 - K12 + K14);
        double K18 = 36*K5*K5*K5 + 36*K2*K4*K5*n*n2 + 9*K1*K2*K4*K4*n*n2*n2
                        - Math.sqrt(81*(4*K5*K5*K5+4*K2*K4*K5*n*n2+K1*K2*K4*K4*n*n2*n2)*
                                (4*K5*K5*K5+4*K2*K4*K5*n*n2+K1*K2*K4*K4*n*n2*n2)
                                -48*Math.pow(3*K5*K5+2*K2*K4*n*n2 ,3));
        double K19 = -K5/(K1*K4*n2) - Math.pow(2, 2/3)*(3*K5*K5 + 2*K2*K4*n*n2)/
                           (Math.pow(3*K18, 1.0/3.0)*K1*K4*n2) - Math.pow(K18, 1.0/3.0)/
                           (Math.pow(6, 2.0/3.0)*K1*K4*n2);
        double K20 = 6*K1*K3*K4*K5 + 4*K2*K3*K3*n - K1*K1*K4*K4*K4*n2;
        double K21 = K11 + K12 + K5/(2*n*K1*K3);
        double K22 = Math.sqrt(3*K4*K4/(4*K3*K3) - 3*K5/(K1*K3*n2) + 
                Math.sqrt(4*K21*K21 - n*K2/(n2*K1*K1*K3)));
        double f = 0;
        if(n3 < 3*n2/2)f = K13 + K15 - K17;
        else if(n3 == 3*n2/2)f = K19;
        else if(n3 > 3*n2/2 && K20 > 0)f = -K13 + K15 + K16;
        else if(K20 == 0)f = K15 + K22;
        else if(K20 < 0)f = K13 + K15 + K17;
        else System.out.println("Error when computing f: undefined");
        
        System.out.println("f: "+f);
        return f;
    }

    /**
	 * Computes the minimum number of phases needed 
	 * to represent the tuple of the normalized moments
	 * @param n2 Second normalized moment
	 * @param n3 Third normalized moment
	 * @return Minimum number of phases needed to represent the 
	 * tuple of the normalized moments 
	 */
	public int getSize(double n2, double n3){
		int n = 0;
		if(n3<=n2 || n2<=1){
			System.out.println("The set of moments (n2, n3) " +
					"is not representable: ("+n2+","+n3+")");
		}else{
			
			if(n2>=2 && n3>=3  && n3==1.5*n2){
				System.out.println("The set of moments (n2, n3) is representable " +
						"by an exponential distribution: ("+n2+","+n3+")");
				return 1;
			}
				
			
			n=2;
			
			while(n <= 100){
				if(n2 >= 1.0+1.0/n){
					System.out.println("Second moment "+n2+" " +
							"representable with "+n+" phases");
					
					double un = 1/(n*n*n2)*( 2*(n-2)*(n*n2-n-1)*Math.sqrt(1+n*(n2-2)/(n-1)) + (n+2)*(3*n*n2-2*n-2) );
					double pn = (n+1)*(n2-2)/(3*n2*(n-1)) * (-2*Math.sqrt(n+1)/( Math.sqrt(4*(n+1)-3*n*n2) ) - 1);
					double an = (n2-2)/(pn*(1-n2) + Math.sqrt(pn*pn + pn*n*(n2-2)/(n-1))  );
					double ln = ((3+an)*(n-1)+2*an)/((n-1)*(1+an*pn)) - 2*an*(n+1)/(2*(n-1)+an*pn*(n*an+2*n-2));
					
					//lower bound
					if(n2 >= (n+1.0)/n && n2 < (n+4.0)/(n+1.0)){
						if(n3 >= ln){
							//upper bound
							if(n2 >= (n+1.0)/n && n2 <= n/(n-1.0)){
								if(n3<=un)return n;
							}else if(n2 > n/(n-1.0)){
								return n;
							}
						}
					}else if (n2 >= (n+4.0)/(n+1.0)){
						if(n3 > (n+1.0)*n2/n){
							//upper bound
							if(n2 >= (n+1.0)/n && n2 <= n/(n-1.0)){
								if(n3<=un)return n;
							}else if(n2 > n/(n-1.0)){
								return n;
							}
						}
					}
					System.out.println("Third moment "+n3+" " +
								"not representable with "+n+" phases");
					n++;						
				}else{
					n++;
				}
			}
		}
		
		return n;
	}
}
