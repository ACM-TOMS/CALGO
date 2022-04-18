
package examples.jphase;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.StringTokenizer;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;

import jphase.ContPhaseVar;
import jphase.DenseContPhaseVar;
import jphase.PhaseVar;
import jphase.fit.EMHyperErlangFit;
import jphase.fit.MomentsACPHFit;
import jphase.generator.NeutsContPHGenerator;

/**
 * This class contains some basic examples that illustrate 
 * the use of jPhase 
 * @author Juan F. Pérez
 *
 */
public class BasicExamples {

	/**
	 * Main method to run the examples
	 * @param args
	 */
	public static void main(String[] args){
		BasicExamples myExample = new BasicExamples();
		myExample.example0();
		myExample.example1();
		myExample.example2();
		myExample.example3();
		myExample.example4();
		myExample.example5();
		myExample.example6();  
    }
	
	/**
	 * Example about closure properties and 
     * probabilities computation
	 */
	private void example0(){
        System.out.println("\nEXAMPLE 0");
        ContPhaseVar v1 = DenseContPhaseVar.expo(3);
		ContPhaseVar v2 = DenseContPhaseVar.Erlang(1.5, 2);
		ContPhaseVar v3 = v1.max(v2);
		System.out.println("P( v3 <= 2.0 ):\t" +v3.cdf(2.0));
	}
	
	
	/**
	 * Example about closure properties and output printing
	 */
	private void example1(){
        System.out.println("\nEXAMPLE 1");
		ContPhaseVar v1 = DenseContPhaseVar.Erlang(0.8, 3);
		ContPhaseVar v2 = DenseContPhaseVar.Erlang(1.5, 2);

		ContPhaseVar v3 = v1.sum(v2);
		System.out.println("v3:\n"+v3.toString());
	}
	
	/**
	 * Example about closure properties and construction
     * from double arrays
	 */
	private void example2(){
        System.out.println("\nEXAMPLE 2");
		DenseMatrix A = new DenseMatrix(new double[][] { {-4,2,1} , {1,-3,1} , {2, 1,-5} } );
		DenseVector alpha = new DenseVector(new double[] {0.1, 0.2, 0.2});
		DenseContPhaseVar v1 = new DenseContPhaseVar(alpha, A);
		double rho = 0.5;
		PhaseVar v2 = v1.waitingQ(rho);
		System.out.println("v2:\n"+v2.toString());
		
	}
	
	/**
	 * Example about fitting with maximum-likelihood
     * algorithms
	 */
	private void example3(){
        System.out.println("\nEXAMPLE 3");
		double[] data = readTextFile("examples/jphase/W2.txt");
        EMHyperErlangFit fitter = new EMHyperErlangFit(data); 
		ContPhaseVar v1 = fitter.fit(4);
		if(v1!=null){
            System.out.println("v1:\n"+v1.toString());
            System.out.println("logLH:\t"+fitter.getLogLikelihood());
        }
	}
	
	/**
	 * Example about fitting with moment-matching techniques
	 *
	 */
	private void example4(){
        System.out.println("\nEXAMPLE 4");
		MomentsACPHFit fitter = new MomentsACPHFit(2, 6, 25); 
		ContPhaseVar v1 = fitter.fit();
		System.out.println("v1:\n"+v1.toString());
		System.out.println("m1:\t"+v1.moment(1));
		System.out.println("m2:\t"+v1.moment(2));
		System.out.println("m3:\t"+v1.moment(3));
	}
	
    /**
     * Example about random variates generation and
     * construction from DenseMatrix and DenseVector
     */
	private void example5(){
        System.out.println("\nEXAMPLE 5");
		DenseMatrix A = new DenseMatrix(new double[][] { {-4,2,1} , {1,-3,1} , {2, 1,-5} } );
		DenseVector alpha = new DenseVector(new double[] {0.3, 0.3, 0.4});
		DenseContPhaseVar v1 = new DenseContPhaseVar(alpha, A);
		NeutsContPHGenerator gen = new NeutsContPHGenerator(v1); 
		double[] variates = new double[10];
		variates = gen.getRandom(10);
		for(int i = 0; i < 10; i++)
			System.out.println("var["+i+"]:\t"+variates[i]);
		
	}
	
	/**
     * Example about closure properties: residual time 
     * 
     */
    private void example6(){
        double[] lambdas = { .00030, 1, 2 };
        double[] alpha = { .50 ,.50, 0}; 
        ContPhaseVar v1 = DenseContPhaseVar.HyperExpo(lambdas,alpha);

        ContPhaseVar v2 = v1.copy();
        
        for (int i=0;i<1000;i++){
            System.out.println(v2.expectedValue());
            v2= v2.eqResidualTime();
        }
        System.out.println("v2:\n"+v2.toString());
    }
	
	/**
	 * This method reads a text file and stores the information
     * in a data array 
	 * @param fileName File name
	 * @return data array
	 */
	private double[] readTextFile(String fileName){
		ArrayList<Double> data = new ArrayList<Double>();
		
		try{
			FileReader archivo = new FileReader(fileName);
			BufferedReader entrada = new BufferedReader(archivo);
			String s;
			StringTokenizer str;
			System.out.println("Data file found");
											
			while (entrada.ready()){ 
				s=entrada.readLine();
				str=new StringTokenizer (s);
				while(str.countTokens() == 0){
					s=entrada.readLine();
					str=new StringTokenizer (s);
				}
				if (str.countTokens() != 1)throw new Exception ();
				data.add( new Double(Double.parseDouble(str.nextToken())) );
			}
			entrada.close();
			
		}catch(Exception e){
			System.out.println("File could not be read");
			return null;
		}
		double[] datos = new double[data.size()];
		for(int i = 0; i < data.size();i++)datos[i] = data.get(i).doubleValue();
		return datos;
	}
	
	
}
