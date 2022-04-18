DenseMatrix A = new DenseMatrix(new double[][] { {-4,2,1} , 
					{1,-3,1} , {2, 1,-5} } );
DenseVector alpha = new DenseVector(new double[] {0.3, 0.3, 0.4});
DenseContPhaseVar v1 = new DenseContPhaseVar(alpha, A);

NeutsContPHGenerator gen = new NeutsContPHGenerator(v1); 
double[] variates = new double[10];
variates = gen.getRandom(10);

for(int i = 0; i < 10; i++)
	System.out.println("var["+i+"]:"+variates[i]);
		

		
		