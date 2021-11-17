DenseMatrix A = new DenseMatrix(
			new double[][] { {-4,2,1} , 
				{1,-3,1} , {2, 1,-5} } );
DenseVector alpha = new DenseVector(new double[] 
				{0.1, 0.2, 0.2});

DenseContPhaseVar v1 = new DenseContPhaseVar(alpha, A);

double rho = 0.5;
PhaseVar v2 = v1.waitingQ(rho);
System.out.println("v2:\n"+v2.toString());
		
		