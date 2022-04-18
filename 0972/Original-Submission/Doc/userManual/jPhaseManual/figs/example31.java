double[] data = readTextFile("data/W2.txt");

EMHyperErlangFit fitter = new EMHyperErlangFit(data); 

ContPhaseVar v1 = fitter.fit(4);

System.out.println("v1:\n"+v1.toString());

System.out.println("logLH:\t"+fitter.getLogLikelihood());

		
		