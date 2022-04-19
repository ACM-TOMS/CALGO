ContPhaseVar v1 = DenseContPhaseVar.expo(3);
ContPhaseVar v2 = DenseContPhaseVar.Erlang(1.5, 2);
ContPhaseVar v3 = v1.max(v2);
System.out.println("P(v3<=2.0):\t" +v3.cdf(2.0));
		
		